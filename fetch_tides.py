#!/usr/bin/env python3
"""
fetch_tides.py - the PUGET-TIDES data fetcher. Python standard library only.

It fetches and tidies the data. All the science is done in FORTRAN.

  python3 fetch_tides.py             daily run
  python3 fetch_tides.py --history   hourly water levels 1991-2020, for NORMALS

Writes (all in feet above mean lower low water, MLLW, times in UTC):
  harcon.csv          NOAA's 37 harmonic constituents for each station
  datums.csv          each station's tidal datums, 1983-2001 epoch
  water_level.csv     the last 30 days of 6-minute water levels. Not committed.
  pressure.csv        the last 30 days of air pressure where a station has a
                      barometer. Not committed.
  monthly_mean.csv    every monthly mean sea level on record. Not committed.
  noaa_hilo.csv       NOAA's own predicted highs and lows for the next week,
                      so FORTRAN can check itself against them. Not committed.
  noaa_trend.csv      NOAA's published sea-level trends, for comparison only.
  fetch_status.csv    what was fetched, what failed, and today's date.

Source: NOAA CO-OPS (Center for Operational Oceanographic Products and
Services), api.tidesandcurrents.noaa.gov. No key.

If something can't be fetched, the last good copy is kept and the failure is
recorded. Nothing is ever filled in with made-up values.
"""
import csv
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import date, datetime, timedelta, timezone

DG = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"
MD = "https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations"
DP = "https://api.tidesandcurrents.noaa.gov/dpapi/prod/webapi/product"
UA = {"User-Agent": "puget-tides (github.com/bdgroves/puget-tides)"}
DAYS = 30


def pacific_today(now_utc):
    """Today's date in Pacific time (DST from the second Sunday of March to
    the first Sunday of November), without needing tz data."""
    y = now_utc.year
    mar = date(y, 3, 8) + timedelta(days=(6 - date(y, 3, 8).weekday()) % 7)
    nov = date(y, 11, 1) + timedelta(days=(6 - date(y, 11, 1).weekday()) % 7)
    start = datetime(y, mar.month, mar.day, 10, tzinfo=timezone.utc)
    end = datetime(y, nov.month, nov.day, 9, tzinfo=timezone.utc)
    off = 7 if start <= now_utc < end else 8
    return (now_utc - timedelta(hours=off)).date()


def http(url, tries=3, timeout=60):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            last = e
            if e.code in (400, 404):
                raise                    # asking again won't help
            # 403 and 429 are NOAA asking us to slow down: back right off
            time.sleep((30 if e.code in (403, 429) else 4) * (i + 1))
        except Exception as e:          # noqa: BLE001
            last = e
            time.sleep(4 * (i + 1))
    raise last


def read_csv(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, cols, rows):
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)               # atomic: never a half-written file


def datagetter(station, product, begin, end, **extra):
    """Rows of a CO-OPS CSV product; [] when NOAA says there are no data."""
    q = dict(station=station, product=product, begin_date=begin, end_date=end, datum="MLLW",
             units="english", time_zone="gmt", format="csv", application="puget-tides")
    q.update(extra)
    text = http(f"{DG}?{urllib.parse.urlencode(q)}").decode("utf-8", "replace")
    if text.lstrip().startswith("Error") or "No data was found" in text[:300]:
        return []
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines or "Error" in lines[0]:
        raise ValueError(lines[0][:120] if lines else "empty response")
    head = [h.strip() for h in lines[0].split(",")]
    return [dict(zip(head, (v.strip() for v in ln.split(",")))) for ln in lines[1:]]


def ft(v):
    try:
        return f"{float(v):.3f}"
    except (TypeError, ValueError):
        return ""


def stamp(d):
    return d.strftime("%Y%m%d %H:%M")


def station_parts(st, now, status):
    """Everything for one station; each part separately, so one failure
    doesn't lose the rest."""
    sid, out = st["id"], {}
    begin = now - timedelta(days=DAYS)

    def part(name, fn):
        try:
            out[name] = fn()
        except Exception as e:           # noqa: BLE001
            out[name] = None
            status[f"{sid}_{name}"] = f"FAILED {str(e)[:70]}"

    def harcon():
        js = json.loads(http(f"{MD}/{sid}/harcon.json?units=english"))
        rows = [dict(station=sid, name=c["name"], amp_ft=c["amplitude"], phase_gmt=c["phase_GMT"],
                     speed=c["speed"]) for c in js["HarmonicConstituents"]]
        if len(rows) < 30:
            raise ValueError(f"only {len(rows)} constituents")
        return rows

    def datums():
        js = json.loads(http(f"{MD}/{sid}/datums.json?units=english"))
        rows = [dict(station=sid, name=d["name"], value_ft=d["value"]) for d in js["datums"]]
        for k in ("HAT", "LAT", "max", "min"):
            if js.get(k) is not None:
                d = js.get(k + "date", "") or ""
                when = f"{d[:4]}-{d[4:6]}-{d[6:8]} {js.get(k + 'time', '')}".strip() if len(d) == 8 else ""
                rows.append(dict(station=sid, name=k, value_ft=js[k], when=when))
        return rows

    def water():
        rows = datagetter(sid, "water_level", stamp(begin), stamp(now))
        return [dict(station=sid, time=r["Date Time"], ft=ft(r.get("Water Level")))
                for r in rows if ft(r.get("Water Level"))]

    def pressure():
        rows = datagetter(sid, "air_pressure", stamp(begin), stamp(now))
        return [dict(station=sid, time=r["Date Time"], hpa=r.get("Pressure", "")) for r in rows
                if r.get("Pressure", "").replace(".", "").isdigit()]

    def monthly():
        rows = datagetter(sid, "monthly_mean", "18900101", now.strftime("%Y%m%d"))   # one request covers it all
        return [dict(station=sid, year=r["Year"], month=r["Month"], msl=ft(r.get("MSL")),
                     highest=ft(r.get("Highest")), lowest=ft(r.get("Lowest"))) for r in rows]

    def hilo():
        a, b = now - timedelta(days=1), now + timedelta(days=8)
        rows = datagetter(sid, "predictions", a.strftime("%Y%m%d"), b.strftime("%Y%m%d"), interval="hilo")
        return [dict(station=sid, time=r["Date Time"], ft=ft(r.get("Prediction")), type=r.get("Type", ""))
                for r in rows]

    def trend():
        try:
            js = json.loads(http(f"{DP}/sealvltrends.json?station={sid}"))
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return []                # NOAA publishes no trend for this gauge
            raise
        if not js.get("SeaLvlTrends"):
            return []                    # NOAA publishes no trend for this gauge
        t = js["SeaLvlTrends"][0]
        return [dict(station=sid, mm_per_yr=f"{float(t['trend']) * 25.4 / 10:.2f}",
                     err_mm_per_yr=f"{float(t['trendError']) * 25.4 / 10:.2f}",
                     start=t.get("startDate", ""), end=t.get("endDate", ""))]

    for name, fn in (("harcon", harcon), ("datums", datums), ("water", water), ("pressure", pressure),
                     ("monthly", monthly), ("hilo", hilo), ("trend", trend)):
        part(name, fn)
    return sid, out


def merge(path, cols, fresh, sids_ok):
    """New rows for the stations that answered, the last good rows for the rest."""
    keep = [r for r in read_csv(path) if r.get("station") not in sids_ok]
    write_csv(path, cols, keep + fresh)


def daily():
    now = datetime.now(timezone.utc).replace(second=0, microsecond=0)
    today = pacific_today(now)
    print(f"PUGET-TIDES FETCH  {now:%Y-%m-%d %H:%M} UTC  (Pacific date {today})")
    stations = read_csv("stations.csv")
    status = {"run_utc": now.strftime("%Y-%m-%d %H:%M"), "pacific_date": str(today)}
    with ThreadPoolExecutor(4) as pool:
        results = dict(pool.map(lambda s: station_parts(s, now, status), stations))
    files = {
        "harcon": ("harcon.csv", ["station", "name", "amp_ft", "phase_gmt", "speed"]),
        "datums": ("datums.csv", ["station", "name", "value_ft", "when"]),
        "water": ("water_level.csv", ["station", "time", "ft"]),
        "pressure": ("pressure.csv", ["station", "time", "hpa"]),
        "monthly": ("monthly_mean.csv", ["station", "year", "month", "msl", "highest", "lowest"]),
        "hilo": ("noaa_hilo.csv", ["station", "time", "ft", "type"]),
        "trend": ("noaa_trend.csv", ["station", "mm_per_yr", "err_mm_per_yr", "start", "end"]),
    }
    failed = 0
    for part, (path, cols) in files.items():
        ok = {sid for sid, r in results.items() if r.get(part) is not None}
        fresh = [row for sid in ok for row in results[sid][part]]
        merge(path, cols, fresh, ok)
        if part in ("harcon", "datums", "water", "hilo"):
            failed += len(stations) - len(ok)
    for s in stations:
        r = results[s["id"]]
        n = len(r.get("water") or [])
        last = (r.get("water") or [{}])[-1].get("time", "none")
        bad = [k for k in files if r.get(k) is None]
        print(f"  {s['id']}  {s['name']:<14} {n:5d} water levels, latest {last}"
              + (f"; FAILED: {' '.join(bad)}" if bad else "")
              + ("" if r.get("pressure") else "; no barometer"))
        status[f"{s['id']}_water"] = str(n)
        status[f"{s['id']}_latest"] = last
    status["failures"] = str(failed)
    write_csv("fetch_status.csv", ["key", "value"], [dict(key=k, value=v) for k, v in status.items()])
    return 4 if failed else 0


def history(y0=1991, y1=2020):
    """Hourly water levels, one file per station and year, listed in files.txt."""
    os.makedirs("history", exist_ok=True)
    jobs = [(s["id"], y) for s in read_csv("stations.csv") for y in range(y0, y1 + 1)]

    def one(job):
        sid, y = job
        path = f"history/{sid}_{y}.csv"
        if os.path.exists(path):
            return path, "cached"
        time.sleep(0.5)
        try:
            rows = datagetter(sid, "hourly_height", f"{y}0101", f"{y}1231 23:00")
        except Exception as e:           # noqa: BLE001
            return None, f"{sid} {y} FAILED {str(e)[:60]}"
        rows = [dict(time=r["Date Time"], ft=ft(r.get("Water Level"))) for r in rows if ft(r.get("Water Level"))]
        if not rows:
            return None, f"{sid} {y} no data"
        write_csv(path, ["time", "ft"], rows)
        return path, f"{len(rows)} hours"

    with ThreadPoolExecutor(2) as pool:          # gently: NOAA throttles bursts
        done = list(pool.map(one, jobs))
    got = [p for p, _ in done if p]
    for p, note in done:
        if not p:
            print("  " + note)
    print(f"  {len(got)} station-years of hourly water levels")
    with open("history/files.txt", "w") as f:
        f.write("\n".join(os.path.basename(p) for p in got) + "\n")
    return 0 if got else 8


if __name__ == "__main__":
    if "--history" in sys.argv:
        print("PUGET-TIDES HISTORY  NOAA hourly water levels 1991-2020")
        sys.exit(history())
    sys.exit(daily())
