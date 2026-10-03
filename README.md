```
╔══════════════════════════════════════════════════════════════════════════╗
║                                                                          ║
║   ██████╗ ██╗   ██╗ ██████╗ ███████╗████████╗                            ║
║   ██╔══██╗██║   ██║██╔════╝ ██╔════╝╚══██╔══╝                            ║
║   ██████╔╝██║   ██║██║  ███╗█████╗     ██║     T I D E S                 ║
║   ██╔═══╝ ██║   ██║██║   ██║██╔══╝     ██║                               ║
║   ██║     ╚██████╔╝╚██████╔╝███████╗   ██║                               ║
║   ╚═╝      ╚═════╝  ╚═════╝ ╚══════╝   ╚═╝                               ║
║                                                                          ║
║   TIDES OF PUGET SOUND  ·  FORTRAN  ·  NOAA CO-OPS                       ║
║   NEAH BAY → PORT ANGELES → FRIDAY HARBOR → CHERRY POINT →               ║
║   PORT TOWNSEND → BREMERTON → SEATTLE → TACOMA                           ║
║                                                                          ║
╚══════════════════════════════════════════════════════════════════════════╝
```

**[→ The live page: brooksgroves.com/puget-tides](https://brooksgroves.com/puget-tides/)**

Eight NOAA tide gauges, from Neah Bay on the open Pacific to Tacoma at the far end of Puget Sound, run through a FORTRAN batch job every six hours. **PUGET-TIDES.f90** predicts the tide at each gauge from NOAA's 37 harmonic constituents, the way tide tables are made, and checks its highs and lows against NOAA's published predictions every run. Then it measures the surge (water level minus tide), ranks it against 1991–2020, finds the king tides and the daylight minus tides of the coming year, and fits the sea-level trend over each gauge's whole record. It prints a 132-column report.

It is the third in a family with [SIERRA-FLOW](https://github.com/bdgroves/sierra-flow-cobol) (Sierra rivers in COBOL) and [CASCADIA-WX](https://github.com/bdgroves/cascadia-wx) (Northwest weather balloons in FORTRAN).

---

## The batch job

At 00:23, 06:23, 12:23 and 18:23 UTC, GitHub Actions runs `run_job.sh` on Ubuntu 24.04:

```
STEP FETCH     python3 fetch_tides.py      NOAA CO-OPS → harcon, datums, 30 days of water levels, monthly means
STEP COMPILE   gfortran -O2                TIDE_PHYS.f90 + PUGET-TIDES.f90, GNU Fortran 13
STEP PUGETIDE  ./puget-tides               → puget-tides-report.txt and the CSVs below
```

When `normals.csv` is missing, or the workflow is run by hand with *rebuild normals* ticked, two more steps run first:

```
STEP HISTORY   python3 fetch_tides.py --history   hourly water levels 1991-2020, every gauge → history/
STEP NORMALS   ./normals                          TIDE_PHYS.f90 + NORMALS.f90 → normals.csv, normals-log.txt
```

| RC | Meaning |
|---|---|
| 0 | Normal: every gauge reported within 3 hours and FORTRAN matches NOAA |
| 4 | Warning: a gauge late or missing (`TIDE015W`, `TIDE016W`), no normals (`TIDE018W`), or FORTRAN and NOAA differ by more than 3 minutes or 0.05 ft (`TIDE020W`) |
| 8 | No harmonic constituents at all. Nothing is written |
| 12 | An input file can't be opened |

If NOAA doesn't answer for a gauge, the fetcher keeps that gauge's last good data and the FORTRAN flags it. Nothing is ever filled in.

## What it computes

**The tide.** h(t) = Z₀ + Σ f·H·cos(V + u − κ), summed over 37 constituents (Schureman 1958). H and κ are NOAA's amplitude and Greenwich phase for each constituent at each gauge; Z₀ is mean sea level above mean lower low water from the 1983–2001 datum epoch. FORTRAN computes everything else: the mean longitudes of the moon, sun, lunar perigee, lunar node and solar perigee (Meeus 1998), the equilibrium arguments V from Doodson numbers, and the 18.6-year nodal corrections f and u, taken at the middle of each year as NOAA does. Highs and lows are found by scanning the analytic rate of change every 6 minutes and bisecting each turn to about a second.

**Checked.** Over the whole of 2020 and 2026, FORTRAN's hourly tide matches NOAA's published hourly predictions to an RMS of 0.002 ft at Seattle and Neah Bay, at most 0.011 ft. On the first full run, all 257 of NOAA's highs and lows for the coming week matched FORTRAN's at all eight gauges, within 0.8 minutes and 0.010 ft. Every run repeats the check and warns if it slips.

One finding on the way: NOAA's M1 (a small lunar diurnal constituent, about an inch at Seattle) starts each year from T − s + h − 90° without the lunar perigee p, but runs at a speed that includes p's motion. Matching NOAA means subtracting p at January 1 of each year from M1's argument; without that, FORTRAN was off by up to 0.3 ft in 2020.

**Surge.** Measured water level minus FORTRAN's tide, for every 6-minute reading of the last 30 days: the latest, the 24-hour mean, and the 30-day highest, lowest and mean. Where a gauge has a barometer, the inverse-barometer effect, (1013.25 − p) × 0.0328 ft per hPa, shows how much of the surge air pressure explains. Surge also holds any change in mean sea level since the 1983–2001 epoch, and anything the constituents miss.

**Normals.** NORMALS.f90 reads every hourly water level from 1991 through 2020, subtracts the tide it predicts for that hour, and keeps each Pacific-standard-time day's mean surge and highest water (days with at least 20 hours). For each day of the year it keeps the 10th, 25th, 50th, 75th and 90th percentiles of surge from every day within 15 days of the date, and of highest water within 7 days (Weibull plotting positions, a 366-day calendar). A gauge needs 20 years. Classes follow the USGS's WaterWatch: below the 10th percentile is much below normal, 25th–75th normal, above the 90th much above.

**King tides.** The ten highest predicted highs of the next 365 days, no two within three days of each other, so each is a different spell.

**Daylight minus tides.** Every predicted low below 0 ft MLLW in the next year that falls while the sun is above −0.833° (sunrise and sunset), with the sun's position from NOAA's solar equations, computed here.

**The tide comes in.** The phase of the main lunar tide, M2, at each gauge against Neah Bay's, divided by its speed (28.984° per hour), is how many hours later it arrives; the ratio of amplitudes is how much it has grown or shrunk.

**Sea level.** Least squares on every monthly mean on record, with a separate mean for each calendar month so the seasons don't bias the slope. The 95% interval is widened for the months' lag-one autocorrelation; NOAA's model gives a narrower interval, but the trends agree to within 0.05 mm/yr at every gauge NOAA publishes. A gauge needs 240 monthly means.

## The eight gauges

| Gauge | NOAA ID | Where | Notes |
|---|---|---|---|
| Neah Bay | 9443090 | Cape Flattery | The ocean tide. Sea level here is falling: the land is rising faster |
| Port Angeles | 9444090 | Strait of Juan de Fuca | The smallest main tide of the eight |
| Friday Harbor | 9449880 | San Juan Islands | |
| Cherry Point | 9449424 | Strait of Georgia | |
| Port Townsend | 9444900 | Admiralty Inlet | The door to Puget Sound |
| Bremerton | 9445958 | Sinclair Inlet | No record between 1978 and 2021, so no normals or trend |
| Seattle | 9447130 | Elliott Bay | Monthly means since 1898 |
| Tacoma | 9446484 | Commencement Bay | Since 1996; the biggest tide of the eight |

## Files

| File | What |
|---|---|
| `TIDE_PHYS.f90` | Shared module: dates, Pacific time, the tide from its constituents, highs and lows, the sun, the trend. Kept to 72 columns |
| `PUGET-TIDES.f90` | The batch program |
| `NORMALS.f90` | Builds `normals.csv` from 30 years of hourly water levels |
| `fetch_tides.py` | The fetcher. Python standard library only |
| `stations.csv` | The eight gauges, ocean to inland |
| `harcon.csv`, `datums.csv` | NOAA's constituents and datums, refreshed every run |
| `now.csv`, `series.csv`, `hilo.csv` | The water now, the week, the highs and lows with NOAA's beside them |
| `kingtides.csv`, `minustides.csv` | The year ahead |
| `travel.csv`, `trend.csv`, `annual.csv` | The tide's travel, the sea-level trend, annual means |
| `summary.csv`, `puget-tides-report.txt`, `job-log.txt` | Summary, the printed report, the job log |
| `index.html` | The page. It reads the files above and draws them; it calculates nothing |

## Run it yourself

```bash
sudo apt install gfortran        # Ubuntu / WSL;  macOS: brew install gcc
./run_job.sh                     # the first run also builds normals.csv, a few minutes
cat puget-tides-report.txt
python3 -m http.server           # then open http://localhost:8000
```

On Windows, run it in WSL.

## Sources

- NOAA CO-OPS Data API and Metadata API: water levels, harmonic constituents, datums, monthly means, predictions, sea-level trends. api.tidesandcurrents.noaa.gov, no key.
- Schureman, P. (1958). *Manual of Harmonic Analysis and Prediction of Tides.* U.S. Coast and Geodetic Survey Special Publication 98.
- Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed.
- Zervas, C. (2009). *Sea Level Variations of the United States 1854–2006.* NOAA Technical Report NOS CO-OPS 053.

Predictions are for the gauges, not for navigation.

---

```
  PUGET-TIDES  V1.0
  NORMAL TERMINATION.  RETURN CODE 0.
```
