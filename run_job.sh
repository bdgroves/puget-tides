#!/usr/bin/env bash
# The PUGET-TIDES batch job: fetch, compile, (normals), run. Writes job-log.txt.
#   ./run_job.sh              regular run
#   ./run_job.sh --normals    also rebuild normals.csv from 30 years of water levels
# Every step's output and return code goes in the log. The job carries on
# through warnings (RC 4) and stops on RC >= 8, like a mainframe job step
# with COND=(8,LE).
set -u
export TZ=UTC
LOG=job-log.txt
REBUILD=0; [ "${1:-}" = "--normals" ] && REBUILD=1; [ -f normals.csv ] || REBUILD=1
START=$(date +%s.%N)
MAXRC=0
fcv=$(gfortran --version 2>/dev/null | head -1)
{
  echo "JOB PUGETIDE   $(date '+%Y-%m-%d %H:%M:%S') UTC"
  echo "HOST           ${RUNNER_OS:-$(uname -s)} $(uname -m) · ${fcv:-gfortran not found} · $(python3 --version 2>&1)"
  [ -n "${GITHUB_RUN_NUMBER:-}" ] && echo "RUN            GitHub Actions run ${GITHUB_RUN_NUMBER} (${GITHUB_EVENT_NAME:-})"
  echo
} > "$LOG"

step () {   # step NAME command...
  local name=$1; shift
  local t0=$(date +%s.%N)
  local out; out=$("$@" 2>&1); local rc=$?
  local t=$(echo "$(date +%s.%N) - $t0" | bc)
  printf "STEP %-10s RC=%04d  %6.2fs  %s\n" "$name" "$rc" "$t" "$*" >> "$LOG"
  echo "$out" | grep -v -e '^$' | sed 's/^/    /' >> "$LOG"
  echo >> "$LOG"
  echo "$out"
  [ $rc -gt $MAXRC ] && MAXRC=$rc
  return $rc
}

step FETCH timeout 900 python3 -u fetch_tides.py; rc=$?; [ $rc -ge 8 ] && exit 1
step COMPILE gfortran -O2 -o puget-tides TIDE_PHYS.f90 PUGET-TIDES.f90 || exit 1
if [ $REBUILD = 1 ]; then
  step HISTORY timeout 2400 python3 -u fetch_tides.py --history || exit 1
  step COMPILE gfortran -O2 -o normals TIDE_PHYS.f90 NORMALS.f90 || exit 1
  step NORMALS timeout 1800 ./normals; rc=$?; [ $rc -ge 8 ] && exit 1
  { echo "NORMALS BUILT  $(date '+%Y-%m-%d %H:%M:%S') UTC  from NOAA hourly water levels, 1991-2020"; echo;
    sed -n '/STEP HISTORY/,$p' "$LOG" | grep -v -e '^STEP COMPILE'; } > normals-log.txt
  rm -rf history
else
  echo "STEP NORMALS    skipped: normals.csv is current (see normals-log.txt)" >> "$LOG"
  echo >> "$LOG"
fi
step PUGETIDE timeout 300 ./puget-tides; rc=$?; [ $rc -ge 8 ] && exit 1
T=$(echo "$(date +%s.%N) - $START" | bc)
printf "JOB PUGETIDE   ENDED  MAXCC=%04d  %.1fs\n" "$MAXRC" "$T" >> "$LOG"
rm -f *.mod
exit 0
