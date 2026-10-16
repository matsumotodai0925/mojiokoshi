#!/bin/bash
# GeoClaw を 1 ケース実行する。
#   使い方: scripts/04_run.sh CASE MAXLEVEL
#     CASE     = A | B | C | V2025
#     MAXLEVEL = 1 (粗格子 4 分のみ) | 3 (4 分→1 分→15 秒の入れ子)
#   MAXLEVEL=3 のときは、同じケースの粗格子結果(L1)があれば
#   その到達時間の 30 分前から 15 秒格子を強制する（計算時間短縮）。
set -euo pipefail
CASE=$1; MAXLEVEL=${2:-1}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
source "$ROOT/env/env.sh"
GEO=$ROOT/geoclaw
[ -x "$GEO/xgeoclaw" ] || (cd "$GEO" && make .exe > build.log 2>&1)
RUN=$ROOT/work/runs/${CASE}_L${MAXLEVEL}
rm -rf "$RUN"; mkdir -p "$RUN/_output"; cd "$RUN"
if [ "$MAXLEVEL" -ge 3 ] && [ -d "$ROOT/work/runs/${CASE}_L1/_output" ]; then
  python3 "$ROOT/scripts/04a_l3_start_times.py" "$CASE" "$RUN/t_l3.json"
  export T_L3_JSON=$RUN/t_l3.json
fi
CASE=$CASE MAXLEVEL=$MAXLEVEL python3 "$GEO/setrun.py"
mv *.data _output/
cd _output
echo "start $(date -Is) CASE=$CASE MAXLEVEL=$MAXLEVEL OMP=$OMP_NUM_THREADS" | tee run_info.txt
T0=$(date +%s)
"$GEO/xgeoclaw" > run.log 2>&1
echo "elapsed_s $(( $(date +%s) - T0 ))" | tee -a run_info.txt
echo "end $(date -Is)" | tee -a run_info.txt
