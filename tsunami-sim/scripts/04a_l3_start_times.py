"""粗格子(L1)結果の到達時間から、15 秒格子を強制し始める時刻(到達 30 分前)を JSON に書く。
使い方: python3 04a_l3_start_times.py CASE OUT.json"""
import os
import sys
import json
import numpy as np
from clawpack.geoclaw import fgmax_tools
import common as C

case = sys.argv[1]
out = os.path.join(C.RUNS, f'{case}_L1', '_output')
fg = fgmax_tools.FGmaxGrid()
fg.read_fgmax_grids_data(1, os.path.join(out, 'fgmax_grids.data'))
fg.read_output(outdir=out, verbose=False)
res = {}
for name, (x1, x2, y1, y2) in C.L3_BOXES.items():
    # 箱の周囲 1 度まで含めて最も早い到達を採る（粗格子では湾内が陸になるため）
    m = (fg.X >= x1 - 1) & (fg.X <= x2 + 1) & (fg.Y >= y1 - 1) & (fg.Y <= y2 + 1)
    t = np.ma.masked_invalid(np.ma.array(fg.arrival_time)[m])
    res[name] = max(0.0, float(t.min()) - 1800.0) if t.count() else 0.0
with open(sys.argv[2], 'w') as f:
    json.dump(res, f, indent=1)
print(json.dumps(res, indent=1))
