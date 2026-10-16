"""
A/B/C の比較（幅で示す）。
  results/compare_maxmap_L<lv>.png : 3 ケースの最大水位図を並べたもの（依頼範囲 135-180E, 30-65N）
  results/compare_range_L<lv>.md   : 主要地点の到達時間・最大水位の範囲（最小〜最大）
使い方: python3 09_compare.py LEVEL
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import common as C
import viz

sys.path.insert(0, os.path.dirname(__file__))
pm = __import__('05_plot_maps')


def main(level):
    fig, axs = plt.subplots(1, 3, figsize=(18, 5.6))
    cm, nm = pm.cmap_eta()
    w = C.MAP_WIDE
    for ax, case in zip(axs, ['A', 'B', 'C']):
        rd = os.path.join(C.RESULTS, f'{case}_L{level}')
        d = np.load(os.path.join(rd, 'fgmax_wide.npz'))
        a = np.load(os.path.join(rd, 'arrival_wide.npz'))
        X, Y, B, eta = d['X'], d['Y'], d['B'], d['eta']
        ax.pcolormesh(X, Y, np.where(B >= 0, 1.0, np.nan), cmap=ListedColormap(['#d9d9d9']), shading='auto')
        pc = ax.pcolormesh(X, Y, np.ma.masked_less(eta, pm.LEVELS_ETA[0]), cmap=cm, norm=nm, shading='auto')
        ax.contour(X, Y, B, [0], colors='k', linewidths=0.4)
        tt = np.where(a['B'] < -20, a['t'], np.nan) / 3600
        cs = ax.contour(a['X'], a['Y'], tt, np.arange(1, 13, 1), colors='#1f4e9c', linewidths=0.7)
        ax.clabel(cs, fmt='%g h', fontsize=7)
        viz.draw_fault(ax, case)
        ax.set_xlim(w['x1'], w['x2']); ax.set_ylim(w['y1'], w['y2'])
        viz.set_lonlat_aspect(ax, 47.5)
        sc = C.SCENARIOS[case]
        slip = C.mean_slip(sc['L'], sc['W'])
        ax.set_title(f"ケース{case}: {sc['L']/1e3:.0f}×{sc['W']/1e3:.0f} km, 平均すべり {slip:.1f} m"
                     + ('（60m超・極端）' if slip > C.SLIP_WARN else ''), fontsize=10)
    cb = fig.colorbar(pc, ax=axs, shrink=0.8, ticks=pm.LEVELS_ETA, extend='max')
    cb.set_label('最大水位 [m]（海域のみ・遡上高ではない）')
    grid = {1: '4分格子', 3: '4分→1分→15秒 入れ子'}[level]
    fig.suptitle(f'M9.9 仮想シナリオ A/B/C の比較（最大水位と到達時間[1時間毎], {grid}）\n{viz.DISCLAIMER}', fontsize=11)
    fn = os.path.join(C.RESULTS, f'compare_maxmap_L{level}.png')
    fig.savefig(fn, dpi=130, bbox_inches='tight'); plt.close(fig); print('wrote', fn)

    # 地点別の幅
    rows = {}
    for line in open(os.path.join(C.RESULTS, f'gauge_summary_L{level}.csv')).readlines()[1:]:
        p = line.strip().split(',')
        rows.setdefault(int(p[1]), {})[p[0]] = dict(t=float(p[6]), mx=float(p[8]), tmax=float(p[9]))

    def hm(s):
        return '—' if not np.isfinite(s) else f'{int(s // 3600)}:{int(s % 3600 // 60):02d}'
    fn = os.path.join(C.RESULTS, f'compare_range_L{level}.md')
    with open(fn, 'w') as f:
        f.write('| 地点 | 到達 A / B / C | 到達の幅 | 最大水位 A / B / C [m] | 最大水位の幅 [m] |\n|---|---|---|---|---|\n')
        for gid in C.MAIN_GAUGES:
            r = rows.get(gid, {})
            ts = [r.get(c, {}).get('t', np.nan) for c in 'ABC']
            ms = [r.get(c, {}).get('mx', np.nan) for c in 'ABC']
            ok_t = [x for x in ts if np.isfinite(x)]
            ok_m = [x for x in ms if np.isfinite(x)]
            tr = f'{hm(min(ok_t))}〜{hm(max(ok_t))}' if ok_t else '計算できず'
            mr = f'{min(ok_m):.1f}〜{max(ok_m):.1f}' if ok_m else '計算できず'
            f.write(f"| {viz.GAUGE_JA[gid]} | {' / '.join(hm(x) for x in ts)} | {tr} | "
                    f"{' / '.join('—' if not np.isfinite(x) else f'{x:.1f}' for x in ms)} | {mr} |\n")
    print(open(fn).read())


if __name__ == '__main__':
    main(int(sys.argv[1]))
