"""
地点別の水位時系列（A/B/C を重ねて比較）と、到達時間・最大水位の表。
  results/gauges_L<lv>.png
  results/gauge_summary_L<lv>.csv / .md
  results/<case>_L<lv>/gauges_timeseries.csv
使い方: python3 06_gauges.py LEVEL [CASES...]
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import common as C
import viz

COL = {'A': '#d62728', 'B': '#1f77b4', 'C': '#2ca02c', 'V2025': 'k'}


def main(level, cases):
    gauges = {g['id']: g for g in C.load_adjusted_gauges()}
    fig, axs = plt.subplots(len(C.MAIN_GAUGES), 1, figsize=(11, 2.3 * len(C.MAIN_GAUGES)), sharex=True)
    rows = []
    for case in cases:
        if not os.path.isdir(viz.outdir(case, level)):
            print('skip', case); continue
        ts = {}
        for ax, gid in zip(axs, C.MAIN_GAUGES):
            t, eta, h = viz.read_gauge(case, level, gid)
            st = viz.gauge_stats(t, eta, h, tol=viz.arr_tol(case))
            g = gauges[gid]
            lo, la = viz.gauge_location(case, level, gid)
            rows.append(dict(case=case, gid=gid, name=viz.GAUGE_JA[gid], lon=lo, lat=la,
                             depth=g['depth'] if abs(lo - g['lon']) + abs(la - g['lat']) < 1e-3 else np.nan, **st))
            ax.plot(t / 3600, np.where(h > 0.01, eta, np.nan), color=COL[case], lw=0.9, label=f'ケース{case}')
            ts[gid] = (t, eta)
        # 時系列 CSV（15 秒間隔、地点ごとに時刻がわずかに異なるので補間して揃える）
        tt = np.arange(0, C.TFINAL + 1, 15.0)
        with open(os.path.join(viz.resdir(case, level), 'gauges_timeseries.csv'), 'w') as f:
            f.write('t_s,' + ','.join(gauges[g]['name'] + '_eta_m' for g in C.MAIN_GAUGES) + '\n')
            cols = [np.interp(tt, *ts[g]) for g in C.MAIN_GAUGES]
            for i, x in enumerate(tt):
                f.write(f'{x:.0f},' + ','.join(f'{c[i]:.3f}' for c in cols) + '\n')
    for ax, gid in zip(axs, C.MAIN_GAUGES):
        g = dict(gauges[gid])
        lo, la = viz.gauge_location(cases[0], level, gid)
        if abs(lo - g['lon']) + abs(la - g['lat']) > 1e-3:
            g.update(lon=lo, lat=la, depth=np.nan)    # 旧設定（水深10m以深）で計算した結果
        ax.axhline(0, color='0.6', lw=0.5)
        ax.set_ylabel('水位 [m]')
        lon = f"{g['lon']:.3f}E" if g['lon'] <= 180 else f"{360 - g['lon']:.3f}W"
        dep = f"GEBCO水深 {g['depth']:.0f} m" if np.isfinite(g['depth']) else '水深10m以深の点'
        ax.set_title(f"{viz.GAUGE_JA[gid]}（計算点 {lon}, {g['lat']:.3f}N, {dep}）",
                     fontsize=9, loc='left')
        ax.grid(alpha=0.3)
    axs[0].legend(ncol=4, fontsize=8, loc='upper right')
    axs[-1].set_xlabel('地震発生からの時間 [h]')
    axs[-1].set_xlim(0, C.TFINAL / 3600)
    grid = {1: '4分格子のみ（粗格子、沿岸点は信頼度低い）', 3: '4分→1分→15秒 入れ子'}[level]
    fig.suptitle(f'地点別の水位時系列（{grid}）\n{viz.DISCLAIMER}', fontsize=11)
    fig.tight_layout()
    fn = os.path.join(C.RESULTS, f'gauges_L{level}.png')
    fig.savefig(fn, dpi=140); plt.close(fig); print('wrote', fn)

    fn = os.path.join(C.RESULTS, f'gauge_summary_L{level}.csv')
    keys = ['case', 'gid', 'name', 'lon', 'lat', 'depth', 't_arrival', 't_arrival_pos', 'eta_max', 't_max', 'eta_min', 'eta0', 'dry_frac']
    with open(fn, 'w') as f:
        f.write(','.join(keys) + '\n')
        for r in rows:
            f.write(','.join(f'{r[k]:.4g}' if isinstance(r[k], (float, np.floating)) else str(r[k]) for k in keys) + '\n')
    fn = os.path.join(C.RESULTS, f'gauge_summary_L{level}.md')

    def hm(s):
        return '—' if not np.isfinite(s) else f'{int(s // 3600)}:{int(s % 3600 // 60):02d}'
    with open(fn, 'w') as f:
        f.write('| 地点 | ケース | 到達 (h:mm) | 最大水位 [m]（平均海面基準） | 最大の時刻 | 最低水位 [m] | 地殻変動による初期水位変化 [m] | 最大水位−初期変化 [m]（変動後の海底・地盤基準の目安） |\n|---|---|---|---|---|---|---|---|\n')
        for gid in C.MAIN_GAUGES:
            for r in [r for r in rows if r['gid'] == gid]:
                if not np.isfinite(r['eta_max']):
                    f.write(f"| {r['name']} | {r['case']} | 計算できず（この格子では計算点が陸） | — | — | — | — | — |\n")
                    continue
                f.write(f"| {r['name']} | {r['case']} | {hm(r['t_arrival'])} | {r['eta_max']:.2f} | {hm(r['t_max'])} | "
                        f"{r['eta_min']:.2f}{'（引き波で一時干出）' if r['dried'] else ''} | {r['eta0']:+.2f} | {r['eta_max'] - r['eta0']:.2f} |\n")
        f.write('\n到達: 地殻変動完了後、直近10分間で水位が 0.2 m 以上変化した最初の時刻（有意な津波の到達。地図の到達時間と同じ定義）。それより前に数cm〜15cm程度の長周期の先行変動が現れる場合がある。初期水位変化は地震5分後の水位。'
                '水位は静水面(平均海面)基準、潮汐は考慮しない。沖合〜浅海域の計算点の値で、遡上高ではない。「一時干出」は引き波で計算点の水深がなくなったことを示し、その最低水位は計算セルの海底高さ。\n')
    print('wrote', fn)


if __name__ == '__main__':
    main(int(sys.argv[1]), sys.argv[2:] or ['A', 'B', 'C'])
