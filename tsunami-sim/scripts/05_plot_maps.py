"""
最大水位（沖合〜浅海域）と到達時間の地図を作る。
  results/<case>_L<lv>/map_wide.png       : 依頼範囲 135-180E, 30-65N
  results/<case>_L<lv>/map_pacific.png    : 計算領域全体（ヒロを含む）
  results/<case>_L<lv>/map_hokkaido_tohoku.png, map_hakodate.png : 入れ子格子(L3)のとき
  results/<case>_L<lv>/fgmax_*.npz        : 数値（再描画用）
使い方: python3 05_plot_maps.py CASE LEVEL
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import common as C
import viz

LEVELS_ETA = [0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 50]
LEVELS_LOCAL = [0.2, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 7, 10, 15]   # 函館・北海道東北の拡大図用


def cmap_eta(levels=LEVELS_ETA):
    cols = plt.get_cmap('YlOrRd')(np.linspace(0.08, 1.0, len(levels) - 1))
    return ListedColormap(cols), BoundaryNorm(levels, len(levels) - 1)


def plot_map(X, Y, B, eta, arr, case, title, fn, extent, t_levels, t_label_unit='h',
             gauges=True, fault=True, coast_lw=0.5, figsize=(9, 8), note='', levels=LEVELS_ETA):
    fig, ax = plt.subplots(figsize=figsize)
    cm, nm = cmap_eta(levels)
    land = np.where(B >= 0, 1.0, np.nan)
    ax.pcolormesh(X, Y, land, cmap=ListedColormap(['#d9d9d9']), shading='auto')
    pc = ax.pcolormesh(X, Y, np.ma.masked_less(eta, levels[0]), cmap=cm, norm=nm, shading='auto')
    ax.contour(X, Y, B, [0], colors='k', linewidths=coast_lw)
    Xa, Ya, ta, Ba = arr
    # 水深20m未満の沿岸は等値線が密集するため描かない
    tt = np.where(Ba < -20, ta, np.nan) / (3600.0 if t_label_unit == 'h' else 60.0)
    cs = ax.contour(Xa, Ya, tt, t_levels, colors='#1f4e9c', linewidths=0.8)
    ax.clabel(cs, fmt=('%g h' if t_label_unit == 'h' else '%g分'), fontsize=7)
    if fault:
        viz.draw_fault(ax, case)
    if gauges:
        for g in C.load_adjusted_gauges():
            if g['id'] in C.MAIN_GAUGES and extent[0] <= g['lon'] <= extent[1] and extent[2] <= g['lat'] <= extent[3]:
                ax.plot(g['lon'], g['lat'], 'o', mfc='cyan', mec='k', ms=5)
                ax.annotate(viz.GAUGE_JA[g['id']], (g['lon'], g['lat']), xytext=(4, 4),
                            textcoords='offset points', fontsize=8)
    ax.set_xlim(extent[:2]); ax.set_ylim(extent[2:])
    viz.set_lonlat_aspect(ax, np.mean(extent[2:]))
    xt = ax.get_xticks()
    ax.set_xticks(xt); ax.set_xticklabels([viz.lon_label(x, 1 if extent[1] - extent[0] < 4 else 0) for x in xt])
    ax.set_xlim(extent[:2])
    ax.set_ylabel('緯度 (°N)')
    cb = fig.colorbar(pc, ax=ax, shrink=0.75, extend='max', ticks=levels)
    cb.set_label('最大水位 [m]（静水面基準、海域のみ。遡上高ではない）')
    ax.set_title(title, fontsize=11)
    fig.text(0.01, 0.01, viz.disclaimer(case) + f'｜白: 最大水位{levels[0]:g} m未満の海域｜青線: 到達時間（直近10分で水位が' + ('5 cm' if case == 'V2025' else '0.2 m') + '以上変化した最初の時刻）' + note, fontsize=7.5)
    fig.savefig(fn, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('wrote', fn)


def main(case, level):
    rd = viz.resdir(case, level)
    lab = {'A': 'A: 1200×200 km', 'B': 'B: 1500×200 km', 'C': 'C: 2000×250 km',
           'V2025': '検証: 2025年7月 M8.8（USGS有限断層モデル）'}[case]
    grid = {1: '4分格子のみ（粗）', 3: '4分→1分→15秒 入れ子'}[level]
    arr1 = viz.arrival_from_fgout(case, level, 1)
    viz.save_npz(os.path.join(rd, 'arrival_wide.npz'), X=arr1[0], Y=arr1[1], t=arr1[2], B=arr1[3])
    X, Y, B, eta, t = viz.read_fgmax(case, level, 1)
    viz.save_npz(os.path.join(rd, 'fgmax_wide.npz'), X=X, Y=Y, B=B, eta=eta)
    w = C.MAP_WIDE
    plot_map(X, Y, B, eta, arr1, case, f'最大水位と到達時間 — {lab}（{grid}）',
             os.path.join(rd, 'map_wide.png'), [w['x1'], w['x2'], w['y1'], w['y2']],
             np.arange(0.5, 12.1, 0.5), note='（0.5時間毎）')
    d = C.DOMAIN
    plot_map(X, Y, B, eta, arr1, case, f'最大水位と到達時間（計算領域全体）— {lab}',
             os.path.join(rd, 'map_pacific.png'), [d['x1'], d['x2'], d['y1'], d['y2']],
             np.arange(1, 13, 1), figsize=(12, 7.5), coast_lw=0.3, note='（1時間毎）')
    if level >= 3:
        X, Y, B, eta, t = viz.read_fgmax(case, level, 2)
        viz.save_npz(os.path.join(rd, 'fgmax_hokkaido_tohoku.npz'), X=X, Y=Y, B=B, eta=eta)
        b = C.L2_BOXES['hokkaido_tohoku']
        plot_map(X, Y, B, eta, arr1, case, f'北海道・東北太平洋側（1分格子）— {lab}',
                 os.path.join(rd, 'map_hokkaido_tohoku.png'), b, np.arange(0, 241, 15),
                 t_label_unit='min', note='（15分毎, 到達時間は6分格子）。沿岸の値は格子が粗く信頼度が低い', levels=LEVELS_LOCAL)
        arr2 = viz.arrival_from_fgout(case, level, 2)
        viz.save_npz(os.path.join(rd, 'arrival_hakodate.npz'), X=arr2[0], Y=arr2[1], t=arr2[2], B=arr2[3])
        X, Y, B, eta, t = viz.read_fgmax(case, level, 3)
        viz.save_npz(os.path.join(rd, 'fgmax_hakodate.npz'), X=X, Y=Y, B=B, eta=eta)
        b = C.L3_BOXES['hakodate']
        t0 = np.nanmin(arr2[2]) / 60 if np.isfinite(arr2[2]).any() else 0
        plot_map(X, Y, B, eta, arr2, case, f'函館周辺（津軽海峡〜函館湾, 15秒格子）— {lab}',
                 os.path.join(rd, 'map_hakodate.png'), b,
                 np.arange(np.floor(t0 / 10) * 10, 721, 10), t_label_unit='min', fault=False,
                 figsize=(10, 6.5), note='（10分毎）。GEBCO 15秒（約460m）格子：港湾内・沿岸の値は参考程度', levels=LEVELS_LOCAL)


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]))
