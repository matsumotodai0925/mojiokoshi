"""
津波伝播アニメーション（mp4）。fgout 出力を使う。
  results/<case>_L<lv>/animation_pacific.mp4   : 計算領域全体（6分格子・5分毎）
  results/<case>_L<lv>/animation_hakodate.mp4  : 函館周辺（15秒格子・2分毎, L3 のみ）
使い方: python3 07_animation.py CASE LEVEL
"""
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.colors import SymLogNorm, ListedColormap
from clawpack.geoclaw import fgout_tools
import common as C
import viz


def make(case, level, fgno, fname, vmax, figsize, title):
    out = viz.outdir(case, level)
    frames = sorted(glob.glob(os.path.join(out, f'fgout{fgno:04d}.t*')))
    if not frames:
        print('no fgout', fgno); return
    fg = fgout_tools.FGoutGrid(fgno, out, 'binary32')
    fg.read_fgout_grids_data(data_file=os.path.join(out, 'fgout_grids.data'))
    fr = fg.read_frame(1)
    X, Y = fr.X, fr.Y
    fl = fg.read_frame(len(frames))       # 陸の表示は最終時刻（細格子が有効）の地形を使う
    B0 = fl.eta - fl.h
    fig, ax = plt.subplots(figsize=figsize)
    ax.pcolormesh(X, Y, np.where(B0 >= 0, 1, np.nan), cmap=ListedColormap(['#cfcfcf']), shading='auto')
    norm = SymLogNorm(linthresh=vmax / 50, vmin=-vmax, vmax=vmax, base=10)
    eta = np.ma.masked_where((fr.h < 0.01) | (B0 >= 0), fr.eta)
    pc = ax.pcolormesh(X, Y, eta, cmap='RdBu_r', norm=norm, shading='auto')
    ax.contour(X, Y, B0, [0], colors='k', linewidths=0.4)
    viz.draw_fault(ax, case, lw=0.8)
    for g in C.load_adjusted_gauges():
        if g['id'] in C.MAIN_GAUGES and X.min() <= g['lon'] <= X.max() and Y.min() <= g['lat'] <= Y.max():
            ax.plot(g['lon'], g['lat'], 'o', mfc='yellow', mec='k', ms=4)
    viz.set_lonlat_aspect(ax, Y.mean())
    xt = ax.get_xticks(); ax.set_xticks(xt); ax.set_xticklabels([viz.lon_label(x, 1 if X.max() - X.min() < 4 else 0) for x in xt])
    ax.set_xlim(X.min(), X.max()); ax.set_ylim(Y.min(), Y.max())
    cb = fig.colorbar(pc, ax=ax, shrink=0.7, extend='both')
    cb.set_label('水位 [m]（対数目盛）')
    ttl = ax.set_title('')
    fig.text(0.01, 0.01, viz.disclaimer(case), fontsize=8)

    def update(k):
        f = fg.read_frame(k)
        e = np.ma.masked_where((f.h < 0.01) | (B0 >= 0), f.eta)
        pc.set_array(e.ravel())
        ttl.set_text(f'{title}  地震発生から {f.t / 3600:5.2f} 時間')
        return pc, ttl

    anim = animation.FuncAnimation(fig, update, frames=range(1, len(frames) + 1), blit=False)
    writer = animation.FFMpegWriter(fps=15, bitrate=1200)
    anim.save(fname, writer=writer, dpi=110)
    plt.close(fig)
    print('wrote', fname, len(frames), 'frames')


def main(case, level):
    rd = viz.resdir(case, level)
    lab = {'A': 'ケースA', 'B': 'ケースB', 'C': 'ケースC', 'V2025': '検証2025 M8.8'}[case]
    vmax = 1.0 if case == 'V2025' else 10.0
    make(case, level, 1, os.path.join(rd, 'animation_pacific.mp4'), vmax, (12, 7.5), lab)
    if level >= 3:
        make(case, level, 2, os.path.join(rd, 'animation_hakodate.mp4'), vmax / 2, (10, 6.5), lab + ' 函館周辺')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]))
