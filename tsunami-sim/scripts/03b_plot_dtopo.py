"""断層位置と初期海底上下変位（Okada）の図: results/fault_dtopo.png"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from clawpack.geoclaw import dtopotools
import common as C
import viz

cases = ['A', 'B', 'C', 'V2025']
fig, axs = plt.subplots(1, 4, figsize=(18, 6.2))
tr = np.loadtxt(os.path.join(C.DATA, 'trench_trace.csv'), delimiter=',', skiprows=1)
from clawpack.geoclaw import topotools
_t = topotools.Topography(os.path.join(C.TOPO_DIR, 'topo_4min.tt3'), topo_type=3)
topo = dict(X=_t.X, Y=_t.Y, B=_t.Z)
for ax, case in zip(axs, cases):
    d = dtopotools.DTopography(os.path.join(C.DTOPO_DIR, f'dtopo_{case}.tt3'), dtopo_type=3)
    dz = d.dZ[-1]
    lim = max(abs(dz.min()), dz.max())
    pc = ax.pcolormesh(d.X, d.Y, dz, cmap='RdBu_r', norm=TwoSlopeNorm(0, -lim, lim), shading='auto')
    ax.contour(d.X, d.Y, dz, [-1, 1], colors=['b', 'r'], linewidths=0.5)
    if topo is not None:
        ax.contour(topo['X'], topo['Y'], topo['B'], [0], colors='k', linewidths=0.5)
    ax.plot(tr[:, 0], tr[:, 1], color='purple', lw=0.8, label='海溝軸(GEBCO)')
    viz.draw_fault(ax, case)
    if case == 'V2025':
        ev = (160.2395, 52.4948)   # USGS 震央
        ax.plot(*ev, '*', color='gold', mec='k', ms=14, label='2025年震央(USGS)')
    ax.plot(159.779, 52.623, 'P', color='lime', mec='k', ms=9, label='1952年震央(USGS ComCat/ISC-GEM)')
    ax.set_xlim(140, 170); ax.set_ylim(38, 60)
    viz.set_lonlat_aspect(ax, 49)
    t = {'V2025': '検証: 2025 M8.8 USGS FFM'}.get(case, f"ケース{case}: {C.SCENARIOS.get(case, {}).get('L', 0)/1e3:.0f}×"
                                                  f"{C.SCENARIOS.get(case, {}).get('W', 0)/1e3:.0f} km, "
                                                  f"すべり{C.mean_slip(**C.SCENARIOS[case]) if case in C.SCENARIOS else 0:.1f} m")
    ax.set_title(t + f'\n上下変位 最大 {dz.max():.1f} m / 最小 {dz.min():.1f} m', fontsize=10)
    fig.colorbar(pc, ax=ax, shrink=0.6, orientation='horizontal', pad=0.06).set_label('海底上下変位 [m]')
axs[0].legend(fontsize=7, loc='lower right')
axs[3].legend(fontsize=7, loc='lower right')
fig.suptitle('断層モデルと初期海底上下変位（Okada 1985）  ' + viz.DISCLAIMER, fontsize=11, y=1.03)
fig.tight_layout()
fn = os.path.join(C.RESULTS, 'fault_dtopo.png')
fig.savefig(fn, dpi=130, bbox_inches='tight')
print('wrote', fn)
