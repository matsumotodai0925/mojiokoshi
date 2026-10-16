"""
GEBCO 2026 (15秒) から GeoClaw 用の地形ファイルを作る。

  work/topo/topo_4min.tt3      : 計算領域全体 4分（16x16 ブロック平均）
  work/topo/topo_1min_*.tt3    : 入れ子 Level2 範囲 1分（4x4 ブロック平均）
  work/topo/topo_15s_*.tt3     : 入れ子 Level3 範囲 15秒（GEBCO原データ）
  data/trench_trace.csv        : GEBCO から抽出した千島・カムチャッカ海溝軸
  data/gauges_adjusted.csv     : 観測点を最寄りの海域格子点に移した位置

使い方: python3 02_make_topo.py
"""
import os
import numpy as np
from netCDF4 import Dataset
from clawpack.geoclaw import topotools
import common as C

ds = Dataset(C.GEBCO_NC)
LAT = ds.variables['lat'][:].data
LON = ds.variables['lon'][:].data
ELEV = ds.variables['elevation']
D = 1.0 / 240.0   # 15 arc-second


def read_gebco(x1, x2, y1, y2):
    """[x1,x2]x[y1,y2]（経度 0-360 可）の GEBCO セル中心値を返す。"""
    j1 = int(round((y1 - (-90.0)) / D))
    j2 = int(round((y2 - (-90.0)) / D))
    lats = LAT[j1:j2]

    def lon_idx(x):
        return int(round((x - (-180.0)) / D))

    pieces, lons = [], []
    xa = x1
    while xa < x2 - 1e-9:
        xs = ((xa + 180.0) % 360.0) - 180.0     # -180..180
        seg_end = min(x2, xa + (180.0 - xs))      # 日付変更線で分割
        i1 = lon_idx(xs)
        i2 = i1 + int(round((seg_end - xa) / D))
        pieces.append(ELEV[j1:j2, i1:i2].astype(np.float32))
        lons.append(xa + (np.arange(i2 - i1) + 0.5) * D)
        xa = seg_end
    Z = np.hstack(pieces)
    x = np.concatenate(lons)
    return x, lats, Z


def block_mean(x, y, Z, n):
    ny, nx = Z.shape
    ny2, nx2 = ny // n, nx // n
    Zb = Z[:ny2 * n, :nx2 * n].reshape(ny2, n, nx2, n).mean(axis=(1, 3))
    xb = x[:nx2 * n].reshape(nx2, n).mean(axis=1)
    yb = y[:ny2 * n].reshape(ny2, n).mean(axis=1)
    return xb, yb, Zb


def write_topo(fname, x, y, Z):
    t = topotools.Topography()
    t.set_xyZ(x, y, Z)
    t.write(fname, topo_type=3, Z_format='%.0f')
    print('wrote', fname, Z.shape)


def make_topos():
    d = C.DOMAIN
    # 4 分: 巨大なので緯度帯ごとに読んでブロック平均
    xs, ys, Zs = None, [], []
    step = 5.0
    for y0 in np.arange(d['y1'], d['y2'], step):
        x, y, Z = read_gebco(d['x1'], d['x2'], y0, min(y0 + step, d['y2']))
        xb, yb, Zb = block_mean(x, y, Z, 16)
        xs = xb; ys.append(yb); Zs.append(Zb)
    write_topo(os.path.join(C.TOPO_DIR, 'topo_4min.tt3'), xs, np.concatenate(ys), np.vstack(Zs))

    for name, (x1, x2, y1, y2) in C.L2_BOXES.items():
        # 外側に余白（1度）をとる
        x, y, Z = read_gebco(x1 - 1, x2 + 1, y1 - 1, y2 + 1)
        xb, yb, Zb = block_mean(x, y, Z, 4)
        write_topo(os.path.join(C.TOPO_DIR, f'topo_1min_{name}.tt3'), xb, yb, Zb)

    for name, (x1, x2, y1, y2) in C.L3_BOXES.items():
        x, y, Z = read_gebco(x1 - 0.2, x2 + 0.2, y1 - 0.2, y2 + 0.2)
        write_topo(os.path.join(C.TOPO_DIR, f'topo_15s_{name}.tt3'), x, y, Z)


def trench_trace():
    """各緯度で概略位置±1.5度の範囲の最深点を海溝軸とみなし、平滑化する。"""
    guide = np.array([(56.2, 164.6), (55.5, 164.1), (54.0, 163.2), (52.0, 160.9), (50.0, 158.3),
                      (48.0, 155.6), (46.0, 153.0), (44.0, 150.3), (42.5, 147.6),
                      (41.5, 145.3), (40.5, 144.4)])
    x, y, Z = read_gebco(140.0, 167.0, 40.0, 57.0)
    xb, yb, Zb = block_mean(x, y, Z, 4)     # 1 分
    pts = []
    for j, lat in enumerate(yb):
        if lat < 40.8 or lat > 55.4:
            continue
        lg = np.interp(lat, guide[::-1, 0], guide[::-1, 1])
        hw = 1.5 if lat < 54.0 else 0.6     # 北部は会合部の地形に引かれないよう狭く
        m = (xb > lg - hw) & (xb < lg + hw)
        i = np.argmin(Zb[j, m])
        pts.append((xb[m][i], lat, Zb[j, m][i]))
    pts = np.array(pts)
    # 最深点は谷底の凹凸で数km単位でばらつくため、平滑化スプラインで軸を近似
    # （目標の残差 ~0.1度）。55.4N 以北は会合部で海溝軸が不明瞭なので使わない。
    from scipy.interpolate import UnivariateSpline
    spl = UnivariateSpline(pts[:, 1], pts[:, 0], k=3, s=len(pts) * 0.1 ** 2)
    lon_s = spl(pts[:, 1])
    fn = os.path.join(C.DATA, 'trench_trace.csv')
    with open(fn, 'w') as f:
        f.write('lon_smoothed,lat,lon_deepest,depth_m\n')
        for (lo, la, z), ls in zip(pts, lon_s):
            f.write(f'{ls:.4f},{la:.4f},{lo:.4f},{z:.0f}\n')
    print('wrote', fn, len(pts))


def adjust_gauges():
    """観測点を、15秒 GEBCO で水深 >= MIN_GAUGE_DEPTH の最寄り点へ移動。"""
    fn = os.path.join(C.DATA, 'gauges_adjusted.csv')
    with open(fn, 'w') as f:
        f.write('id,name,lon,lat,depth_m,lon_nominal,lat_nominal,shift_km,kind\n')
        for gid, name, lon, lat, kind in C.GAUGES:
            x, y, Z = read_gebco(lon - 0.3, lon + 0.3, lat - 0.3, lat + 0.3)
            X, Y = np.meshgrid(x, y)
            dist = np.hypot((X - lon) * np.cos(np.radians(lat)), Y - lat) * 111.2
            dist = np.where(-Z >= C.MIN_GAUGE_DEPTH, dist, np.inf)
            j, i = np.unravel_index(np.argmin(dist), dist.shape)
            f.write(f'{gid},{name},{X[j,i]:.5f},{Y[j,i]:.5f},{-Z[j,i]:.0f},'
                    f'{lon},{lat},{dist[j,i]:.2f},{kind}\n')
            print(f'{name:18s} -> ({X[j,i]:.4f},{Y[j,i]:.4f}) depth {-Z[j,i]:.0f} m, shift {dist[j,i]:.2f} km')


if __name__ == '__main__':
    adjust_gauges()
    trench_trace()
    make_topos()
