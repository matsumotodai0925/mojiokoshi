"""可視化・読み込みの共通関数。"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from clawpack.geoclaw import fgmax_tools  # noqa: E402
from clawpack.pyclaw.gauges import GaugeSolution  # noqa: E402
import common as C  # noqa: E402

for f in ('Noto Sans CJK JP', 'IPAexGothic', 'WenQuanYi Zen Hei'):
    if any(f == x.name for x in font_manager.fontManager.ttflist):
        plt.rcParams['font.family'] = f
        break
plt.rcParams['axes.unicode_minus'] = False

DISCLAIMER = 'M9.9仮想シナリオの数値計算（予測・警報の代替ではない）'
DISCLAIMER_V = '2025年7月M8.8の再現計算（モデル検証用。予測・警報の代替ではない）'


def disclaimer(case):
    return DISCLAIMER_V if case == 'V2025' else DISCLAIMER
GAUGE_JA = {1: 'ペトロパブロフスク', 2: '根室(花咲)', 3: '函館', 4: '宮古', 5: '東京湾口',
            6: 'ヒロ(ハワイ)', 11: '釧路', 12: '大船渡', 13: '布良'}


def outdir(case, level):
    return os.path.join(C.RUNS, f'{case}_L{level}', '_output')


def resdir(case, level):
    d = os.path.join(C.RESULTS, f'{case}_L{level}')
    os.makedirs(d, exist_ok=True)
    return d


def read_fgmax(case, level, fgno):
    out = outdir(case, level)
    fg = fgmax_tools.FGmaxGrid()
    fg.read_fgmax_grids_data(fgno, os.path.join(out, 'fgmax_grids.data'))
    fg.read_output(outdir=out, verbose=False)
    B = np.ma.filled(np.ma.array(fg.B, dtype=float), np.nan)
    h = np.ma.filled(np.ma.array(fg.h, dtype=float), np.nan)
    t = np.ma.filled(np.ma.array(fg.arrival_time, dtype=float), np.nan)
    # 海域のみ: 地震前・地震後とも海底が平均海面より下のセル。
    # 地殻変動で沈降して新たに海になった元の陸地（浸水・遡上に相当）は除く。
    B0 = B - dz_at(case, fg.X, fg.Y)
    eta = np.where((B < 0) & (B0 < 0), h + B, np.nan)
    return fg.X, fg.Y, B, eta, t


_DZ = {}


def dz_at(case, X, Y):
    """初期変位(dtopo)の最終時刻の上下変位を (X, Y) に双一次補間。"""
    from scipy.interpolate import RegularGridInterpolator
    from clawpack.geoclaw import dtopotools
    if case not in _DZ:
        d = dtopotools.DTopography(os.path.join(C.DTOPO_DIR, f'dtopo_{case}.tt3'), dtopo_type=3)
        _DZ[case] = RegularGridInterpolator((d.y, d.x), d.dZ[-1], bounds_error=False, fill_value=0.0)
    return _DZ[case]((Y, X))


def read_gauge(case, level, gid):
    g = GaugeSolution(gid, path=outdir(case, level))
    return g.t, g.q[3], g.q[0]   # t[s], eta[m], h[m]


def gauge_location(case, level, gid):
    """その計算で実際に使った計算点の位置（gauge ファイルのヘッダ）。"""
    import re
    with open(os.path.join(outdir(case, level), f'gauge{gid:05d}.txt')) as f:
        m = re.search(r'location=\(\s*(\S+)\s+(\S+)\s*\)', f.readline())
    return float(m.group(1)), float(m.group(2))


def gauge_stats(t, eta, h, tol=0.05, t_ref=120.0):
    """到達時間（初めて |eta-eta0| が tol を超えた時刻; eta0 は地殻変動直後の値）、最大・最小水位。"""
    i0 = np.searchsorted(t, 300.0)            # 地殻変動(<=300s)後を基準水位に
    eta0 = eta[min(i0, len(eta) - 1)] if t[-1] > 300 else 0.0
    wet = h > 0.01
    if wet.mean() < 0.5:
        # この格子解像度では計算点が陸（常時または大半が乾燥）: 値を出さない
        return dict(eta0=np.nan, t_arrival=np.nan, t_arrival_pos=np.nan, eta_max=np.nan,
                    t_max=np.nan, eta_min=np.nan, dry_frac=1 - wet.mean(), dried=True)
    d = np.where(wet, eta - eta0, 0.0)
    # 到達: 直近 10 分間で tol 以上変化した最初の時刻（arrival_from_fgout と同じ定義）
    t_arr = np.nan
    i0 = int(np.searchsorted(t, t_ref))    # 地殻変動完了後から判定
    j = i0
    for i in range(i0, len(t)):
        while t[j] < t[i] - ARR_WIN:
            j += 1
        if wet[i] and np.ptp(eta[j:i + 1][wet[j:i + 1]]) > tol:
            t_arr = t[i]
            break
    k2 = np.where(d > tol)[0]
    t_pos = t[k2[0]] if len(k2) else np.nan
    imax = np.argmax(np.where(wet, eta, -np.inf))
    return dict(eta0=eta0, t_arrival=t_arr, t_arrival_pos=t_pos, eta_max=eta[imax],
                t_max=t[imax], eta_min=np.min(np.where(wet, eta, np.inf)), dry_frac=1 - wet.mean(),
                dried=bool((h[i0:] < 0.1).any()))   # 水深 10 cm 未満 = 実質的に干出


def fault_outline(case):
    """断層の地表投影の外形（上端線と下端線）。"""
    fn = os.path.join(C.DATA, f'fault_{case}.csv')
    if not os.path.exists(fn):
        return None
    a = np.genfromtxt(fn, delimiter=',', names=True)
    top = a[a['j_dip'] == 0]
    last = a[a['j_dip'] == a['j_dip'].max()]
    prev = a[a['j_dip'] == a['j_dip'].max() - 1]
    bot = last.copy()          # 最下列の上端 + 1 列分 = 断層下端
    for k in ('lon_topcenter', 'lat_topcenter'):
        bot[k] = 2 * last[k] - prev[k]
    return top, bot


def draw_fault(ax, case, **kw):
    o = fault_outline(case)
    if o is None:
        return
    top, bot = o
    xs = np.concatenate([top['lon_topcenter'], bot['lon_topcenter'][::-1], top['lon_topcenter'][:1]])
    ys = np.concatenate([top['lat_topcenter'], bot['lat_topcenter'][::-1], top['lat_topcenter'][:1]])
    ax.plot(xs, ys, **({'color': 'k', 'lw': 1.2, 'ls': '--'} | kw))


def set_lonlat_aspect(ax, lat0):
    ax.set_aspect(1.0 / np.cos(np.radians(lat0)))


def lon_label(x, nd=0):
    return f'{x:.{nd}f}°E' if x <= 180 else f'{360 - x:.{nd}f}°W'


ARR_WIN = 600.0     # s


def arr_tol(case):
    """到達判定のしきい値 [m]。M9.9 シナリオは 0.2 m（気象庁 津波注意報の下限に相当する「有意な津波」）。
    5〜10 cm では、Okada 静的変位の緩和による長周期の先行変動（〜10-15 cm）を到達と誤判定するため。
    検証(M8.8)は静的変位が小さいので 0.05 m。"""
    return 0.05 if case == 'V2025' else 0.2


def arrival_from_fgout(case, level, fgno, tol=None, win=ARR_WIN):
    """fgout スナップショット（2 分毎）から到達時間を求める。
    定義: 水位が直近 10 分間で tol 以上変化した最初の時刻（tol は arr_tol(case)）。
    M9.9 では Okada の静的変位が 1000 km 先でも 10〜20 cm 残り、それが非常に長周期の
    ゆっくりした水位変化（〜0.1 m/時）として先行するため、平均海面や地殻変動直後の水位からの
    差では「到達」を誤判定する。そこで変化の速さで波の先端を捉える。
    判定は地殻変動完了後のスナップショットから開始する。返り値: X, Y, t_arr[s](未到達・陸は nan), B"""
    import glob
    from collections import deque
    tol = arr_tol(case) if tol is None else tol
    from clawpack.geoclaw import fgout_tools
    out = outdir(case, level)
    n = len(glob.glob(os.path.join(out, f'fgout{fgno:04d}.t*')))
    fg = fgout_tools.FGoutGrid(fgno, out, 'binary32')
    fg.read_fgout_grids_data(data_file=os.path.join(out, 'fgout_grids.data'))
    buf = deque()
    t_arr = None
    t_ref = 360.0 if case == 'V2025' else 120.0   # 地殻変動の完了後から判定（瞬時破壊 / 検証は破壊継続 ~4.5 分）
    for k in range(1, n + 1):
        f = fg.read_frame(k)
        if f.t < t_ref - 1:
            continue
        if t_arr is None:
            B = f.eta - f.h
            wet0 = f.h > 0.01
            t_arr = np.full(f.eta.shape, np.nan)
        buf.append((f.t, f.eta.copy()))
        while buf[0][0] < f.t - win - 1:
            buf.popleft()
        ch = np.zeros_like(f.eta)
        for tb, eb in buf:
            ch = np.maximum(ch, np.abs(f.eta - eb))
        hit = np.isnan(t_arr) & wet0 & (f.h > 0.01) & (ch > tol)
        t_arr[hit] = f.t
    return f.X, f.Y, t_arr, B


def save_npz(fn, **arrs):
    """図の再描画用の数値を float32 で圧縮保存（リポジトリを軽くするため）。"""
    np.savez_compressed(fn, **{k: np.asarray(v, dtype=np.float32) for k, v in arrs.items()})
