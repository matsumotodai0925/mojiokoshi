"""
共通設定（パス・計算領域・格子・観測点・シナリオ）。
他のスクリプトはすべてここを import して使う。
"""
import os
import numpy as np

# ---------------------------------------------------------------- paths
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
DATA = os.path.join(ROOT, 'data')
WORK = os.path.join(ROOT, 'work')          # 大きい中間ファイル（git管理外）
TOPO_DIR = os.path.join(WORK, 'topo')
DTOPO_DIR = os.path.join(WORK, 'dtopo')
RUNS = os.path.join(WORK, 'runs')
RESULTS = os.path.join(ROOT, 'results')
GEBCO_NC = os.environ.get('GEBCO_NC', '/opt/gebco/GEBCO_2026.nc')

for d in (WORK, TOPO_DIR, DTOPO_DIR, RUNS, RESULTS):
    os.makedirs(d, exist_ok=True)

# ---------------------------------------------------------------- domain
# 広域計算領域。依頼の表示範囲(135-180E, 30-65N)を含み、
# ヒロ(ハワイ, 204.9E=155.1W)まで届くよう東・南へ拡張している。
# 経度は 0-360 表記（日付変更線をまたぐため）。
DOMAIN = dict(x1=128.0, x2=210.0, y1=15.0, y2=65.0)
MAP_WIDE = dict(x1=135.0, x2=180.0, y1=30.0, y2=65.0)   # 依頼の広域図範囲

# 格子: Level1 = 4分, Level2 = 1分, Level3 = 15秒(GEBCO原解像度)
DX1_MIN = 4.0
REFINE = [4, 4]

# 入れ子（細格子）を許す範囲 [x1, x2, y1, y2]
L2_BOXES = {
    'hokkaido_tohoku': [139.0, 146.8, 38.6, 44.2],
    'kanto':           [138.8, 141.2, 34.3, 36.0],
    'petropavlovsk':   [157.6, 160.4, 51.8, 53.8],
    'hilo':            [203.6, 206.2, 18.6, 20.9],
}
L3_BOXES = {
    'hakodate':      [139.95, 141.45, 41.30, 42.05],   # 津軽海峡東口〜函館湾
    'nemuro':        [145.30, 145.90, 43.05, 43.45],
    'miyako':        [141.85, 142.25, 39.50, 39.80],
    'tokyo_bay_mouth': [139.55, 140.05, 34.85, 35.35],
    'petropavlovsk': [158.45, 159.05, 52.75, 53.15],
    'hilo':          [204.80, 205.25, 19.60, 19.95],
    # 検証用（2025年 M8.8 再現で検潮記録と比較）
    'kushiro':       [144.20, 144.55, 42.85, 43.05],
    'ofunato':       [141.65, 141.95, 38.90, 39.12],
}

# ---------------------------------------------------------------- gauges
# (id, 名称, 公称経度, 公称緯度, 種別)  経度は 0-360
# 沿岸の検潮所位置は粗い格子では陸になることがあるため、
# 02_make_topo.py で「最寄りの水深 >= MIN_GAUGE_DEPTH の点」に移して
# data/gauges_adjusted.csv に保存し、計算ではその位置を使う。
GAUGES = [
    (1, 'Petropavlovsk',   158.650, 53.017, 'tide'),   # IOC: petr
    (2, 'Nemuro(Hanasaki)', 145.570, 43.280, 'tide'),  # IOC: hana (根室市花咲)
    (3, 'Hakodate',        140.720, 41.780, 'tide'),   # IOC: hako
    (4, 'Miyako',          141.975, 39.642, 'tide'),   # 気象庁 宮古検潮所
    (5, 'TokyoBayMouth',   139.740, 35.040, 'point'),  # 東京湾口（浦賀水道南口 洲崎-剱崎間）
    (6, 'Hilo',            204.940, 19.733, 'tide'),   # IOC: hilo
    # 検証・参考用
    (11, 'Kushiro',        144.370, 42.980, 'tide'),   # IOC: kush
    (12, 'Ofunato',        141.750, 39.020, 'tide'),   # IOC: ofun
    (13, 'Mera',           139.830, 34.920, 'tide'),   # IOC: mera (布良, 東京湾口外側)
    (21416, 'DART21416',   163.355, 48.126, 'dart'),
    (21419, 'DART21419',   155.653, 44.401, 'dart'),
    (21418, 'DART21418',   148.753, 38.840, 'dart'),
    (21413, 'DART21413',   152.124, 30.487, 'dart'),
    (21415, 'DART21415',   171.867, 50.120, 'dart'),
    (21414, 'DART21414',   178.199, 48.977, 'dart'),
    (46408, 'DART46408',   190.175, 49.677, 'dart'),
]
MAIN_GAUGES = [1, 2, 3, 4, 5, 6]
MIN_GAUGE_DEPTH = 20.0   # m（M9.9 の引き波で計算点が干上がらないよう 20 m 以深にとる）
IOC_CODE = {1: 'petr', 2: 'hana', 3: 'hako', 6: 'hilo', 11: 'kush', 12: 'ofun', 13: 'mera'}

# ---------------------------------------------------------------- scenarios
MW = 9.9
M0 = 10 ** (1.5 * MW + 9.1)     # N m
MU = 4.0e10                      # 剛性率 Pa（基準値）
MU_RANGE = [3.0e10, 4.0e10, 5.0e10]
DIP = 15.0                       # deg（依頼範囲 15-25 の下限。下端深さを抑えるため）
RAKE = 90.0
TOP_DEPTH = 5.0e3                # 断層上端深さ m
NORTH_END_LAT = 55.0             # 北端：カムチャッカ・アリューシャン会合部の手前（海溝軸が明瞭な北限）
SUBFAULT_LEN = 25.0e3            # 走向方向の分割長 m
N_DOWNDIP = 8                    # 傾斜方向の分割数

SCENARIOS = {
    'A': dict(L=1200e3, W=200e3),
    'B': dict(L=1500e3, W=200e3),
    'C': dict(L=2000e3, W=250e3),
}
SLIP_WARN = 60.0   # m

TFINAL = float(os.environ.get('TFINAL', 12 * 3600.0))   # 試験用に環境変数で短縮可


def mean_slip(L, W, mu=MU, m0=M0):
    return m0 / (mu * L * W)


def load_adjusted_gauges():
    """02_make_topo.py が書いた調整後の観測点位置を読む。"""
    fn = os.path.join(DATA, 'gauges_adjusted.csv')
    out = []
    with open(fn) as f:
        next(f)
        for line in f:
            p = line.strip().split(',')
            out.append(dict(id=int(p[0]), name=p[1], lon=float(p[2]), lat=float(p[3]),
                            depth=float(p[4]), lon0=float(p[5]), lat0=float(p[6]),
                            shift_km=float(p[7]), kind=p[8]))
    return out
