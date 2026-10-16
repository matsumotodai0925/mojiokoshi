"""
初期海底変位（Okada 1985, clawpack.geoclaw.dtopotools）を作る。

  シナリオ A/B/C : M9.9 仮想断層。海溝軸（data/trench_trace.csv）に沿って
                   北端 NORTH_END_LAT から南西へ長さ L、幅 W、傾斜 DIP、すべり角 90°、
                   一様すべり D = M0/(mu L W)。瞬時破壊。
  検証 V2025     : 2025-07-29 M8.8 の USGS 有限断層モデル(complete_inversion.fsp)、
                   破壊開始時刻・立ち上がり時間を考慮したキネマティック破壊。

出力: work/dtopo/dtopo_<case>.tt3, data/fault_<case>.csv, results/fault_parameters.md
使い方: python3 03_make_dtopo.py [A B C V2025] [--report-only]
"""
import os
import sys
import numpy as np
from clawpack.geoclaw import dtopotools
import common as C

R_EARTH = 6371.0e3


def destination(lon, lat, az_deg, dist_m):
    """球面上で (lon,lat) から方位 az に dist 進んだ点。"""
    p1, l1, az = np.radians(lat), np.radians(lon), np.radians(az_deg)
    d = dist_m / R_EARTH
    p2 = np.arcsin(np.sin(p1) * np.cos(d) + np.cos(p1) * np.sin(d) * np.cos(az))
    l2 = l1 + np.arctan2(np.sin(az) * np.sin(d) * np.cos(p1), np.cos(d) - np.sin(p1) * np.sin(p2))
    return np.degrees(l2), np.degrees(p2)


def azimuth(lon1, lat1, lon2, lat2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(lon2 - lon1)
    y = np.sin(dl) * np.cos(p2)
    x = np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl)
    return np.degrees(np.arctan2(y, x)) % 360.0


def haversine(lon1, lat1, lon2, lat2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    return 2 * R_EARTH * np.arcsin(np.sqrt(a))


def load_trench():
    a = np.loadtxt(os.path.join(C.DATA, 'trench_trace.csv'), delimiter=',', skiprows=1)
    a = a[np.argsort(-a[:, 1])]          # 北→南
    lon, lat = a[:, 0], a[:, 1]
    s = np.concatenate([[0], np.cumsum(haversine(lon[:-1], lat[:-1], lon[1:], lat[1:]))])
    s0 = np.interp(-C.NORTH_END_LAT, -lat, s)   # 北端の弧長
    return lon, lat, s - s0


def scenario_fault(name):
    L, W = C.SCENARIOS[name]['L'], C.SCENARIOS[name]['W']
    slip = C.mean_slip(L, W)
    lon_t, lat_t, s_t = load_trench()
    if s_t[-1] < L:
        raise ValueError(f'trench trace too short for L={L/1e3} km')
    nseg = int(round(L / C.SUBFAULT_LEN))
    dw = W / C.N_DOWNDIP
    dip = C.DIP
    fault = dtopotools.Fault()
    fault.subfaults = []
    rows = []
    for k in range(nseg):
        sa, sb = k * C.SUBFAULT_LEN, (k + 1) * C.SUBFAULT_LEN
        la_lon, la_lat = np.interp(sa, s_t, lon_t), np.interp(sa, s_t, lat_t)
        lb_lon, lb_lat = np.interp(sb, s_t, lon_t), np.interp(sb, s_t, lat_t)
        strike = azimuth(la_lon, la_lat, lb_lon, lb_lat)
        mid_lon, mid_lat = destination(la_lon, la_lat, strike, C.SUBFAULT_LEN / 2)
        for r in range(C.N_DOWNDIP):
            off = r * dw * np.cos(np.radians(dip))
            lon, lat = destination(mid_lon, mid_lat, strike + 90.0, off)
            sf = dtopotools.SubFault()
            sf.coordinate_specification = 'top center'
            sf.longitude, sf.latitude = lon, lat
            sf.depth = C.TOP_DEPTH + r * dw * np.sin(np.radians(dip))
            sf.strike, sf.dip, sf.rake = strike, dip, C.RAKE
            sf.length, sf.width, sf.slip = C.SUBFAULT_LEN, dw, slip
            sf.mu = C.MU
            fault.subfaults.append(sf)
            rows.append((k, r, lon, lat, sf.depth, strike, dip, C.RAKE, C.SUBFAULT_LEN, dw, slip))
    return fault, rows, dict(L=L, W=W, slip=slip, nseg=nseg)


def usgs_fault():
    """USGS FFM (.fsp, 複数セグメント) を読み込む。"""
    fn = os.path.join(C.DATA, 'usgs_ffm', 'complete_inversion.fsp')
    fault = dtopotools.Fault()
    fault.subfaults = []
    rows = []
    strike = dip = dx = dz = None
    for line in open(fn):
        if line.startswith('% SEGMENT'):
            p = line.replace('=', ' ').split()
            strike = float(p[p.index('STRIKE') + 1]); dip = float(p[p.index('DIP') + 1])
        elif line.startswith('% Dx') and strike is not None:
            p = line.replace('=', ' ').split()
            dx = float(p[2]); dz = float(p[5])
        elif not line.startswith('%') and line.strip():
            lat, lon, x, y, z, slip, rake, trup, rise, mo = map(float, line.split())
            sf = dtopotools.SubFault()
            sf.coordinate_specification = 'centroid'
            sf.longitude, sf.latitude, sf.depth = lon, lat, z * 1e3
            sf.strike, sf.dip, sf.rake, sf.slip = strike, dip, rake, slip
            sf.length, sf.width = dx * 1e3, dz * 1e3
            sf.rupture_time, sf.rise_time = trup, max(rise, 2.0)
            fault.subfaults.append(sf)
            rows.append((lon, lat, z * 1e3, strike, dip, rake, dx * 1e3, dz * 1e3, slip, trup, rise))
    return fault, rows


def make_dtopo(fault, fname, times=(0.0, 1.0), pad=12.0, dx=3.0 / 60, taper=2.0):
    """Okada 変位を断層の周囲 pad 度の範囲で計算する。
    M9.9 では 3 度離れても変位が数十 cm 残り、格子の端で打ち切ると偽の段差波源になるため、
    範囲を広く取り、さらに外周 taper 度でコサイン型に 0 へ減衰させる（端の値は data/dtopo_edge_check.txt に記録）。"""
    lons = np.array([s.longitude for s in fault.subfaults])
    lats = np.array([s.latitude for s in fault.subfaults])
    d = C.DOMAIN
    x1, x2 = max(np.floor(lons.min() - pad), d['x1']), min(np.ceil(lons.max() + pad), d['x2'])
    y1, y2 = max(np.floor(lats.min() - pad), d['y1']), min(np.ceil(lats.max() + pad), d['y2'])
    x = np.arange(x1, x2 + dx / 2, dx)
    y = np.arange(y1, y2 + dx / 2, dx)
    dtopo = fault.create_dtopography(x, y, times=list(times), verbose=False)
    zf = dtopo.dZ[-1]
    edge = max(abs(zf[0]).max(), abs(zf[-1]).max(), abs(zf[:, 0]).max(), abs(zf[:, -1]).max())

    def w(v, a, b):
        r = np.clip(np.minimum(v - a, b - v) / taper, 0, 1)
        return 0.5 - 0.5 * np.cos(np.pi * r)
    W = np.outer(w(y, y1, y2), w(x, x1, x2))
    dtopo.dZ = dtopo.dZ * W[None, :, :]
    dtopo.write(fname, dtopo_type=3)
    dz = dtopo.dZ[-1]
    print(f'wrote {fname}: [{x1},{x2}]x[{y1},{y2}] grid {dz.shape}, dz max {dz.max():.2f} m, '
          f'min {dz.min():.2f} m, |dz| at edge before taper {edge:.3f} m')
    with open(os.path.join(C.DATA, 'dtopo_edge_check.txt'), 'a') as f:
        f.write(f'{os.path.basename(fname)} extent [{x1},{x2}]x[{y1},{y2}] dx={dx*60:.0f}min '
                f'max|dz| at edge before taper = {edge:.4f} m\n')
    return dtopo


def main(cases):
    report = []
    for case in cases:
        if case in C.SCENARIOS:
            fault, rows, info = scenario_fault(case)
            with open(os.path.join(C.DATA, f'fault_{case}.csv'), 'w') as f:
                f.write('i_strike,j_dip,lon_topcenter,lat_topcenter,depth_top_m,strike,dip,rake,length_m,width_m,slip_m\n')
                for r in rows:
                    f.write(','.join(f'{v:.4f}' if isinstance(v, float) else str(v) for v in r) + '\n')
            fn = os.path.join(C.DTOPO_DIR, f'dtopo_{case}.tt3')
            if REPORT_ONLY and os.path.exists(fn):
                dt = dtopotools.DTopography(fn, dtopo_type=3)
            else:
                dt = make_dtopo(fault, fn)
            strikes = np.array([r[5] for r in rows])
            first, last = rows[0], rows[(info['nseg'] - 1) * C.N_DOWNDIP]
            report.append(dict(case=case, Mw=fault.Mw(), **info,
                               strike_min=strikes.min(), strike_max=strikes.max(),
                               north=(first[2], first[3]), south=(last[2], last[3]),
                               bottom_depth=C.TOP_DEPTH + info['W'] * np.sin(np.radians(C.DIP)),
                               dzmax=dt.dZ[-1].max(), dzmin=dt.dZ[-1].min()))
        elif case == 'V2025':
            fault, rows = usgs_fault()
            fault.rupture_type = 'kinematic'
            tmax = max(r[9] + r[10] for r in rows)
            times = np.arange(0.0, tmax + 20.0, 10.0)
            dt = make_dtopo(fault, os.path.join(C.DTOPO_DIR, 'dtopo_V2025.tt3'), times=times)
            mo = sum(s.slip * s.length * s.width * 4.0e10 for s in fault.subfaults)
            print(f'V2025: {len(rows)} subfaults, max slip {max(r[8] for r in rows):.1f} m, '
                  f'rupture duration {tmax:.0f} s, Mo(mu=4e10)={mo:.3e}')
    if report:
        write_report(report)


def write_report(report):
    fn = os.path.join(C.RESULTS, 'fault_parameters.md')
    with open(fn, 'w') as f:
        f.write('# 断層パラメータ（M9.9 仮想シナリオ）\n\n')
        f.write(f'- M0 = 10^(1.5×{C.MW} + 9.1) = {C.M0:.3e} N·m\n')
        f.write(f'- 計算に使う剛性率 μ = {C.MU:.1e} Pa、傾斜 {C.DIP}°、すべり角 {C.RAKE}°、上端深さ {C.TOP_DEPTH/1e3:.0f} km\n')
        f.write(f'- 断層上端は GEBCO 2026 から抽出した海溝軸に沿わせ、北端は北緯 {C.NORTH_END_LAT}° '
                f'（カムチャッカ・アリューシャン会合部付近）から南西へ配置\n')
        f.write(f'- 分割: 走向方向 {C.SUBFAULT_LEN/1e3:.0f} km × 傾斜方向 {C.N_DOWNDIP} 列、一様すべり、瞬時破壊\n\n')
        f.write('| ケース | 長さ×幅 (km) | 平均すべり μ=3e10 | μ=4e10（採用） | μ=5e10 | 走向の範囲 | 北端 | 南端 | 下端深さ | 海底変位 最大/最小 | 判定 |\n')
        f.write('|---|---|---|---|---|---|---|---|---|---|---|\n')
        for r in report:
            s = [C.mean_slip(r['L'], r['W'], mu) for mu in C.MU_RANGE]
            flag = f'**{C.SLIP_WARN:.0f} m 超：物理的に極端**' if r['slip'] > C.SLIP_WARN else f'{C.SLIP_WARN:.0f} m 以下'
            f.write(f"| {r['case']} | {r['L']/1e3:.0f}×{r['W']/1e3:.0f} | {s[0]:.1f} m | {s[1]:.1f} m | {s[2]:.1f} m | "
                    f"{r['strike_min']:.0f}–{r['strike_max']:.0f}° | {r['north'][1]:.2f}N {r['north'][0]:.2f}E | "
                    f"{r['south'][1]:.2f}N {r['south'][0]:.2f}E | {r['bottom_depth']/1e3:.0f} km | "
                    f"{r['dzmax']:.1f} / {r['dzmin']:.1f} m | {flag} |\n")
        f.write('\n（Mw は dtopotools で各小断層の M0 を合計して再計算した値: ' +
                ', '.join(f"{r['case']}={r['Mw']:.2f}" for r in report) + '）\n')
    print('wrote', fn)


REPORT_ONLY = '--report-only' in sys.argv   # 既存の dtopo を読み、表だけ作り直す

if __name__ == '__main__':
    main([a for a in sys.argv[1:] if not a.startswith('--')] or ['A', 'B', 'C', 'V2025'])
