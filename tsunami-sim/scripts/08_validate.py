"""
2025-07-29 M8.8 カムチャッカ地震の再現計算（ケース V2025）と観測値の比較。
  観測: NOAA DART（海底水圧計）、IOC 経由の検潮記録（日本の地点は気象庁データ）
  潮汐除去: 2025-07-26〜08-02 の記録に主要 8 分潮＋一次トレンドを最小二乗で当てはめ、
            地震前後（-1h〜+20h）を除いた区間で係数を決めて差し引く。
出力: results/V2025_L<lv>/validation_dart.png, validation_tide.png, validation_summary.md
使い方: python3 08_validate.py LEVEL
"""
import os
import sys
import gzip
import json
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import common as C
import viz

ORIGIN = dt.datetime(2025, 7, 29, 23, 24, 52, 483000, tzinfo=dt.timezone.utc)   # USGS ComCat us6000qw60
PERIODS_H = {'M2': 12.4206012, 'S2': 12.0, 'N2': 12.65834751, 'K2': 11.96723606,
             'K1': 23.93447213, 'O1': 25.81933871, 'P1': 24.06588766, 'Q1': 26.868350}
DARTS = [21416, 21419, 21418, 21413, 21415, 21414, 46408]
TIDES = [1, 11, 2, 3, 12, 13, 6]


def detide(t_h, z):
    """t_h: 地震発生からの時間[h]。潮汐+トレンドを除いた残差を返す。"""
    cols = [np.ones_like(t_h), t_h]
    for p in PERIODS_H.values():
        w = 2 * np.pi / p
        cols += [np.cos(w * t_h), np.sin(w * t_h)]
    A = np.array(cols).T
    fit = (t_h < -1) | (t_h > 20)
    coef, *_ = np.linalg.lstsq(A[fit], z[fit], rcond=None)
    return z - A @ coef


def read_dart(sid):
    fn = os.path.join(C.DATA, 'dart', f'{sid}t2025.txt.gz')
    rows = []
    with gzip.open(fn, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.split()
            h = float(p[7])
            if h > 9000:
                continue
            t = dt.datetime(int(p[0]), int(p[1]), int(p[2]), int(p[3]), int(p[4]), int(p[5]), tzinfo=dt.timezone.utc)
            rows.append(((t - ORIGIN).total_seconds() / 3600, int(p[6]), h))
    a = np.array(rows)
    a = a[(a[:, 0] > -84) & (a[:, 0] < 84)]
    # 同時刻の重複は細かいサンプリング(type 2,3)を優先
    a = a[np.lexsort((-a[:, 1], a[:, 0]))]
    _, idx = np.unique(np.round(a[:, 0] * 3600), return_index=True)
    a = a[idx]
    return a[:, 0], a[:, 2]


def read_ioc(code):
    d = json.load(open(os.path.join(C.DATA, 'tide_ioc', f'{code}.json')))
    t = np.array([(dt.datetime.fromisoformat(x['stime']).replace(tzinfo=dt.timezone.utc) - ORIGIN).total_seconds() / 3600
                  for x in d])
    z = np.array([x['slevel'] for x in d], dtype=float)
    ok = np.isfinite(z)
    t, z = t[ok], z[ok]
    # 明らかな外れ値（中央値から 5 m 以上）を除外
    z_med = np.median(z)
    ok = np.abs(z - z_med) < 5
    return t[ok], z[ok]


def first_exceed(t, z, tol, t0=0.0):
    k = np.where((t > t0) & (np.abs(z) > tol))[0]
    return t[k[0]] if len(k) else np.nan


def compare(level, ids, reader, fname, title, tol, tmax_h=12):
    gauges = {g['id']: g for g in C.load_adjusted_gauges()}
    fig, axs = plt.subplots(len(ids), 1, figsize=(11, 2.1 * len(ids)), sharex=True)
    rows = []
    for ax, gid in zip(axs, ids):
        t_o, z_o = reader(gid)
        r_o = detide(t_o, z_o)
        w = (t_o >= -0.5) & (t_o <= tmax_h)
        r_o -= np.median(r_o[(t_o > -3) & (t_o < 0)])
        t_m, eta_m, h_m = viz.read_gauge('V2025', level, gid)
        t_mh = t_m / 3600
        eta_m = np.where(h_m > 0.01, eta_m, np.nan)
        ax.plot(t_o[w], r_o[w], color='0.3', lw=0.9, label='観測（潮汐除去）')
        ax.plot(t_mh, eta_m, color='#d62728', lw=0.9, label='計算 (GeoClaw)')
        name = gauges[gid]['name']
        ax.set_title(f"{viz.GAUGE_JA.get(gid, name)}", fontsize=9, loc='left')
        ax.axhline(0, color='0.7', lw=0.5); ax.grid(alpha=0.3); ax.set_ylabel('[m]')
        # 指標: 観測は地震動による水圧計ノイズ・検潮所の短周期振動を除くため 5 分移動平均。
        # 計算の到達（tol 超過）30 分前〜3 時間後の窓で、最大水位とその時刻を比べる。
        ok = np.isfinite(eta_m)
        if ok.sum() > 10:
            d_m = eta_m - eta_m[ok][min(np.searchsorted(t_mh[ok], 0.1), ok.sum() - 1)]
            t_arr = first_exceed(t_mh[ok], d_m[ok], tol, 0.1)
        else:
            t_arr = np.nan        # 計算点が陸（この格子解像度では海域にならない）
        row = dict(id=gid, name=viz.GAUGE_JA.get(gid, name), mod_arr=t_arr, obs_max=np.nan, mod_max=np.nan,
                   obs_tmax=np.nan, mod_tmax=np.nan)
        if np.isfinite(t_arr):
            ww = (t_o >= t_arr - 0.5) & (t_o <= t_arr + 3)
            tt, rr = t_o[ww], np.convolve(r_o[ww], np.ones(5) / 5, mode='same')
            wm = ok & (t_mh >= t_arr - 0.5) & (t_mh <= t_arr + 3)
            if len(rr) > 10:
                i, j = np.argmax(rr), np.argmax(eta_m[wm])
                row.update(obs_max=rr[i], obs_tmax=tt[i], mod_max=eta_m[wm][j], mod_tmax=t_mh[wm][j])
                ax.axvspan(t_arr - 0.5, t_arr + 3, color='#ffe9a8', alpha=0.35, lw=0)
        rows.append(row)
    axs[0].legend(fontsize=8, ncol=2, loc='upper right')
    axs[-1].set_xlabel('地震発生 (2025-07-29 23:24:52 UTC) からの時間 [h]')
    axs[-1].set_xlim(-0.5, tmax_h)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(fname, dpi=130); plt.close(fig)
    print('wrote', fname)
    return rows


def main(level):
    rd = viz.resdir('V2025', level)
    r1 = compare(level, DARTS, read_dart, os.path.join(rd, 'validation_dart.png'),
                 f'検証: 2025年7月 M8.8 — DART 観測 vs 計算（格子 L{level}）', 0.02)
    r2 = compare(level, TIDES, lambda g: read_ioc(C.IOC_CODE[g]), os.path.join(rd, 'validation_tide.png'),
                 f'検証: 2025年7月 M8.8 — 検潮所 観測 vs 計算（格子 L{level}）\n'
                 '計算点は検潮所近傍の水深20m以深の点。港内の検潮記録とは位置・格子解像度が異なる', 0.05)
    fn = os.path.join(rd, 'validation_summary.md')
    with open(fn, 'w') as f:
        f.write('| 地点 | 観測 最大水位 [m] | 計算 最大水位 [m] | 計算/観測 | 観測 最大の時刻 (h) | 計算 最大の時刻 (h) | 計算 第1波到達 (h) |\n|---|---|---|---|---|---|---|\n')
        for r in r1 + r2:
            ratio = r['mod_max'] / r['obs_max'] if r['obs_max'] > 0 else np.nan
            f.write(f"| {r['name']} | {r['obs_max']:.2f} | {r['mod_max']:.2f} | {ratio:.2f} | {r['obs_tmax']:.2f} | "
                    f"{r['mod_tmax']:.2f} | {r['mod_arr']:.2f} |\n")
        f.write('\n比較窓（図の黄色帯）: 計算の第1波到達（DART 2 cm / 検潮所 5 cm 超過）の 30 分前〜3 時間後。'
                '観測は潮汐除去後に 5 分移動平均。計算点が陸（乾燥）になった地点は空欄(nan)。\n')
    print(open(fn).read())


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
