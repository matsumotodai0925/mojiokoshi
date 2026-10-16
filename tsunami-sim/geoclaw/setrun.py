"""
GeoClaw 実行パラメータ。環境変数で切り替える:
  CASE      : A / B / C / V2025
  MAXLEVEL  : 1 = 粗格子のみ(4分)、3 = 入れ子(4分→1分→15秒)
  T_L3_JSON : (任意) 細格子を強制する開始時刻 {box名: 秒} を書いた JSON。
              粗格子計算の到達時間から 04_run.sh が作る。無ければ 0 秒から。
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import common as C  # noqa: E402
from clawpack.geoclaw import fgmax_tools, fgout_tools  # noqa: E402

CASE = os.environ.get('CASE', 'A')
MAXLEVEL = int(os.environ.get('MAXLEVEL', '1'))
T_L3_JSON = os.environ.get('T_L3_JSON', '')

ARRIVAL_TOL = 0.05   # m: 到達時間の判定しきい値（静水面から +5 cm を初めて超えた時刻）


def setrun(claw_pkg='geoclaw'):
    from clawpack.clawutil import data
    rundata = data.ClawRunData(claw_pkg, 2)
    cd = rundata.clawdata
    d = C.DOMAIN
    cd.lower[0], cd.upper[0] = d['x1'], d['x2']
    cd.lower[1], cd.upper[1] = d['y1'], d['y2']
    dx1 = C.DX1_MIN / 60.0
    cd.num_cells[0] = int(round((d['x2'] - d['x1']) / dx1))
    cd.num_cells[1] = int(round((d['y2'] - d['y1']) / dx1))
    cd.num_eqn, cd.num_aux, cd.capa_index = 3, 3, 2
    cd.t0 = 0.0
    cd.restart = False
    cd.output_style = 2
    cd.output_times = [float(h * 3600) for h in range(2, 13, 2) if h * 3600 < C.TFINAL] + [float(C.TFINAL)]
    cd.tfinal = C.TFINAL
    cd.output_format = 'binary32'
    cd.output_q_components = 'all'
    cd.output_aux_components = 'none'
    cd.verbosity = 0
    cd.dt_variable = True
    cd.dt_initial = 0.2
    cd.dt_max = 1e99
    cd.cfl_desired = 0.75
    cd.cfl_max = 1.0
    cd.steps_max = 1000000
    cd.order = 2
    cd.dimensional_split = 'unsplit'
    cd.transverse_waves = 2
    cd.num_waves = 3
    cd.limiter = ['mc', 'mc', 'mc']
    cd.use_fwaves = True
    cd.source_split = 'godunov'
    cd.num_ghost = 2
    cd.bc_lower = ['extrap', 'extrap']
    cd.bc_upper = ['extrap', 'extrap']
    cd.checkpt_style = 0

    amr = rundata.amrdata
    amr.max1d = 60
    amr.amr_levels_max = MAXLEVEL
    amr.refinement_ratios_x = C.REFINE
    amr.refinement_ratios_y = C.REFINE
    amr.refinement_ratios_t = C.REFINE
    amr.aux_type = ['center', 'capacity', 'yleft']
    amr.flag_richardson = False
    amr.flag2refine = True
    amr.regrid_interval = 3
    amr.regrid_buffer_width = 2
    amr.clustering_cutoff = 0.7
    amr.verbosity_regrid = 0

    # ---- 細格子の範囲 [minlevel, maxlevel, t1, t2, x1, x2, y1, y2]
    regions = rundata.regiondata.regions
    regions.append([1, 1, 0., 1e10, d['x1'] - 1, d['x2'] + 1, d['y1'] - 1, d['y2'] + 1])
    if MAXLEVEL >= 2:
        for b in C.L2_BOXES.values():
            regions.append([1, 2, 0., 1e10] + b)
    if MAXLEVEL >= 3:
        t_l3 = json.load(open(T_L3_JSON)) if T_L3_JSON and os.path.exists(T_L3_JSON) else {}
        for name, b in C.L3_BOXES.items():
            if CASE != 'V2025' and name in ('kushiro', 'ofunato'):
                continue
            t1 = float(t_l3.get(name, 0.0))
            regions.append([1, 3, 0., 1e10] + b)          # 波があれば 15 秒格子まで許可
            regions.append([3, 3, t1, 1e10] + b)          # 到達前から 15 秒格子を強制
        # 沿岸観測点のごく近傍は t=0 から 15 秒格子に固定（格子レベルの切替による水位の跳びを防ぐ）
        for g in C.load_adjusted_gauges():
            if g['kind'] == 'tide' and (CASE == 'V2025' or g['id'] not in (11, 12)):
                regions.append([3, 3, 0., 1e10, g['lon'] - 0.04, g['lon'] + 0.04, g['lat'] - 0.04, g['lat'] + 0.04])

    # ---- 観測点
    for g in C.load_adjusted_gauges():
        if CASE != 'V2025' and g['id'] in (11, 12):
            continue
        rundata.gaugedata.gauges.append([g['id'], g['lon'], g['lat'], 0., 1e10])
    rundata.gaugedata.min_time_increment = 15.0
    rundata.gaugedata.file_format = 'ascii'

    geo = rundata.geo_data
    geo.gravity = 9.81
    geo.coordinate_system = 2
    geo.earth_radius = 6367.5e3
    geo.coriolis_forcing = False
    geo.sea_level = 0.0
    geo.dry_tolerance = 1.e-3
    geo.friction_forcing = True
    geo.manning_coefficient = 0.025
    geo.friction_depth = 1e6

    ref = rundata.refinement_data
    ref.variable_dt_refinement_ratios = True
    ref.wave_tolerance = 0.02

    topo = rundata.topo_data
    topo.topofiles.append([3, os.path.join(C.TOPO_DIR, 'topo_4min.tt3')])
    if MAXLEVEL >= 2:
        for name in C.L2_BOXES:
            topo.topofiles.append([3, os.path.join(C.TOPO_DIR, f'topo_1min_{name}.tt3')])
    if MAXLEVEL >= 3:
        for name in C.L3_BOXES:
            topo.topofiles.append([3, os.path.join(C.TOPO_DIR, f'topo_15s_{name}.tt3')])

    rundata.dtopo_data.dtopofiles.append([3, os.path.join(C.DTOPO_DIR, f'dtopo_{CASE}.tt3')])
    rundata.dtopo_data.dt_max_dtopo = 1.0 if CASE != 'V2025' else 5.0
    rundata.qinit_data.qinit_type = 0

    # ---- fgmax（最大水位・到達時間）
    rundata.fgmax_data.num_fgmax_val = 1
    fgs = rundata.fgmax_data.fgmax_grids

    def fgmax(x1, x2, y1, y2, dx, level):
        fg = fgmax_tools.FGmaxGrid()
        fg.point_style = 2
        fg.x1, fg.x2 = x1 + dx / 2, x2 - dx / 2
        fg.y1, fg.y2 = y1 + dx / 2, y2 - dx / 2
        fg.dx = dx
        fg.tstart_max = 0.
        fg.tend_max = 1e10
        fg.dt_check = 30.
        fg.min_level_check = level
        fg.arrival_tol = ARRIVAL_TOL
        fg.interp_method = 0
        fgs.append(fg)

    fgmax(d['x1'], d['x2'], d['y1'], d['y2'], dx1, 1)                 # fgmax 1: 広域 4 分
    if MAXLEVEL >= 3:
        fgmax(*C.L2_BOXES['hokkaido_tohoku'], 1 / 60., 2)              # fgmax 2: 北海道・東北 1 分
        fgmax(*C.L3_BOXES['hakodate'], 1 / 240., 3)                     # fgmax 3: 函館 15 秒

    # ---- fgout（アニメーション用）
    def fgout(no, x1, x2, y1, y2, dx, dt):
        fo = fgout_tools.FGoutGrid()
        fo.fgno = no
        fo.point_style = 2
        fo.output_format = 'binary32'
        fo.q_out_vars = [1, 4]          # h, eta
        fo.x1, fo.x2, fo.y1, fo.y2 = x1, x2, y1, y2
        fo.nx = int(round((x2 - x1) / dx))
        fo.ny = int(round((y2 - y1) / dx))
        fo.tstart, fo.tend = 0., C.TFINAL
        fo.nout = int(C.TFINAL / dt) + 1
        rundata.fgout_data.fgout_grids.append(fo)

    fgout(1, d['x1'], d['x2'], d['y1'], d['y2'], 0.1, 120.)          # 6 分格子・2 分毎（到達時間の算出にも使う）
    if MAXLEVEL >= 3:
        fgout(2, *C.L3_BOXES['hakodate'], 1 / 240., 120.)               # 函館 15 秒・2 分毎

    amr.dprint = amr.eprint = amr.edebug = amr.gprint = False
    amr.nprint = amr.pprint = amr.rprint = amr.sprint = False
    amr.tprint = False
    amr.uprint = False
    return rundata


if __name__ == '__main__':
    setrun(*sys.argv[1:]).write()
