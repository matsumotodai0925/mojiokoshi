# source env/env.sh で計算環境の変数を設定する
export CLAW=${CLAW:-/opt/clawpack_src/clawpack}
export PYTHONPATH=$CLAW${PYTHONPATH:+:$PYTHONPATH}
export FC=gfortran
export FFLAGS="-O2 -fopenmp"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$(nproc)}
export GEBCO_NC=${GEBCO_NC:-/opt/gebco/GEBCO_2026.nc}
