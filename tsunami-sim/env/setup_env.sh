#!/bin/bash
# 計算環境の導入（Ubuntu 24.04 / Python 3.13 で確認）
set -e
apt-get install -y gfortran ffmpeg || (apt-get update && apt-get install -y gfortran ffmpeg)
pip3 install numpy scipy matplotlib netCDF4 xarray six
# clawpack (GeoClaw) v5.13.1 を git から取得し、Python パッケージはシンボリックリンクで参照
mkdir -p /opt/clawpack_src && cd /opt/clawpack_src
[ -d clawpack ] || git clone --branch v5.13.1 https://github.com/clawpack/clawpack.git
cd clawpack && git submodule init && git submodule update
for p in clawutil visclaw geoclaw amrclaw riemann pyclaw; do ln -sfn ../$p/src/python/$p clawpack/$p; done
