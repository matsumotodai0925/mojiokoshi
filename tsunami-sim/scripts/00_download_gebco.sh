#!/bin/bash
# GEBCO 2026 全球 15 秒格子 (7.47 GB, netCDF) を取得する。
# 長時間転送はプロキシで切断されることがあるため、レジュームしながら繰り返す。
set -u
U=https://dap.ceda.ac.uk/bodc/gebco/global/gebco_2026/ice_surface_elevation/netcdf/GEBCO_2026.nc
F=${GEBCO_NC:-/opt/gebco/GEBCO_2026.nc}
TOTAL=7466018396
mkdir -p "$(dirname "$F")"
for i in $(seq 1 200); do
  S=$(stat -c %s "$F" 2>/dev/null || echo 0); [ "$S" -ge "$TOTAL" ] && break
  curl -sS -C - --max-time 600 -o "$F" "$U" || sleep 3
done
ls -l "$F"
