#!/bin/zsh
# usage: launch2.sh <kappa> <tmax> <L...>   — resolvable-regime sweep (2(L-1)/xi <~ 4), chi=128, dt=0.1
cd "$(dirname "$0")/../.."
PY=${PY:-python}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
kappa=$1; tmax=$2; shift 2
for L in "$@"; do for V in 0.0 0.5 1.0 2.0; do for cfg in f r iso; do
  nohup $PY scripts/run_beyond_nn_tebd.py run --L $L --V $V --cfg $cfg --kappa $kappa --chi 128 --tmax $tmax --tmin ${TMIN:-40} --comp SVD > results/beyond_nn_runs/L${L}_V${V}_${cfg}_k${kappa}.log 2>&1 &
done; done; done
