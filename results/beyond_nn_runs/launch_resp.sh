#!/bin/zsh
# usage: launch_resp.sh <kappa> <L...>  — response-method runs, chi=128, dt=0.1, V in {0,0.5,1,2}
cd "$(dirname "$0")/../.."
PY=${PY:-python}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
kappa=$1; shift
for L in "$@"; do for V in 0.0 0.5 1.0 2.0; do
  nohup $PY scripts/run_response.py run --L $L --V $V --kappa $kappa --chi 128 --dt 0.05 --tbase 60 --tresp 60 > results/beyond_nn_runs/resp_L${L}_V${V}_k${kappa}.log 2>&1 &
done; done
