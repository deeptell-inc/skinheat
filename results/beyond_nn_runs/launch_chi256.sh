#!/bin/zsh
# chi=256 precision runs: kappa=0.05, SVD, dt=0.05, tmax=60; usage: launch_chi256.sh <L> <Vs comma>
cd "$(dirname "$0")/../.."
PY=${PY:-python}
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2
L=$1; Vs=$2; kappa=${3:-0.05}
for V in ${(s:,:)Vs}; do for cfg in f r iso; do
  nohup $PY scripts/run_beyond_nn_tebd.py run --L $L --V $V --cfg $cfg --kappa $kappa --chi 256 --dt 0.05 --tmax 60 --tmin 40 --comp SVD --tag _chi256_dt0.05 > results/beyond_nn_runs/L${L}_V${V}_${cfg}_k${kappa}_chi256_dt0.05.log 2>&1 &
done; done
