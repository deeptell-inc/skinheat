#!/bin/zsh
# dt-convergence series: usage launch_dt.sh <dt> <Vs(comma)> <L...>   (kappa=0.05, chi=128, SVD, f/r/iso, tag _dt<dt>)
cd "$(dirname "$0")/../.."
PY=${PY:-python}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
dt=$1; Vs=$2; shift 2
for L in "$@"; do for V in ${(s:,:)Vs}; do for cfg in f r iso; do
  nohup $PY scripts/run_beyond_nn_tebd.py run --L $L --V $V --cfg $cfg --kappa 0.05 --chi 128 --dt $dt --tmax 80 --tmin 40 --comp SVD --tag _dt$dt > results/beyond_nn_runs/L${L}_V${V}_${cfg}_k0.05_dt$dt.log 2>&1 &
done; done; done
