#!/bin/zsh
# chi-convergence series: usage launch_chi.sh <chi> <L...>   (kappa=0.05, V in {0,2}, f/r/iso, SVD, dt=0.1, tag _chi<chi>)
cd "$(dirname "$0")/../.."
PY=${PY:-python}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
chi=$1; shift
for L in "$@"; do for V in 0.0 2.0; do for cfg in f r iso; do
  nohup $PY scripts/run_beyond_nn_tebd.py run --L $L --V $V --cfg $cfg --kappa 0.05 --chi $chi --tmax 80 --tmin 40 --comp SVD --tag _chi$chi > results/beyond_nn_runs/L${L}_V${V}_${cfg}_k0.05_chi$chi.log 2>&1 &
done; done; done
