#!/bin/zsh
# all 24 production runs + one chi-convergence check; logs next to the JSONs
cd "$(dirname "$0")/../.."
PY=${PY:-python}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
for L in 20 40; do for V in 0.0 0.5 1.0 2.0; do for cfg in f r iso; do
  nohup $PY scripts/run_beyond_nn_tebd.py run --L $L --V $V --cfg $cfg --chi 64 > results/beyond_nn_runs/L${L}_V${V}_${cfg}.log 2>&1 &
done; done; done
nohup $PY scripts/run_beyond_nn_tebd.py run --L 20 --V 1.0 --cfg f --chi 32 --tag _chi32 > results/beyond_nn_runs/L20_V1.0_f_chi32.log 2>&1 &
nohup $PY scripts/run_beyond_nn_tebd.py run --L 20 --V 1.0 --cfg f --chi 128 --tag _chi128 > results/beyond_nn_runs/L20_V1.0_f_chi128.log 2>&1 &
