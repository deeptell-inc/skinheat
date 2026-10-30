# skinheat v0.1.0 — code and data for "Exact exponential heat rectification from the Liouvillian skin effect"

H. Wakaura and T. Tanimae (QIRI, Tokyo). Companion deposit for the Letter submitted to Physical Review Letters.

## Contents of `skinheat-v0.1.0.tar.gz`
| Path | What it is |
|---|---|
| `skinheat/` | Python package: quadratic-Lindbladian (Gaussian) solver (`model.py`, `ness.py`), non-Hermitian Landauer transmissions (`transmission.py`), skin-effect diagnostics (`skin.py`), superfermion matrix-product-state Lindblad solver on TeNPy (`lindblad_mps.py`) |
| `scripts/` | Drivers that produced every result file: `run_deltaT_rectification.py` (Fig. 1), `run_beyond_nn.py` (Fig. 2 and the χ=128 series of Fig. 3), `run_beyond_nn_tebd.py` (MPS runs and merge), `analyze_chi.py` (χ/dt convergence), `run_thermal_bulk.py` (Fig. 4a,b), `run_disorder.py` (Fig. 4c,d) |
| `tests/` | pytest suite (exact law, brute-force Fock-space check, Lyapunov vs Landauer, MPS vs exact, thermal bath, disorder) |
| `results/*.json` | Data underlying every figure |
| `results/beyond_nn_runs/*.json` | Time series and final values of every MPS run (file name encodes L, V, configuration, κ, χ, dt) and the launch scripts |
| `paper/make_figures.py`, `paper/figures/` | Regenerates Figs. 1–4 from `results/*.json` alone; the figures as submitted |

## Reproduce
```
pip install -e .            # numpy, scipy, matplotlib; MPS runs additionally need physics-tenpy
python -m pytest tests -q -m "not slow"
python paper/make_figures.py
```
Exact (Gaussian) results regenerate in minutes; each χ=256 MPS run takes 4–16 hours (`results/beyond_nn_runs/launch_chi256.sh`).

## License
Code: MIT (`LICENSE`). Data and figures: CC BY 4.0.
