# skinheat

Repository: https://github.com/deeptell-inc/skinheat (private until the Letter is posted).

Code and data for **"Exact exponential heat rectification from the Liouvillian skin effect"** (H. Wakaura and T. Tanimae, submitted to Physical Review Letters).

A reciprocal tight-binding chain with directional bond loss `L_j = √κ (c_j + e^{iφ} c_{j+1})` between two thermal leads rectifies heat with the exact ratio

    δQ_f / δQ_r = |t_-/t_+|^{2(L-1)} = exp[2(L-1)/ξ],      ξ^{-1} = ln|t_-/t_+|,   t_± = -J - (iκ/2) e^{±iφ},

independent of energy, temperature and on-site disorder. This repository contains the exact (Gaussian) solver that proves and checks the law, the matrix-product-state solver used for the interacting chain, every result file behind the figures, and the manuscript.

## Install

```bash
pip install -e ".[dev]"        # numpy, scipy, matplotlib, pytest
pip install -e ".[dev,mps]"    # additionally TeNPy, for the interacting (MPS) runs
```

Python ≥ 3.10. `make test` runs the fast suite (~4 min); `make test-all` adds the slow MPS-vs-exact comparison.

## Package

| Module | Contents |
|---|---|
| `skinheat/model.py` | `Model` (quadratic Lindbladian with labelled jump operators), `build_chain_with_leads` (chain + mesoscopic leads, optional on-site disorder and finite-temperature bulk bath), `fermi` |
| `skinheat/ness.py` | Lyapunov steady state `ness_covariance`, per-reservoir energy/particle currents, `entropy_production` |
| `skinheat/transmission.py` | `chain_heff` / `chain_heff_bonds` (effective non-Hermitian matrix, NNN hopping, bond-resolved disorder), wide-band `transmissions`, ΔT responses `delta_heat`, `rectification` |
| `skinheat/skin.py` | skin profile, localisation length, spectral winding number |
| `skinheat/lindblad_mps.py` | superfermion vectorisation of the Lindbladian on TeNPy (`VectorizedLindbladian`, `LindbladMPS`, W^II propagation, lead currents from two-point functions) |

## Reproduce the figures

Every figure is drawn from `results/*.json` by `paper/make_figures.py` (`make figures`). The JSON files were produced by the scripts below; the exact ones rerun in minutes, the MPS runs take hours each.

| Figure | Content | Script → result file |
|---|---|---|
| Fig. 1 | exact law, scaling collapse, rectification vs throughput | `scripts/run_deltaT_rectification.py` → `results/deltaT_rectification.json` |
| Fig. 2 | next-nearest-neighbour hopping: energy-resolved asymmetry, non-Bloch exponent, winding | `scripts/run_beyond_nn.py` → `results/beyond_nn.json["nnn"]` |
| Fig. 3 | interacting chain (MPS, χ = 128/256) | `scripts/run_beyond_nn_tebd.py run …` (launchers in `results/beyond_nn_runs/launch*.sh`), `… merge`, `scripts/analyze_chi.py` → `results/beyond_nn.json["interacting_*", "chi256_*", "gain_bath_V"]` |
| Fig. 4 | finite-temperature bath and entropy production; on-site and bond disorder | `scripts/run_thermal_bulk.py` → `results/thermal_bulk.json`; `scripts/run_disorder.py` → `results/disorder.json` |

`results/beyond_nn_runs/*.json` hold the time series and final values of every MPS run (file names encode L, V, configuration f/r/iso, κ, χ, dt). `results/NOTES.md` is the lab notebook (Japanese). `scripts/run_response.py` implements a direct linear-response method that was tested and not used for the paper.

## Tests

| Test | What it checks |
|---|---|
| `test_exact_law.py` | the ratio law to machine precision for many (L, κ, φ) |
| `test_brute_force.py` | Lyapunov steady state against the full Fock-space Lindbladian (4 modes) |
| `test_transmission_vs_lyapunov.py` | wide-band Landauer responses against mesoscopic-lead Lyapunov currents |
| `test_tebd_vs_gauss.py` | superfermion mapping exact on 4 modes with interaction; MPS steady state vs exact at V = 0 (slow) |
| `test_thermal_bulk.py` | ratio independent of bath state, σ ≥ 0, σ = 0 at equilibrium, mirror-transpose identity |
| `test_disorder.py` | exact invariance under on-site disorder, bond sum rule |

## Manuscript

`paper/main.tex` (Letter), `paper/sm.tex` (Supplemental Material), `paper/refs.bib`; `make paper` compiles both. `paper/submission-bundle/` holds the submission package (cover letter, declarations, arXiv single-file version, Zenodo deposit files); see its `MANIFEST.md`.

## Citation

See `paper/submission-bundle/zenodo/CITATION.cff`. Code is MIT-licensed (`LICENSE`); data and figures CC BY 4.0.
