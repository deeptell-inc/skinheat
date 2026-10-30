"""Disorder: (1) on-site disorder leaves T_RL/T_LR and dQf/dQr exactly invariant (tridiagonal identity does not
involve the diagonal); (2) bond disorder (kappa_j, phi_j) gives ln(dQf/dQr) = 2 sum_j ln|t_-,j / t_+,j| exactly."""
import numpy as np
from skinheat import chain_heff, delta_heat
from skinheat.transmission import chain_heff_bonds


def test_onsite_disorder_invariance():
    rng = np.random.default_rng(1); L, kappa = 30, 0.1
    tp, tm = -1 - 0.5j * kappa * 1j, -1 + 0.5j * kappa * 1j
    for Wd in [0.5, 2.0, 5.0]:
        f, r, e, TRL, TLR = delta_heat(chain_heff(L, kappa=kappa, onsite=rng.uniform(-Wd / 2, Wd / 2, L)), 1.0, 0.5, ne=801)
        assert abs(np.log(f / r) - 2 * (L - 1) * np.log(abs(tm / tp))) < 1e-9
        assert np.ptp(np.log(TRL / TLR)) < 1e-9


def test_bond_disorder_sum_rule():
    rng = np.random.default_rng(2); L = 30
    kap = 0.1 * (1 + 0.5 * rng.uniform(-1, 1, L - 1)); phi = np.pi / 2 + 0.8 * rng.uniform(-1, 1, L - 1)
    phi[::3] *= -1                                     # mixed chirality
    f, r, e, TRL, TLR = delta_heat(chain_heff_bonds(L, kap, phi, onsite=rng.uniform(-1, 1, L)), 1.0, 0.5, ne=801)
    pred = 2 * sum(np.log(abs((-1 - 0.5j * k * np.exp(-1j * p)) / (-1 - 0.5j * k * np.exp(1j * p)))) for k, p in zip(kap, phi))
    assert abs(np.log(f / r) - pred) < 1e-9
