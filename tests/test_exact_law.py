"""Exact law: for the nearest-neighbour lossy chain, dQ_f/dQ_r = |t_-/t_+|^{2(L-1)} with
t_+/- = -J - (i kappa/2) e^{+/- i phi}, independent of energy and temperature (tridiagonal Green's-function identity)."""
import numpy as np
from skinheat import chain_heff, delta_heat


def test_exact_nonreciprocity():
    J = 1.0
    for kappa, phi, L in [(0.05, np.pi / 2, 40), (0.3, np.pi / 4, 40), (0.7, 2.0, 30)]:
        tp = -J - 0.5j * kappa * np.exp(1j * phi); tm = -J - 0.5j * kappa * np.exp(-1j * phi)
        pred = 2 * (L - 1) * np.log(abs(tm / tp))
        dQf, dQr, eps, TRL, TLR = delta_heat(chain_heff(L, kappa=kappa, phi=phi), 1.0, 0.5)
        assert abs(np.log(dQf / dQr) - pred) < 1e-8
        assert np.ptp(np.log(TRL / TLR)) < 1e-8
