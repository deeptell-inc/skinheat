"""Non-Hermitian (lossy) Landauer-Buettiker transmission for the chain with wide-band leads.

Single-particle propagator of c-amplitudes is exp(-Z t), Z = X^* = (Ml^T + Mg)/2 + i h,
so the effective non-Hermitian Hamiltonian is  h_eff = h - (i/2)(Ml^T + Mg).
With wide-band leads (self-energy -i Gamma/2 on the end sites):
  T_RL(e) = Gamma^2 |G_{L-1,0}(e)|^2   (left -> right),   T_LR(e) = Gamma^2 |G_{0,L-1}(e)|^2.
"""
import numpy as np
from .model import fermi


def chain_heff(L, J=1.0, kappa=0.3, phi=np.pi / 2, gain_ratio=0.0, Gamma=1.0, onsite=0.0, t2=0.0):
    """t2: reciprocal next-nearest-neighbour hopping (+t2 c^dag_j c_{j+2} + h.c.); breaks tridiagonality."""
    h = np.zeros((L, L), complex)
    for j in range(L - 1):
        h[j, j + 1] = h[j + 1, j] = -J
    for j in range(L - 2):
        h[j, j + 2] = h[j + 2, j] = t2
    h[np.arange(L), np.arange(L)] = onsite
    Ml = np.zeros((L, L), complex); Mg = np.zeros((L, L), complex)
    for j in range(L - 1):
        l = np.zeros(L, complex); l[j] = 1; l[j + 1] = np.exp(1j * phi)
        Ml += kappa * np.outer(l, l.conj())
        if gain_ratio > 0:
            g = np.zeros(L, complex); g[j] = 1; g[j + 1] = np.exp(-1j * phi)
            Mg += kappa * gain_ratio * np.outer(g, g.conj())
    heff = h - 0.5j * (Ml.T + Mg)
    heff[0, 0] -= 0.5j * Gamma
    heff[L - 1, L - 1] -= 0.5j * Gamma
    return heff


def chain_heff_bonds(L, kappa, phi, J=1.0, Gamma=1.0, onsite=0.0, t2=0.0):
    """Bond-resolved version of chain_heff: kappa[j], phi[j] for bond j (bond disorder / mixed chirality)."""
    kappa, phi = np.broadcast_to(kappa, (L - 1,)), np.broadcast_to(phi, (L - 1,))
    h = np.zeros((L, L), complex)
    for j in range(L - 1):
        h[j, j + 1] = h[j + 1, j] = -J
    for j in range(L - 2):
        h[j, j + 2] = h[j + 2, j] = t2
    h[np.arange(L), np.arange(L)] = onsite
    Ml = np.zeros((L, L), complex)
    for j in range(L - 1):
        l = np.zeros(L, complex); l[j] = 1; l[j + 1] = np.exp(1j * phi[j])
        Ml += kappa[j] * np.outer(l, l.conj())
    heff = h - 0.5j * Ml.T
    heff[0, 0] -= 0.5j * Gamma; heff[L - 1, L - 1] -= 0.5j * Gamma
    return heff


def transmissions(heff, eps, Gamma=1.0):
    L = heff.shape[0]
    I = np.eye(L)
    TRL = np.empty(len(eps)); TLR = np.empty(len(eps))
    for n, e in enumerate(eps):
        G = np.linalg.inv(e * I - heff)
        TRL[n] = Gamma ** 2 * abs(G[L - 1, 0]) ** 2
        TLR[n] = Gamma ** 2 * abs(G[0, L - 1]) ** 2
    return TRL, TLR


def delta_heat(heff, Th, Tc, Gamma=1.0, mu=0.0, W=4.0, ne=4001):
    """Linear-in-f response of the heat received by the cold lead when the other lead is heated:
    forward (hot=left): dQ_f = int de/2pi e T_RL (f_h - f_c);  reverse: dQ_r with T_LR."""
    eps = np.linspace(-W, W, ne)
    TRL, TLR = transmissions(heff, eps, Gamma)
    df = fermi(eps, Th, mu) - fermi(eps, Tc, mu)
    w = (eps - mu) * df / (2 * np.pi)
    return np.trapezoid(w * TRL, eps), np.trapezoid(w * TLR, eps), eps, TRL, TLR


def rectification(dQf, dQr):
    return (abs(dQf) - abs(dQr)) / (abs(dQf) + abs(dQr) + 1e-300)
