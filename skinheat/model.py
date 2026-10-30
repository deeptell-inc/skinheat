"""Quadratic Lindbladian model: tight-binding chain + bond loss/gain + mesoscopic thermal leads.

Conventions (fermions, C_ij = <c_i^dag c_j>):
  H       = sum_ij h_ij c_i^dag c_j
  loss    L_mu = sum_j (l_mu)_j c_j
  gain    G_nu = sum_j (g_nu)_j c_j^dag
  dC/dt   = -(X C + C X^dag) + Mg^T,   X = (Ml + Mg^T)/2 - i h^T
  Ml = sum_mu l l^dag,  Mg = sum_nu g g^dag
Each jump vector is tagged with a terminal label so that currents can be resolved per reservoir.
"""
from dataclasses import dataclass, field
import numpy as np


def fermi(e, T, mu=0.0):
    x = (e - mu) / T
    return 0.5 * (1.0 - np.tanh(0.5 * x))  # numerically stable 1/(exp(x)+1)


@dataclass
class Model:
    h: np.ndarray                      # (N,N) Hermitian single-particle matrix (total: chain + leads)
    loss: list = field(default_factory=list)   # list of (vector, label)
    gain: list = field(default_factory=list)   # list of (vector, label)
    chain_sites: np.ndarray = None     # indices of chain sites inside the total index space
    meta: dict = field(default_factory=dict)

    @property
    def N(self):
        return self.h.shape[0]

    def Ml(self):
        M = np.zeros((self.N, self.N), complex)
        for l, _ in self.loss:
            M += np.outer(l, l.conj())
        return M

    def Mg(self):
        M = np.zeros((self.N, self.N), complex)
        for g, _ in self.gain:
            M += np.outer(g, g.conj())
        return M

    def X(self):
        return 0.5 * (self.Ml() + self.Mg().T) - 1j * self.h.T

    def labels(self):
        return sorted({lab for _, lab in self.loss} | {lab for _, lab in self.gain})


def build_chain_with_leads(L, J=1.0, kappa=0.3, phi=np.pi / 2, gain_ratio=0.0,
                           K=40, W=4.0, Gamma=1.0, TL=1.0, TR=0.5, mu=0.0, onsite=0.0, bulk_bath=None):
    """Chain of L sites (hopping -J, reciprocal), bond loss sqrt(kappa)(c_j + e^{i phi} c_{j+1}),
    optional bond gain sqrt(kappa*gain_ratio)(c_j^dag + e^{-i phi} c_{j+1}^dag) [same chirality],
    and two mesoscopic leads of K modes each (band [-W,W]) attached to sites 0 and L-1.
    Lead modes are damped toward Fermi occupation at (T_alpha, mu) with rate equal to the level spacing.

    bulk_bath=dict(K=.., W=.., T=.., mu=..): replace the Lindblad bond loss by a *thermal* mesoscopic bath per bond,
    coupled through the same chiral bond operator A_j = c_j + e^{i phi} c_{j+1} with damping Gamma_b = kappa
    (detailed balance at (T_b, mu_b); mu_b -> -inf recovers the pure-loss model). Bath modes carry label "bulk".
    """
    Kb = 0 if bulk_bath is None else bulk_bath["K"]
    N = L + 2 * K + (L - 1) * Kb
    h = np.zeros((N, N), complex)
    chain = np.arange(L)
    for j in range(L - 1):
        h[j, j + 1] = h[j + 1, j] = -J
    h[chain, chain] = onsite
    loss, gain = [], []
    # bulk directional dissipation
    for j in range(L - 1):
        if bulk_bath is None:
            l = np.zeros(N, complex); l[j] = 1.0; l[j + 1] = np.exp(1j * phi)
            loss.append((np.sqrt(kappa) * l, "bulk"))
            if gain_ratio > 0:
                g = np.zeros(N, complex); g[j] = 1.0; g[j + 1] = np.exp(-1j * phi)
                gain.append((np.sqrt(kappa * gain_ratio) * g, "bulk"))
        else:
            Wb, Tb, mub = bulk_bath["W"], bulk_bath["T"], bulk_bath.get("mu", 0.0)
            eb = np.linspace(-Wb, Wb, Kb); deb = eb[1] - eb[0]
            vb = np.sqrt(kappa * deb / (2 * np.pi))          # Gamma_b = 2 pi v^2 / de = kappa
            base = L + 2 * K + j * Kb
            for k in range(Kb):
                idx = base + k
                h[idx, idx] = eb[k]
                h[idx, j] = vb; h[j, idx] = vb                               # b_k^dag A_j + h.c., A_j = c_j + e^{i phi} c_{j+1}
                h[idx, j + 1] = vb * np.exp(1j * phi); h[j + 1, idx] = vb * np.exp(-1j * phi)
                f = fermi(eb[k], Tb, mub)
                v = np.zeros(N, complex); v[idx] = 1.0
                loss.append((np.sqrt(deb * (1 - f)) * v, "bulk"))
                gain.append((np.sqrt(deb * f) * v, "bulk"))
    # leads
    eps = np.linspace(-W, W, K)
    de = eps[1] - eps[0]
    gamma = de                       # damping = level spacing (smooth spectral function)
    kcoup = np.sqrt(Gamma * de / (2 * np.pi))
    for a, (site, T) in enumerate([(0, TL), (L - 1, TR)]):
        lab = "L" if a == 0 else "R"
        base = L + a * K
        for k in range(K):
            idx = base + k
            h[idx, idx] = eps[k]
            h[idx, site] = h[site, idx] = kcoup
            f = fermi(eps[k], T, mu)
            v = np.zeros(N, complex); v[idx] = 1.0
            loss.append((np.sqrt(gamma * (1 - f)) * v, lab))
            gain.append((np.sqrt(gamma * f) * v, lab))
    return Model(h=h, loss=loss, gain=gain, chain_sites=chain,
                 meta=dict(L=L, J=J, kappa=kappa, phi=phi, gain_ratio=gain_ratio, K=K, W=W,
                           Gamma=Gamma, TL=TL, TR=TR, mu=mu, bulk_bath=bulk_bath,
                           reservoirs={"L": (TL, mu), "R": (TR, mu), **({"bulk": (bulk_bath["T"], bulk_bath.get("mu", 0.0))} if bulk_bath else {})}))
