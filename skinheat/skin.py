"""Skin-effect diagnostics of the effective single-particle generator X (restricted to the chain)."""
import numpy as np


def chain_X(model):
    s = model.chain_sites
    return model.X()[np.ix_(s, s)]


def skin_profile(model):
    """Steady-state density profile on the chain and the mean right-eigenvector weight profile of X_chain."""
    Xc = chain_X(model)
    w, V = np.linalg.eig(Xc)
    P = np.abs(V) ** 2
    P /= P.sum(axis=0, keepdims=True)
    return P.mean(axis=1)


def localization_length(model):
    """Fit log of mean eigenvector weight vs site to exp(-x/xi) (xi>0: left edge, xi<0: right edge)."""
    p = skin_profile(model)
    x = np.arange(len(p))
    slope = np.polyfit(x, np.log(p + 1e-300), 1)[0]
    return -1.0 / slope if abs(slope) > 1e-12 else np.inf


def winding_number(model, E=None, nk=2001):
    """Spectral winding of the PBC Bloch generator of the chain part of -iX? We use the
    non-Hermitian effective Hamiltonian H_eff(k) = i X_bulk(k) (so that eigenvalues ~ energies).
    Winding of det[H_eff(k) - E] around E over the BZ; E defaults to the spectral centroid."""
    Xc = chain_X(model)
    L = Xc.shape[0]
    # extract bulk hoppings from the middle of the chain (translation-invariant part)
    j = L // 2
    t_p = Xc[j, j + 1]   # X_{j,j+1}
    t_m = Xc[j + 1, j]   # X_{j+1,j}
    d = Xc[j, j]
    ks = np.linspace(-np.pi, np.pi, nk)
    xk = d + t_p * np.exp(1j * ks) + t_m * np.exp(-1j * ks)
    Hk = 1j * xk
    if E is None:
        E = np.mean(Hk)
    ang = np.unwrap(np.angle(Hk - E))
    return (ang[-1] - ang[0]) / (2 * np.pi)
