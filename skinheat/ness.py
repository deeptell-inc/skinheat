"""NESS covariance via Lyapunov equation and reservoir-resolved energy/particle currents."""
import numpy as np
from scipy.linalg import solve_continuous_lyapunov


def ness_covariance(model):
    X = model.X()
    Q = model.Mg().T
    C = solve_continuous_lyapunov(X, Q)        # X C + C X^dag = Mg^T
    return 0.5 * (C + C.conj().T)


def _loss_energy(l, h, C):
    # Re <L^dag [H, L]> = -Re[ l^dag C h^T l ]
    return -np.real(l.conj() @ C @ (h.T @ l))


def _gain_energy(g, h, C):
    # Re[ g^dag h g - (h g)^T C g^* ]
    hg = h @ g
    return np.real(g.conj() @ hg - hg @ C @ g.conj())


def energy_current(model, C, label):
    """Energy current flowing from reservoir `label` INTO the total system (= Tr[H D_label(rho)])."""
    h = model.h
    J = 0.0
    for l, lab in model.loss:
        if lab == label:
            J += _loss_energy(l, h, C)
    for g, lab in model.gain:
        if lab == label:
            J += _gain_energy(g, h, C)
    return J


def particle_current(model, C, label):
    I = np.eye(model.N)
    J = 0.0
    for l, lab in model.loss:
        if lab == label:
            J += _loss_energy(l, I, C)
    for g, lab in model.gain:
        if lab == label:
            J += _gain_energy(g, I, C)
    return J


def lead_currents(model, C=None):
    if C is None:
        C = ness_covariance(model)
    out = {}
    for lab in model.labels():
        out[lab] = dict(E=energy_current(model, C, lab), N=particle_current(model, C, lab))
    return out


def entropy_production(model, currents):
    """sigma = -sum_alpha (E_alpha - mu_alpha N_alpha)/T_alpha over reservoirs with a temperature
    (model.meta["reservoirs"][label] = (T, mu)); a T=0 particle sink (pure bond loss) has none and is skipped."""
    s = 0.0
    for lab, (T, mu) in model.meta.get("reservoirs", {}).items():
        s -= (currents[lab]["E"] - mu * currents[lab]["N"]) / T
    return s
