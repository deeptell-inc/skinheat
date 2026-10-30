"""One runnable check: Gaussian Lyapunov NESS + current formulas vs full Fock-space Lindbladian (N=4)."""
import numpy as np
import itertools
from skinheat.model import Model
from skinheat.ness import ness_covariance, energy_current, particle_current


def fock_ops(N):
    dim = 2 ** N
    c = []
    for j in range(N):
        op = np.zeros((dim, dim), complex)
        for s in range(dim):
            if (s >> j) & 1:
                sign = (-1) ** bin(s & ((1 << j) - 1)).count("1")
                op[s ^ (1 << j), s] = sign
        c.append(op)
    return c


def lindblad_ness(H, Ls):
    dim = H.shape[0]
    I = np.eye(dim)
    Lsup = -1j * (np.kron(H, I) - np.kron(I, H.T))
    for Lo in Ls:
        LdL = Lo.conj().T @ Lo
        Lsup += np.kron(Lo, Lo.conj()) - 0.5 * np.kron(LdL, I) - 0.5 * np.kron(I, LdL.T)
    w, v = np.linalg.eig(Lsup)
    k = np.argmin(abs(w))
    rho = v[:, k].reshape(dim, dim)
    rho /= np.trace(rho)
    return rho


def test_ness_and_currents():
    rng = np.random.default_rng(0)
    N = 4
    A = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
    h = 0.5 * (A + A.conj().T)
    loss = [(rng.normal(size=N) + 1j * rng.normal(size=N), "a"),
            ((rng.normal(size=N) + 1j * rng.normal(size=N)) * 0.7, "b")]
    gain = [((rng.normal(size=N) + 1j * rng.normal(size=N)) * 0.6, "a"),
            ((rng.normal(size=N) + 1j * rng.normal(size=N)) * 0.4, "b")]
    m = Model(h=h, loss=loss, gain=gain, chain_sites=np.arange(N))
    C = ness_covariance(m)

    c = fock_ops(N)
    H = sum(h[i, j] * c[i].conj().T @ c[j] for i in range(N) for j in range(N))
    Ls = {}
    for l, lab in loss:
        Ls.setdefault(lab, []).append(sum(l[j] * c[j] for j in range(N)))
    for g, lab in gain:
        Ls.setdefault(lab, []).append(sum(g[j] * c[j].conj().T for j in range(N)))
    rho = lindblad_ness(H, [L for v in Ls.values() for L in v])
    Cbf = np.array([[np.trace(rho @ c[i].conj().T @ c[j]) for j in range(N)] for i in range(N)])
    assert np.allclose(C, Cbf, atol=1e-8), np.abs(C - Cbf).max()

    Nop = sum(ci.conj().T @ ci for ci in c)
    for lab in ["a", "b"]:
        D = sum(L @ rho @ L.conj().T - 0.5 * (L.conj().T @ L @ rho + rho @ L.conj().T @ L) for L in Ls[lab])
        JE_bf = np.real(np.trace(H @ D)); JN_bf = np.real(np.trace(Nop @ D))
        assert abs(energy_current(m, C, lab) - JE_bf) < 1e-8
        assert abs(particle_current(m, C, lab) - JN_bf) < 1e-8
    # conservation in NESS
    assert abs(sum(energy_current(m, C, lab) for lab in "ab")) < 1e-8


if __name__ == "__main__":
    test_ness_and_currents(); print("OK")
