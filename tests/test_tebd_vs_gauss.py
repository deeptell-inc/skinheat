"""Superfermion MPS Lindbladian (skinheat.lindblad_mps) vs exact references:
 (1) exact null vector of the MPO Liouvillian vs full Fock-space Lindbladian, WITH interaction V (N=4 modes);
 (2) W^II time evolution at V=0 vs the Gaussian Lyapunov NESS (skinheat.ness), small chain + mesoscopic leads."""
import sys, os, warnings, numpy as np
warnings.filterwarnings("ignore")
import tenpy.linalg.np_conserved as npc
from tenpy.algorithms.exact_diag import ExactDiag
from skinheat.model import Model, build_chain_with_leads
from skinheat.ness import lead_currents
from skinheat.lindblad_mps import LindbladMPS, layout
sys.path.insert(0, os.path.dirname(__file__))
from test_brute_force import fock_ops, lindblad_ness


def test_mapping_exact_with_interaction():
    N, L, K, V, kappa, phi = 4, 2, 1, 0.7, 0.3, np.pi / 2
    h = np.zeros((N, N), complex); h[0, 1] = h[1, 0] = -1.0
    h[2, 2] = -0.6; h[3, 3] = 0.9; h[2, 0] = h[0, 2] = 0.4; h[3, 1] = h[1, 3] = 0.5
    l = np.zeros(N, complex); l[0] = 1; l[1] = np.exp(1j * phi)
    loss, gain = [(np.sqrt(kappa) * l, "bulk")], []
    for k, lab, f in [(2, "L", 0.7), (3, "R", 0.3)]:
        v = np.zeros(N, complex); v[k] = 1
        loss.append((np.sqrt(0.8 * (1 - f)) * v, lab)); gain.append((np.sqrt(0.8 * f) * v, lab))
    m = Model(h=h, loss=loss, gain=gain, chain_sites=np.arange(L))
    # brute force
    c = fock_ops(N)
    H = sum(h[i, j] * c[i].conj().T @ c[j] for i in range(N) for j in range(N)) + V * (c[0].conj().T @ c[0]) @ (c[1].conj().T @ c[1])
    Ls = {}
    for l_, lab in loss: Ls.setdefault(lab, []).append(sum(l_[j] * c[j] for j in range(N)))
    for g_, lab in gain: Ls.setdefault(lab, []).append(sum(g_[j] * c[j].conj().T for j in range(N)))
    rho = lindblad_ness(H, [Lo for v in Ls.values() for Lo in v])
    Cbf = np.array([[np.trace(rho @ c[i].conj().T @ c[j]) for j in range(N)] for i in range(N)])
    Ebf = {lab: np.real(np.trace(H @ sum(Lo @ rho @ Lo.conj().T - 0.5 * (Lo.conj().T @ Lo @ rho + rho @ Lo.conj().T @ Lo) for Lo in Ls[lab]))) for lab in "LR"}
    # exact null vector of the vectorised Liouvillian (charge sector N_ket - N_bra = 0)
    S = LindbladMPS(m, pos=layout(L, K), V=V, int_bonds=[(0, 1)])
    ed = ExactDiag(S.M, charge_sector=[0]); ed.build_full_H_from_mpo(); Hf = ed.full_H
    w, v = np.linalg.eig(-1j * Hf.to_ndarray()); k0 = np.argmin(abs(w))
    assert abs(w[k0]) < 1e-10
    S.psi = ed.full_to_mps(npc.Array.from_ndarray(v[:, k0], [Hf.legs[0]], labels=['ps']))
    Cm = S.correlations([(i, j) for i in range(N) for j in range(N)])
    assert np.allclose(Cm, Cbf, atol=1e-10)
    cur = S.currents()
    for lab in "LR":
        assert abs(cur[lab]["E"] - Ebf[lab]) < 1e-10
    # time evolution towards the same NESS (dt^2 error)
    S2 = LindbladMPS(m, pos=layout(L, K), V=V, int_bonds=[(0, 1)]).evolve(30, dt=0.05, chi=64)
    assert abs(S2.correlations([(i, j) for i in range(N) for j in range(N)]) - Cbf).max() < 5e-4


def test_free_fermions_vs_gauss():
    L, K = 6, 4
    m = build_chain_with_leads(L, kappa=0.2, phi=np.pi / 2, K=K, TL=1.0, TR=0.5)
    g = lead_currents(m)
    S = LindbladMPS(m, pos=layout(L, K)).evolve(30, dt=0.05, chi=64)
    c = S.currents()
    for lab in "LR":
        assert abs(c[lab]["E"] - g[lab]["E"]) < 2e-3 * max(1.0, abs(g[lab]["E"]) / 0.05), (lab, c[lab], g[lab])
        assert abs(c[lab]["N"] - g[lab]["N"]) < 2e-3


if __name__ == "__main__":
    test_mapping_exact_with_interaction(); test_free_fermions_vs_gauss(); print("OK")
