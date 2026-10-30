"""Interacting Lindblad NESS with matrix product states (tenpy).

Superfermion vectorisation: |rho>> = rho |I>>, |I>> = prod_j (1 + c_j^dag ct_j^dag)|0>, with a second
("tilde"/bra) fermion ct_j per physical mode. For parity-even rho this gives the exact rules
  c_j rho        -> c_j            rho c_j^dag     -> -ct_j
  rho c_j        -> ct_j^dag       rho c_i^dag c_j -> ct_j^dag ct_i
so the Lindbladian becomes a *local* fermionic operator on 2N modes; each physical mode is one
GroupedSite (ket, bra) of dimension 4 and N_ket - N_bra = 0 is conserved.
We evolve exp(dt L) with the W^II MPO propagator (ExpMPOEvolution with H_model := i L).
Expectation values are <<I|O|rho>> / <<I|rho>>.  Lead energy currents only need the lead-mode
two-point functions (the lead jump operators commute with the interaction), so the Gaussian formula
`skinheat.ness.energy_current` is reused with a partially filled correlation matrix.
"""
import numpy as np
import tenpy.linalg.np_conserved as npc
from tenpy.networks.site import FermionSite, GroupedSite
from tenpy.networks.mps import MPS, MPSEnvironment
from tenpy.models.model import CouplingMPOModel
from tenpy.models.lattice import Chain
from tenpy.algorithms.mpo_evolution import ExpMPOEvolution
from .ness import energy_current, particle_current


def layout(L, K):
    """MPS position of each mode of build_chain_with_leads: [chain 0][lead L][chain 1..L-1][lead R].
    Every long-range term then *starts* on a chain site, so the MPO bond dimension stays ~10
    (tenpy keys MPO channels by the starting site of a term)."""
    N = L + 2 * K
    order = [0] + list(range(L, L + K)) + list(range(1, L)) + list(range(L + K, N))
    pos = np.empty(N, int); pos[order] = np.arange(N)
    return pos


def superfermion_site(conserve=True):
    ket = FermionSite(conserve='N' if conserve else None)
    bra = FermionSite(conserve='N' if conserve else None)
    if conserve:   # charge of the bra fermion counts negative -> conserved N_ket - N_bra
        bra.change_charge(npc.LegCharge.from_qflat(bra.leg.chinfo, [[0], [-1]]))
    return GroupedSite([ket, bra], labels=['k', 'b'], charges='same')


class VectorizedLindbladian(CouplingMPOModel):
    """H_model = i * L_vec, so that the tenpy real-time propagator exp(-i dt H_model) = exp(dt L_vec).
    params: model (skinheat.model.Model), pos (mode -> MPS site), V, int_bonds (list of (i,j) mode pairs)."""

    def init_sites(self, p):
        return superfermion_site(p.get('conserve', True))

    def init_lattice(self, p):
        N = p['model'].N if p.get('model', None) is not None else p['N']
        return Chain(N, self.init_sites(p), bc='open', bc_MPS='finite')

    def init_terms(self, p):
        m, pos = p.get('model', None), p['pos']
        V, bonds = p.get('V', 0.0), p.get('int_bonds', [])
        add = lambda c, term: self.add_local_term(c, [(op, [int(pos[i]), 0]) for op, i in term])
        for c, term in p.get('terms', []):      # explicit extra terms (plain operator, no factor i)
            add(c, term)
        if m is None:
            return
        h = m.h
        for i, j in zip(*np.nonzero(h)):
            add(h[i, j], [('Cdk', i), ('Ck', j)])            # i(-i H rho)
            add(-h[i, j], [('Cdb', j), ('Cb', i)])           # i(+i rho H)
        for i, j in bonds:
            add(V, [('Nk', i), ('Nk', j)])
            add(-V, [('Nb', i), ('Nb', j)])
        for l, _ in m.loss:
            for a in np.nonzero(l)[0]:
                for b in np.nonzero(l)[0]:
                    add(-1j * l[a] * np.conj(l[b]), [('Ck', a), ('Cb', b)])           # L rho L^dag
                    add(-0.5j * np.conj(l[a]) * l[b], [('Cdk', a), ('Ck', b)])        # -1/2 L^dag L rho
                    add(-0.5j * np.conj(l[a]) * l[b], [('Cdb', b), ('Cb', a)])        # -1/2 rho L^dag L
        for g, _ in m.gain:
            for a in np.nonzero(g)[0]:
                for b in np.nonzero(g)[0]:
                    add(1j * g[a] * np.conj(g[b]), [('Cdk', a), ('Cdb', b)])          # G rho G^dag
                    add(-0.5j * np.conj(g[a]) * g[b], [('Ck', a), ('Cdk', b)])        # -1/2 G^dag G rho
                    add(-0.5j * np.conj(g[a]) * g[b], [('Cb', b), ('Cdb', a)])        # -1/2 rho G^dag G


def _local_vec(site, n):
    v = np.zeros(4); v[site.state_index('empty_k empty_b')] = 1 - n; v[site.state_index('full_k full_b')] = n
    return v / np.linalg.norm(v)


def product_rho(sites, occupations):
    """|rho>> for rho = prod_j [(1-n_j)|0><0| + n_j|1><1|] (occupations indexed by MPS site)."""
    return MPS.from_product_state(sites, [_local_vec(s, n) for s, n in zip(sites, occupations)], bc='finite')


def identity_bra(sites):
    return product_rho(sites, [0.5] * len(sites))


class LindbladMPS:
    """Bundle model + state; `evolve` towards the NESS, `correlations`/`currents` read it out."""

    def __init__(self, model, pos=None, V=0.0, int_bonds=(), conserve=True, occupations=None):
        self.model = model
        self.pos = np.arange(model.N) if pos is None else np.asarray(pos)
        self.M = VectorizedLindbladian(dict(model=model, pos=self.pos, V=V, int_bonds=list(int_bonds), conserve=conserve))
        self.sites = self.M.lat.mps_sites()
        self.I = identity_bra(self.sites)
        if occupations is None:          # default: Gaussian (V=0) NESS occupations of the same model -> shortest relaxation
            from .ness import ness_covariance
            occupations = np.clip(np.real(np.diag(ness_covariance(model))), 1e-6, 1 - 1e-6)
        self.psi = product_rho(self.sites, np.asarray(occupations)[np.argsort(self.pos)])
        self.t = 0.0
        self.eng = None

    def evolve(self, T, dt=0.1, chi=64, **opts):
        o = dict(dt=dt, N_steps=max(1, int(round(T / dt))), order=2, approximation='II',
                 compression_method='variational', trunc_params=dict(chi_max=chi, svd_min=1e-12),
                 max_sweeps=2, preserve_norm=True)
        o.update(opts)
        if self.eng is None:
            self.eng = ExpMPOEvolution(self.psi, self.M, o)
        else:
            self.eng.options.update(o); self.eng.N_steps = o['N_steps']
        self.eng.run()
        self.psi = self.eng.psi
        self.t = float(self.eng.evolved_time)
        return self

    def trace(self):
        return self.I.overlap(self.psi)

    def expect(self, term):
        """<<I| prod op |rho>> / <<I|rho>>, term = [(ket-op name, mode), ...]."""
        env = MPSEnvironment(self.I, self.psi)
        return env.expectation_value_term([(op, int(self.pos[i])) for op, i in term]) / self.trace()

    def correlations(self, pairs):
        """C_ij = <c_i^dag c_j> for the requested (i, j) mode pairs; other entries left 0."""
        C = np.zeros((self.model.N, self.model.N), complex)
        for i, j in pairs:
            C[i, j] = self.expect([('Cdk', i), ('Ck', j)])
        return C

    def lead_pairs(self, labels=("L", "R")):
        pairs = set()
        for vecs in (self.model.loss, self.model.gain):
            for v, lab in vecs:
                if lab in labels:
                    k = np.nonzero(v)[0][0]
                    pairs.add((k, k))
                    for s in np.nonzero(self.model.h[k])[0]:
                        pairs.add((k, s)); pairs.add((s, k))
        return sorted(pairs)

    def currents(self, labels=("L", "R")):
        C = self.correlations(self.lead_pairs(labels))
        return {lab: dict(E=float(np.real(energy_current(self.model, C, lab))),
                          N=float(np.real(particle_current(self.model, C, lab)))) for lab in labels}


def lead_rate_terms(model, label):
    """{mode k: (loss rate, gain rate)} of the single-mode jump operators of reservoir `label`."""
    rates = {}
    for l, lab in model.loss:
        if lab == label:
            k = int(np.nonzero(l)[0][0]); rates.setdefault(k, [0.0, 0.0])[0] += abs(l[k]) ** 2
    for g, lab in model.gain:
        if lab == label:
            k = int(np.nonzero(g)[0][0]); rates.setdefault(k, [0.0, 0.0])[1] += abs(g[k]) ** 2
    return rates


def source_terms(base, target, labels=("L", "R")):
    """Terms of (L_target - L_base) for reservoirs whose only difference is the lead-mode rates:
    (a_k - b_k) D[c_k] + (g_k - h_k) D[c_k^dag] in the superfermion representation."""
    terms = []
    for lab in labels:
        rb, rt = lead_rate_terms(base, lab), lead_rate_terms(target, lab)
        for k in rt:
            dl, dg = rt[k][0] - rb[k][0], rt[k][1] - rb[k][1]
            if abs(dl) > 1e-15:
                terms += [(-dl, [('Ck', k), ('Cb', k)]), (-0.5 * dl, [('Cdk', k), ('Ck', k)]), (-0.5 * dl, [('Cdb', k), ('Cb', k)])]
            if abs(dg) > 1e-15:
                terms += [(dg, [('Cdk', k), ('Cdb', k)]), (-0.5 * dg, [('Ck', k), ('Cdk', k)]), (-0.5 * dg, [('Cb', k), ('Cdb', k)])]
    return terms


class LindbladResponse:
    """Exact finite-difference response without subtracting two NESS:
        rho_t - rho_b = int_0^inf exp(L_t s) (L_t - L_b) rho_b ds     (L_t rho_t = 0 = L_b rho_b).
    `base` is a converged LindbladMPS (state |rho_b>>), `target` the Model with the modified reservoirs.
    The perturbation |chi(s)>> = exp(L_t s)(L_t - L_b)|rho_b>> is traceless and decays; its lead currents are
    integrated in time (trapezoid + exponential tail) and divided by <<I|rho_b>>."""

    def __init__(self, base, target, chi=128, **evolve_opts):
        self.base, self.target = base, target
        self.T = LindbladMPS.__new__(LindbladMPS)
        self.T.model, self.T.pos = target, base.pos
        self.T.M = VectorizedLindbladian(dict(model=target, pos=base.pos, V=base.M.options['V'],
                                              int_bonds=base.M.options['int_bonds'], conserve=base.M.options.get('conserve', True)))
        self.T.sites, self.T.I, self.T.t, self.T.eng = base.sites, base.I, 0.0, None
        S = VectorizedLindbladian(dict(model=None, N=target.N, pos=base.pos, terms=source_terms(base.model, target)))
        psi = base.psi.copy()
        S.H_MPO.apply(psi, dict(compression_method='SVD', trunc_params=dict(chi_max=chi, svd_min=1e-14)))
        self.T.psi = psi
        self.trace_b = base.trace()
        self.chi, self.evolve_opts = chi, evolve_opts
        self.zero = {lab: energy_current(target, np.zeros((target.N, target.N), complex), lab) for lab in ("L", "R")}
        cb = base.currents(); self.E_b = {lab: cb[lab]["E"] for lab in ("L", "R")}
        self.ts, self.E = [0.0], [self._currents()]

    def _currents(self):
        """Linear functional E_lab[chi] = <<I|J_lab|chi>> per unit trace of rho_b (constant part removed since Tr chi = 0).
        A spurious stationary component (trace leakage tau = <<I|chi>> from truncation / W^II error) does not decay;
        it is removed as tau * E_b (error O(tau * dQ))."""
        C = np.zeros((self.target.N, self.target.N), complex)
        env = MPSEnvironment(self.T.I, self.T.psi)
        for i, j in self.T.lead_pairs(("L", "R")):
            C[i, j] = env.expectation_value_term([('Cdk', int(self.T.pos[i])), ('Ck', int(self.T.pos[j]))])
        tau = np.real(self.T.I.overlap(self.T.psi)) / np.real(self.trace_b)
        return {lab: float(np.real(energy_current(self.target, C, lab) - self.zero[lab]) / np.real(self.trace_b) - tau * (self.E_b[lab] - self.zero[lab]))
                for lab in ("L", "R")}

    def run(self, T, dt=0.1, chunk=1.0, t_fine=6.0):
        """Evolve to time T; sample every step while t < t_fine (fast lead transient), then every `chunk`."""
        while self.T.t < T - 1e-9:
            step = dt if self.T.t < t_fine - 1e-9 else chunk
            self.T.evolve(step, dt=dt, chi=self.chi, preserve_norm=False, max_trunc_err=np.inf,
                          compression_method='SVD', **self.evolve_opts)
            self.ts.append(self.T.t); self.E.append(self._currents())
        return self

    def integral(self, lab, tail=True):
        """int_0^inf E_lab(s) ds: trapezoid over the data, plus an exponential-tail extrapolation from the last 5 points."""
        t = np.array(self.ts); y = np.array([e[lab] for e in self.E])
        I = np.trapezoid(y, t)
        if tail and len(y) > 6 and np.all(np.sign(y[-5:]) == np.sign(y[-1])) and y[-1] != 0:
            p = np.polyfit(t[-5:], np.log(np.abs(y[-5:])), 1)
            if p[0] < 0:
                I += y[-1] / (-p[0])
        return float(I)
