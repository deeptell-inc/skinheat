"""Thermal bulk bath (bond operator A_j coupled to a mesoscopic bath at (T_b, mu_b)) replaces the pure bond loss:
 - the Delta-T responses dQ_f, dQ_r are exactly independent of (T_b, mu_b) (linear response in the lead distributions),
 - entropy production is >= 0 and vanishes at global equilibrium,
 - the empty-bath limit mu_b -> -inf reproduces the pure-loss model (up to band effects of the mesoscopic bath)."""
import numpy as np
from skinheat import build_chain_with_leads, lead_currents, entropy_production

L, kappa, K = 12, 0.2, 8


def responses(bb):
    cur = {cfg: lead_currents(build_chain_with_leads(L, kappa=kappa, phi=np.pi / 2, K=K, TL=TL, TR=TR, bulk_bath=bb))
           for cfg, (TL, TR) in {"f": (1.0, 0.5), "r": (0.5, 1.0), "iso": (0.5, 0.5)}.items()}
    return -(cur["f"]["R"]["E"] - cur["iso"]["R"]["E"]), -(cur["r"]["L"]["E"] - cur["iso"]["L"]["E"]), cur


def test_thermal_bulk_bath():
    bb = lambda T, mu: dict(K=8, W=4.0, T=T, mu=mu)
    f1, r1, c1 = responses(bb(0.5, 0.0))
    f2, r2, c2 = responses(bb(2.0, -1.0))
    assert abs(f1 - f2) < 1e-9 and abs(r1 - r2) < 1e-9          # responses independent of the bath state
    f0, r0, _ = responses(None)
    assert abs(f1 / f0 - 1) < 0.15 and abs(np.log(f1 / r1) - np.log(f0 / r0)) < 0.15 * np.log(f0 / r0)
    # entropy production: >= 0 in all configurations, = 0 at global equilibrium (T_b = T_L = T_R, mu equal)
    for T, mu in [(0.5, 0.0), (2.0, -1.0), (0.25, 0.0)]:
        m = {cfg: build_chain_with_leads(L, kappa=kappa, phi=np.pi / 2, K=K, TL=TL, TR=TR, bulk_bath=bb(T, mu))
             for cfg, (TL, TR) in {"f": (1.0, 0.5), "r": (0.5, 1.0), "iso": (0.5, 0.5)}.items()}
        for cfg, mod in m.items():
            s = entropy_production(mod, lead_currents(mod))
            if cfg == "iso" and T == 0.5 and mu == 0.0:
                assert abs(s) < 1e-8, s
            else:
                assert s > 0, (cfg, T, mu, s)


if __name__ == "__main__":
    test_thermal_bulk_bath(); print("OK")


def test_sigma_mirror_transpose_identity():
    """sigma_f(T_b=T_c) = sigma_f(T_b=T_h) exactly: both equal (1/T_c-1/T_h) * int e (f_h-f_c) [1-R(e)] with the end
    reflection R_LL = R_RR, because mirror reflection maps the uniform chiral chain onto its transpose (P h P = h^T).
    Generic on-site disorder breaks P h P = h^T (and the identity); mirror-symmetric disorder restores it."""
    L, K = 16, 6
    rng = np.random.default_rng(3); dis = rng.uniform(-1, 1, L); sym = dis.copy(); sym[L // 2:] = dis[:L // 2][::-1]
    def sig(onsite, Tb, phi=np.pi / 3):
        m = build_chain_with_leads(L, kappa=0.2, phi=phi, K=K, TL=1.0, TR=0.5, onsite=onsite, bulk_bath=dict(K=6, W=4.0, T=Tb, mu=0.0))
        return entropy_production(m, lead_currents(m))
    assert abs(sig(0.0, 0.5) - sig(0.0, 1.0)) < 1e-9
    assert abs(sig(sym, 0.5) - sig(sym, 1.0)) < 1e-9
    assert abs(sig(dis, 0.5) - sig(dis, 1.0)) > 1e-4
