"""Precision study of the W^II/MPS NESS vs Gauss at V=0: L=6, K=4 (N=14 sites), sweep chi and dt."""
import sys, json, time, warnings, numpy as np
warnings.filterwarnings("ignore")
from skinheat import build_chain_with_leads, lead_currents
from skinheat.lindblad_mps import LindbladMPS, layout
chi, dt = int(sys.argv[1]), float(sys.argv[2])
L, K = 6, 4
m = build_chain_with_leads(L, kappa=0.2, phi=np.pi / 2, K=K, TL=1.0, TR=0.5)
g = lead_currents(m)
S = LindbladMPS(m, pos=layout(L, K)); t0 = time.time(); rows = []
for _ in range(20):
    S.evolve(2.0, dt=dt, chi=chi, max_trunc_err=np.inf); c = S.currents()
    rows.append(dict(t=S.t, EL=c["L"]["E"], ER=c["R"]["E"], errL=c["L"]["E"] - g["L"]["E"], errR=c["R"]["E"] - g["R"]["E"], chi=int(max(S.psi.chi)), wall=time.time() - t0))
    json.dump(dict(chi=chi, dt=dt, gauss=g, rows=rows), open(f"results/beyond_nn_runs/conv/chi{chi}_dt{dt}.json", "w"), indent=1)
    print(f"chi={chi} dt={dt} t={S.t:.0f} errL={rows[-1]['errL']:+.2e} errR={rows[-1]['errR']:+.2e} wall={rows[-1]['wall']:.0f}", flush=True)
