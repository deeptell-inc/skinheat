"""Compression-method study (L=6, K=4, V=0 vs Gauss): is the ~1e-4 floor from the variational compression?"""
import sys, json, time, warnings, numpy as np
warnings.filterwarnings("ignore")
from skinheat import build_chain_with_leads, lead_currents
from skinheat.lindblad_mps import LindbladMPS, layout
tag = sys.argv[1]
opts = {"svd": dict(compression_method="SVD"),
        "var8": dict(compression_method="variational", max_sweeps=8, min_sweeps=4, tol_theta_diff=1e-12),
        "var4": dict(compression_method="variational", max_sweeps=4, min_sweeps=2, tol_theta_diff=1e-12),
        "zip": dict(compression_method="zip_up")}[tag]
L, K = 6, 4
m = build_chain_with_leads(L, kappa=0.2, phi=np.pi / 2, K=K, TL=1.0, TR=0.5)
g = lead_currents(m)
S = LindbladMPS(m, pos=layout(L, K)); t0 = time.time(); rows = []
for _ in range(16):
    S.evolve(2.0, dt=0.1, chi=64, max_trunc_err=np.inf, **opts); c = S.currents()
    rows.append(dict(t=S.t, errL=c["L"]["E"] - g["L"]["E"], errR=c["R"]["E"] - g["R"]["E"], wall=time.time() - t0))
    json.dump(dict(tag=tag, rows=rows), open(f"results/beyond_nn_runs/conv/{tag}.json", "w"), indent=1)
    print(f"{tag} t={S.t:.0f} errL={rows[-1]['errL']:+.2e} errR={rows[-1]['errR']:+.2e} wall={rows[-1]['wall']:.0f}", flush=True)
