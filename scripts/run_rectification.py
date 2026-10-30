"""Sweep chirality phase phi, bulk dissipation kappa and size L; compute heat rectification,
skin localization length, spectral winding, lead entropy production and bulk energy drain."""
import json, sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from skinheat import build_chain_with_leads, ness_covariance, lead_currents, localization_length, winding_number

Th, Tc = 1.0, 0.5

def run(L, kappa, phi, gain_ratio=0.0, K=60):
    out = {}
    for tag, (TL, TR) in {"f": (Th, Tc), "r": (Tc, Th)}.items():
        m = build_chain_with_leads(L, kappa=kappa, phi=phi, gain_ratio=gain_ratio, K=K, TL=TL, TR=TR)
        C = ness_covariance(m)
        cur = lead_currents(m, C)
        cold = "R" if tag == "f" else "L"
        hot = "L" if tag == "f" else "R"
        out[tag] = dict(Q_cold=-cur[cold]["E"], Q_hot=cur[hot]["E"], E_bulk=cur["bulk"]["E"],
                        N_bulk=cur["bulk"]["N"],
                        sigma_leads=-(cur["L"]["E"] / TL + cur["R"]["E"] / TR),
                        n_chain=float(np.real(np.trace(C[np.ix_(m.chain_sites, m.chain_sites)])) / L))
        if tag == "f":
            out["xi"] = localization_length(m); out["w"] = winding_number(m)
    Qf, Qr = out["f"]["Q_cold"], out["r"]["Q_cold"]
    out["R"] = (abs(Qf) - abs(Qr)) / (abs(Qf) + abs(Qr) + 1e-300)
    return out

res = {"params": dict(Th=Th, Tc=Tc, J=1.0, Gamma=1.0, W=4.0, K=60)}
# 1) phi sweep
phis = np.linspace(0, np.pi, 25)
res["phi_sweep"] = [dict(phi=float(p), **{k: v for k, v in run(40, 0.3, p).items()}) for p in phis]
# 2) kappa sweep at phi=pi/2
kaps = np.concatenate([[0.0], np.logspace(-2, 0.5, 16)])
res["kappa_sweep"] = [dict(kappa=float(k), **run(40, k, np.pi / 2)) for k in kaps]
# 3) size dependence
res["L_sweep"] = [dict(L=L, **run(L, 0.3, np.pi / 2)) for L in [10, 20, 40, 80, 160]]
# 4) with bond gain (keeps density up)
res["gain_sweep"] = [dict(gain_ratio=g, **run(40, 0.3, np.pi / 2, gain_ratio=g)) for g in [0.0, 0.25, 0.5, 0.75, 1.0]]
json.dump(res, open("results/rectification.json", "w"), indent=1, default=float)

fig, ax = plt.subplots(2, 2, figsize=(10, 7))
ps = res["phi_sweep"]
ax[0, 0].plot([d["phi"] for d in ps], [d["R"] for d in ps], "o-"); ax[0, 0].set_xlabel(r"$\phi$"); ax[0, 0].set_ylabel("rectification R"); ax[0, 0].set_title(r"L=40, $\kappa$=0.3")
ax2 = ax[0, 0].twinx(); ax2.plot([d["phi"] for d in ps], [d["w"] for d in ps], "r--"); ax2.set_ylabel("winding", color="r")
ks = res["kappa_sweep"]
ax[0, 1].semilogx([d["kappa"] + 1e-3 for d in ks], [d["R"] for d in ks], "o-"); ax[0, 1].set_xlabel(r"$\kappa$"); ax[0, 1].set_ylabel("R"); ax[0, 1].set_title(r"$\phi=\pi/2$")
ax3 = ax[0, 1].twinx(); ax3.semilogx([d["kappa"] + 1e-3 for d in ks], [abs(d["f"]["Q_cold"]) for d in ks], "g--", label="|Q_f|"); ax3.semilogx([d["kappa"] + 1e-3 for d in ks], [abs(d["r"]["Q_cold"]) for d in ks], "m--", label="|Q_r|"); ax3.legend(); ax3.set_ylabel("|Q_cold|")
Ls = res["L_sweep"]
ax[1, 0].plot([d["L"] for d in Ls], [d["R"] for d in Ls], "o-"); ax[1, 0].set_xlabel("L"); ax[1, 0].set_ylabel("R")
ax[1, 1].plot([d["phi"] for d in ps], [d["f"]["sigma_leads"] for d in ps], "o-", label="forward"); ax[1, 1].plot([d["phi"] for d in ps], [d["r"]["sigma_leads"] for d in ps], "s-", label="reverse"); ax[1, 1].set_xlabel(r"$\phi$"); ax[1, 1].set_ylabel(r"$\sigma_{leads}$"); ax[1, 1].legend()
plt.tight_layout(); plt.savefig("results/rectification.png", dpi=130)

print("phi sweep (phi, R, Qf, Qr, xi, w, E_bulk_f):")
for d in ps[::3]:
    print(f"  {d['phi']:.3f}  R={d['R']:+.4f}  Qf={d['f']['Q_cold']:+.5f}  Qr={d['r']['Q_cold']:+.5f}  xi={d['xi']:+.2f}  w={d['w']:+.2f}  Ebulk={d['f']['E_bulk']:+.4f}")
print("kappa sweep (kappa, R, Qf, Qr, xi):")
for d in ks[::2]:
    print(f"  {d['kappa']:.3f}  R={d['R']:+.4f}  Qf={d['f']['Q_cold']:+.5f}  Qr={d['r']['Q_cold']:+.5f}  xi={d['xi']:+.2f}  n={d['f']['n_chain']:.3f}")
print("L sweep:", [(d["L"], round(d["R"], 4), round(d["xi"], 2)) for d in Ls])
print("gain sweep:", [(d["gain_ratio"], round(d["R"], 4), round(d["f"]["n_chain"], 3)) for d in res["gain_sweep"]])
