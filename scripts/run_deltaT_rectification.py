"""Delta-T-response rectification from lossy Landauer-Buettiker transmissions; skin/winding diagnostics;
Lyapunov cross-checks (density, lead entropy production)."""
import json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from skinheat import chain_heff, transmissions, delta_heat, rectification, build_chain_with_leads, ness_covariance, lead_currents
from skinheat.model import fermi

Th, Tc, J, Gam = 1.0, 0.5, 1.0, 1.0

def xi_from_heff(heff, Gamma=Gam):
    H = heff.copy(); H[0, 0] += 0.5j * Gamma; H[-1, -1] += 0.5j * Gamma   # strip lead self-energies
    w, V = np.linalg.eig(H); P = np.abs(V) ** 2; P /= P.sum(0, keepdims=True); p = P.mean(1)
    x = np.arange(len(p)); s = np.polyfit(x, np.log(p + 1e-300), 1)[0]
    return -1 / s if abs(s) > 1e-12 else np.inf

def winding_E(heff, E, nk=4001):
    j = heff.shape[0] // 2
    d, tp, tm = heff[j, j], heff[j, j + 1], heff[j + 1, j]
    k = np.linspace(-np.pi, np.pi, nk); hk = d + tp * np.exp(1j * k) + tm * np.exp(-1j * k)
    a = np.unwrap(np.angle(hk - E)); return (a[-1] - a[0]) / (2 * np.pi)

def point(L, kappa, phi, gain_ratio=0.0):
    heff = chain_heff(L, J=J, kappa=kappa, phi=phi, gain_ratio=gain_ratio, Gamma=Gam)
    dQf, dQr, eps, TRL, TLR = delta_heat(heff, Th, Tc, Gam)
    return dict(dQf=dQf, dQr=dQr, R=rectification(dQf, dQr), xi=xi_from_heff(heff),
                lnratio=float(np.log(abs(dQf) / abs(dQr)))), (eps, TRL, TLR, heff)

res = {"params": dict(Th=Th, Tc=Tc, J=J, Gamma=Gam)}
phis = np.linspace(0, np.pi, 33)
res["phi"] = [dict(phi=float(p), **point(40, 0.3, p)[0]) for p in phis]
kaps = np.logspace(-2.5, 0.7, 20)
res["kappa"] = [dict(kappa=float(k), **point(40, k, np.pi / 2)[0]) for k in kaps]
Ls = [5, 10, 20, 30, 40, 60, 80, 120, 160]
res["L"] = {str(k): [dict(L=L, **point(L, k, np.pi / 2)[0]) for L in Ls] for k in [0.05, 0.1, 0.3]}
res["gain"] = [dict(gain_ratio=g, **point(40, 0.3, np.pi / 2, g)[0]) for g in [0, 0.25, 0.5, 1.0]]
# energy-resolved picture at phi=pi/2, kappa=0.3, L=40
p, (eps, TRL, TLR, heff) = point(40, 0.3, np.pi / 2)
wE = np.array([winding_E(heff, e) for e in eps[::20]])
res["energy_resolved"] = dict(eps=eps[::20].tolist(), TRL=TRL[::20].tolist(), TLR=TLR[::20].tolist(), winding=wE.tolist())
# Lyapunov cross-check: entropy production of the leads and chain density at a few phases
lyap = []
for phi in [0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4]:
    row = dict(phi=float(phi))
    for tag, (TL, TR) in {"f": (Th, Tc), "r": (Tc, Th), "iso": (Tc, Tc)}.items():
        m = build_chain_with_leads(40, kappa=0.3, phi=phi, K=80, TL=TL, TR=TR)
        C = ness_covariance(m); c = lead_currents(m, C)
        row[tag] = dict(sigma_leads=-(c["L"]["E"] / TL + c["R"]["E"] / TR), E_bulk=c["bulk"]["E"],
                        N_bulk=c["bulk"]["N"], n=float(np.real(np.trace(C[:40, :40])) / 40))
    lyap.append(row)
res["lyapunov_check"] = lyap
json.dump(res, open("results/deltaT_rectification.json", "w"), indent=1, default=float)

fig, ax = plt.subplots(2, 3, figsize=(14, 8))
P = res["phi"]; ax[0, 0].plot([d["phi"] for d in P], [d["R"] for d in P], "o-"); ax[0, 0].set_xlabel(r"$\phi$"); ax[0, 0].set_ylabel("R (ΔT response)"); ax[0, 0].set_title(r"L=40, $\kappa$=0.3")
a2 = ax[0, 0].twinx(); a2.plot([d["phi"] for d in P], [1 / d["xi"] if np.isfinite(d["xi"]) else 0 for d in P], "r--"); a2.set_ylabel(r"$1/\xi$", color="r")
Kp = res["kappa"]; ax[0, 1].semilogx([d["kappa"] for d in Kp], [d["R"] for d in Kp], "o-"); ax[0, 1].set_xlabel(r"$\kappa$"); ax[0, 1].set_ylabel("R"); ax[0, 1].set_title(r"$\phi=\pi/2$, L=40")
a3 = ax[0, 1].twinx(); a3.loglog([d["kappa"] for d in Kp], [abs(d["dQf"]) for d in Kp], "g--", label=r"$\delta Q_f$"); a3.loglog([d["kappa"] for d in Kp], [abs(d["dQr"]) for d in Kp], "m--", label=r"$\delta Q_r$"); a3.legend(fontsize=8)
for k, rows in res["L"].items():
    ax[0, 2].plot([d["L"] / abs(d["xi"]) for d in rows], [d["lnratio"] for d in rows], "o-", label=f"κ={k}, ξ={abs(rows[0]['xi']):.1f}")
ax[0, 2].set_xlabel(r"$L/|\xi|$"); ax[0, 2].set_ylabel(r"$\ln(\delta Q_f/\delta Q_r)$"); ax[0, 2].legend(); ax[0, 2].set_title("exponential nonreciprocity")
E = res["energy_resolved"]; ax[1, 0].semilogy(E["eps"], E["TRL"], label=r"$T_{RL}$ (L→R)"); ax[1, 0].semilogy(E["eps"], E["TLR"], label=r"$T_{LR}$ (R→L)"); ax[1, 0].set_xlabel(r"$\varepsilon$"); ax[1, 0].legend(); ax[1, 0].set_title(r"$\phi=\pi/2,\ \kappa=0.3,\ L=40$")
ax[1, 1].plot(E["eps"], np.log10(np.array(E["TRL"]) / np.array(E["TLR"])), label=r"$\log_{10} T_{RL}/T_{LR}$"); ax[1, 1].plot(E["eps"], E["winding"], "r--", label="winding w(ε)"); ax[1, 1].set_xlabel(r"$\varepsilon$"); ax[1, 1].legend()
G = res["gain"]; ax[1, 2].plot([d["gain_ratio"] for d in G], [d["R"] for d in G], "o-"); ax[1, 2].set_xlabel("bond gain / loss"); ax[1, 2].set_ylabel("R"); ax[1, 2].set_title("with bond gain (same chirality)")
plt.tight_layout(); plt.savefig("results/deltaT_rectification.png", dpi=130)

print("phi sweep:"); [print(f"  phi={d['phi']:.3f} R={d['R']:+.4f} dQf={d['dQf']:+.3e} dQr={d['dQr']:+.3e} xi={d['xi']:+.2f}") for d in P[::4]]
print("kappa sweep:"); [print(f"  k={d['kappa']:.4f} R={d['R']:+.4f} dQf={d['dQf']:+.3e} dQr={d['dQr']:+.3e} xi={d['xi']:+.2f}") for d in Kp[::2]]
print("L sweep ln(dQf/dQr) vs L/|xi|:")
for k, rows in res["L"].items(): print(" ", k, [(d["L"], round(d["L"] / abs(d["xi"]), 1), round(d["lnratio"], 2)) for d in rows])
print("gain:", [(d["gain_ratio"], round(d["R"], 4)) for d in G])
print("lyapunov check:"); [print(f"  phi={r['phi']:.3f} " + " ".join(f"{t}: sigma={r[t]['sigma_leads']:+.4f} Eb={r[t]['E_bulk']:+.4f} n={r[t]['n']:.3f}" for t in ["f", "r", "iso"])) for r in lyap]
