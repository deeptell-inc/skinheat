"""Stage 3: finite-temperature bulk bath (detailed balance on the bond operator) and the full entropy balance.
Sweeps kappa x T_b (mu_b=0), mu_b (sink <-> thermal), and the bath discretisation (K_b, W_b) for the exact law.
Writes results/thermal_bulk.json / .png incrementally."""
import json, time, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from skinheat import build_chain_with_leads, ness_covariance, lead_currents, entropy_production

L, K, J, Th, Tc, phi = 40, 10, 1.0, 1.0, 0.5, np.pi / 2
CFG = {"f": (Th, Tc), "r": (Tc, Th), "iso": (Tc, Tc)}
OUT = "results/thermal_bulk.json"
res = dict(params=dict(L=L, K=K, J=J, Th=Th, Tc=Tc, phi=phi, W_lead=4.0))


def point(kappa, bb, profile=False):
    row = dict(kappa=kappa, bulk_bath=bb)
    for cfg, (TL, TR) in CFG.items():
        m = build_chain_with_leads(L, kappa=kappa, phi=phi, K=K, TL=TL, TR=TR, bulk_bath=bb)
        C = ness_covariance(m); c = lead_currents(m, C)
        row[cfg] = dict(E={k: v["E"] for k, v in c.items()}, N={k: v["N"] for k, v in c.items()},
                        sigma=entropy_production(m, c), n_mean=float(np.real(np.trace(C[:L, :L])) / L))
        if profile:
            row[cfg]["n_profile"] = np.real(np.diag(C)[:L]).tolist()
    dQf = -(row["f"]["E"]["R"] - row["iso"]["E"]["R"]); dQr = -(row["r"]["E"]["L"] - row["iso"]["E"]["L"])
    tp, tm = -J - 0.5j * kappa * np.exp(1j * phi), -J - 0.5j * kappa * np.exp(-1j * phi)
    row.update(dQf=dQf, dQr=dQr, lnratio=float(np.log(dQf / dQr)), R=(abs(dQf) - abs(dQr)) / (abs(dQf) + abs(dQr)),
               nn_law=float(2 * (L - 1) * np.log(abs(tm / tp))),
               dsigma_f=row["f"]["sigma"] - row["iso"]["sigma"], dsigma_r=row["r"]["sigma"] - row["iso"]["sigma"],
               Qbulk_f=-row["f"]["E"]["bulk"] if bb else None)
    return row


def save():
    json.dump(res, open(OUT, "w"), indent=1, default=float)


t0 = time.time()
kaps = [0.05, 0.1, 0.2, 0.3, 0.5]
Tbs = [0.25, 0.5, 1.0, 2.0]
res["kappa_Tb"] = [point(k, dict(K=10, W=4.0, T=Tb, mu=0.0), profile=(Tb in (0.25, 2.0))) for k in kaps for Tb in Tbs]
save(); print("kappa x Tb done", time.time() - t0, flush=True)
res["pure_loss"] = [point(k, None) for k in kaps]
save()
res["mu_sweep"] = [point(0.2, dict(K=10, W=4.0, T=0.5, mu=mu), profile=True) for mu in [0.0, -0.5, -1.0, -2.0, -4.0, -30.0]]
save(); print("mu sweep done", time.time() - t0, flush=True)
res["discretisation"] = [point(0.2, dict(K=Kb, W=Wb, T=0.5, mu=0.0)) for Kb, Wb in [(10, 4.0), (10, 6.0), (20, 4.0), (20, 6.0), (40, 6.0)]]
save(); print("discretisation done", time.time() - t0, flush=True)

# ---- figure ----
fig, ax = plt.subplots(2, 3, figsize=(16, 9))
KT = res["kappa_Tb"]
for Tb, m in zip(Tbs, "osd^"):
    rows = [r for r in KT if r["bulk_bath"]["T"] == Tb]
    ax[0, 0].semilogx([r["kappa"] for r in rows], [r["lnratio"] for r in rows], m + "-", label=rf"$T_b$={Tb}")
    ax[0, 1].semilogx([r["kappa"] for r in rows], [r["f"]["sigma"] for r in rows], m + "-", label=rf"$\sigma_f$, $T_b$={Tb}")
    ax[0, 1].semilogx([r["kappa"] for r in rows], [r["iso"]["sigma"] for r in rows], m + ":", mfc="none", label=rf"$\sigma_{{iso}}$, $T_b$={Tb}")
    ax[0, 2].semilogx([r["kappa"] for r in rows], [r["dsigma_f"] / r["dQf"] for r in rows], m + "-", label=rf"$T_b$={Tb}")
    ax[1, 0].semilogx([r["kappa"] for r in rows], [r["Qbulk_f"] for r in rows], m + "-", label=rf"$T_b$={Tb}")
PL = res["pure_loss"]
ax[0, 0].semilogx([r["kappa"] for r in PL], [r["lnratio"] for r in PL], "k--", label="pure loss (sink)")
ax[0, 0].semilogx([r["kappa"] for r in PL], [r["nn_law"] for r in PL], "k:", label=r"$2(L-1)\ln|t_-/t_+|$")
ax[0, 0].set_xlabel(r"$\kappa$"); ax[0, 0].set_ylabel(r"$\ln(\delta Q_f/\delta Q_r)$"); ax[0, 0].legend(fontsize=8); ax[0, 0].set_title(r"(a) nonreciprocity is independent of $T_b$ (exactly)")
ax[0, 1].set_xlabel(r"$\kappa$"); ax[0, 1].set_ylabel(r"entropy production $\sigma$"); ax[0, 1].legend(fontsize=7); ax[0, 1].set_title(r"(b) $\sigma=-\sum_\alpha Q_\alpha/T_\alpha\geq 0$, all reservoirs thermal")
ax[0, 2].set_xlabel(r"$\kappa$"); ax[0, 2].set_ylabel(r"$(\sigma_f-\sigma_{iso})/\delta Q_f$"); ax[0, 2].legend(fontsize=8); ax[0, 2].set_title(r"(c) entropy cost per unit rectified heat")
ax[1, 0].set_xlabel(r"$\kappa$"); ax[1, 0].set_ylabel(r"heat into bulk bath, forward"); ax[1, 0].legend(fontsize=8); ax[1, 0].set_title(r"(d) bulk heat flux $-E_{bulk}$ (forward)")
MU = res["mu_sweep"]
for r in MU:
    ax[1, 1].plot(r["f"]["n_profile"], "-", label=rf"$\mu_b$={r['bulk_bath']['mu']}")
ax[1, 1].set_xlabel("site j"); ax[1, 1].set_ylabel(r"$\langle n_j\rangle$ (forward)"); ax[1, 1].legend(fontsize=7); ax[1, 1].set_title(r"(e) density profile, $\kappa$=0.2, $T_b$=0.5: sink $\to$ thermal")
D = res["discretisation"]
ax[1, 2].bar([f"K={r['bulk_bath']['K']}\nW={r['bulk_bath']['W']}" for r in D], [r["lnratio"] for r in D])
ax[1, 2].axhline(D[0]["nn_law"], color="k", ls=":", label="NN law (wide band)"); ax[1, 2].axhline(PL[2]["lnratio"], color="k", ls="--", label="Lindblad loss, K_lead=10")
ax[1, 2].set_ylabel(r"$\ln(\delta Q_f/\delta Q_r)$"); ax[1, 2].legend(fontsize=8); ax[1, 2].set_title(r"(f) bath discretisation, $\kappa$=0.2")
plt.tight_layout(); plt.savefig("results/thermal_bulk.png", dpi=130)
for r in KT: print(f"kappa={r['kappa']} Tb={r['bulk_bath']['T']}: ln={r['lnratio']:.3f} R={r['R']:.6f} dQf={r['dQf']:.3e} sig_f={r['f']['sigma']:.4f} sig_iso={r['iso']['sigma']:.4f} dsig_f/dQf={r['dsigma_f']/r['dQf']:.3f} Qbulk_f={r['Qbulk_f']:+.4f} n={r['f']['n_mean']:.3f}")
for r in MU: print(f"mu={r['bulk_bath']['mu']}: ln={r['lnratio']:.3f} sig_f={r['f']['sigma']:.4f} n={r['f']['n_mean']:.3f} Nbulk_f={r['f']['N']['bulk']:+.4f}")
for r in D: print(f"Kb={r['bulk_bath']['K']} Wb={r['bulk_bath']['W']}: ln={r['lnratio']:.3f} (NN law {r['nn_law']:.3f})")
