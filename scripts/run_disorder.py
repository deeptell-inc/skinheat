"""Stage 4: robustness of skin rectification to disorder (Gaussian, wide-band leads).
 (a) on-site disorder W_d: ratio exactly invariant; throughput dQf collapses (Anderson)  -> rectification window
 (b) L-dependence of <ln dQf> for several W_d: effective forward decay length vs skin / localisation lengths
 (c) bond disorder (kappa_j, phi_j) and random chirality: ln ratio = 2 sum_j ln|t-_j/t+_j| (sum rule, self-averaging)
 (d) NNN t2 + on-site disorder: ratio becomes sample dependent
Writes results/disorder.json / .png incrementally."""
import json, time, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from skinheat import chain_heff, delta_heat
from skinheat.transmission import chain_heff_bonds

J, Th, Tc, NE = 1.0, 1.0, 0.5, 2001
OUT = "results/disorder.json"
rng = np.random.default_rng(42)
res = dict(params=dict(J=J, Th=Th, Tc=Tc, Gamma=1.0, W=4.0, ne=NE, phi=np.pi / 2))
save = lambda: json.dump(res, open(OUT, "w"), indent=1, default=float)
lnlaw = lambda kappa, L: 2 * (L - 1) * np.log(abs((-1 + 0.5j * kappa * 1j) / (-1 - 0.5j * kappa * 1j)))


def sample(L, kappa, Wd, t2=0.0, ns=100):
    out = []
    for _ in range(ns):
        f, r, *_ = delta_heat(chain_heff(L, kappa=kappa, phi=np.pi / 2, t2=t2, onsite=rng.uniform(-Wd / 2, Wd / 2, L)), Th, Tc, ne=NE)
        out.append((f, r))
    f, r = np.array(out).T
    return dict(lnQf_mean=float(np.mean(np.log(f))), lnQf_std=float(np.std(np.log(f))), Qf_median=float(np.median(f)),
                Qf_p10=float(np.percentile(f, 10)), Qf_p90=float(np.percentile(f, 90)),
                ln_mean=float(np.mean(np.log(f / r))), ln_std=float(np.std(np.log(f / r))), ln_max_dev=float(np.max(np.abs(np.log(f / r) - lnlaw(kappa, L)))),
                R_mean=float(np.mean((f - r) / (f + r))))


t0 = time.time()
# (a) on-site disorder, L=40
Wds = [0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0]
res["onsite"] = [dict(kappa=k, Wd=Wd, L=40, law=lnlaw(k, 40), **sample(40, k, Wd)) for k in [0.05, 0.1, 0.2] for Wd in Wds]
save(); print("onsite done", time.time() - t0, flush=True)
# (b) L dependence of the forward throughput
res["Ldep"] = [dict(kappa=k, Wd=Wd, L=L, **sample(L, k, Wd, ns=60)) for k in [0.0, 0.05, 0.2] for Wd in [0.5, 1.0, 2.0] for L in [10, 20, 40, 80]]
save(); print("Ldep done", time.time() - t0, flush=True)
# (c) bond disorder: kappa_j spread, phi_j spread, random chirality signs
bond = []
for L in [20, 40, 80]:
    for kind in ["kappa", "phi", "sign"]:
        lns, preds = [], []
        for _ in range(100):
            kap = np.full(L - 1, 0.1); phi = np.full(L - 1, np.pi / 2)
            if kind == "kappa": kap = 0.1 * (1 + 0.8 * rng.uniform(-1, 1, L - 1))
            if kind == "phi": phi = np.pi / 2 + 1.0 * rng.uniform(-1, 1, L - 1)
            if kind == "sign": phi = rng.choice([-1, 1], L - 1) * np.pi / 2
            f, r, *_ = delta_heat(chain_heff_bonds(L, kap, phi, onsite=rng.uniform(-0.5, 0.5, L)), Th, Tc, ne=NE)
            lns.append(np.log(f / r)); preds.append(2 * sum(np.log(abs((-1 - 0.5j * k * np.exp(-1j * p)) / (-1 - 0.5j * k * np.exp(1j * p)))) for k, p in zip(kap, phi)))
        bond.append(dict(L=L, kind=kind, ln_mean=float(np.mean(lns)), ln_std=float(np.std(lns)), max_dev_from_sum_rule=float(np.max(np.abs(np.array(lns) - preds))), clean=lnlaw(0.1, L)))
res["bond"] = bond; save(); print("bond done", time.time() - t0, flush=True)
# (d) NNN + on-site disorder
res["nnn"] = [dict(kappa=0.2, t2=t2, Wd=Wd, L=40, **sample(40, 0.2, Wd, t2=t2, ns=60)) for t2 in [0.0, 0.3] for Wd in [0, 0.5, 1.0, 2.0]]
save(); print("nnn done", time.time() - t0, flush=True)

# ---- figure ----
fig, ax = plt.subplots(2, 3, figsize=(16, 9))
for k, m in zip([0.05, 0.1, 0.2], "osd"):
    rows = [r for r in res["onsite"] if r["kappa"] == k]
    ax[0, 0].semilogy([r["Wd"] for r in rows], [r["Qf_median"] for r in rows], m + "-", label=rf"$\kappa$={k}")
    ax[0, 0].fill_between([r["Wd"] for r in rows], [r["Qf_p10"] for r in rows], [r["Qf_p90"] for r in rows], alpha=0.15)
    ax[0, 1].plot([r["Wd"] for r in rows], [r["ln_max_dev"] for r in rows], m + "-", label=rf"$\kappa$={k}")
ax[0, 0].set_xlabel("$W_d$"); ax[0, 0].set_ylabel(r"$\delta Q_f$ (median, 10–90 %)"); ax[0, 0].legend(); ax[0, 0].set_title("(a) throughput collapses with on-site disorder (L=40)")
ax[0, 1].set_yscale("log"); ax[0, 1].set_xlabel("$W_d$"); ax[0, 1].set_ylabel(r"max$_s$ |ln r − 2(L−1)ln|t₋/t₊||"); ax[0, 1].legend(); ax[0, 1].set_title("(b) ratio exactly invariant (100 samples)")
for Wd, c in zip([0.5, 1.0, 2.0], ["C0", "C1", "C2"]):
    for k, ls, m in [(0.0, ":", "x"), (0.05, "--", "o"), (0.2, "-", "s")]:
        rows = [r for r in res["Ldep"] if r["kappa"] == k and r["Wd"] == Wd]
        ax[0, 2].plot([r["L"] for r in rows], [r["lnQf_mean"] for r in rows], ls, marker=m, color=c, label=rf"$W_d$={Wd}, $\kappa$={k}")
ax[0, 2].set_xlabel("L"); ax[0, 2].set_ylabel(r"$\langle\ln\delta Q_f\rangle$"); ax[0, 2].legend(fontsize=7, ncol=2); ax[0, 2].set_title("(c) forward decay vs L: skin effect vs Anderson")
for kind, m in zip(["kappa", "phi", "sign"], "osd"):
    rows = [r for r in res["bond"] if r["kind"] == kind]
    ax[1, 0].errorbar([r["L"] for r in rows], [r["ln_mean"] for r in rows], [r["ln_std"] for r in rows], fmt=m + "-", capsize=3, label=f"bond disorder: {kind}")
ax[1, 0].plot([20, 40, 80], [lnlaw(0.1, L) for L in [20, 40, 80]], "k:", label="clean law")
ax[1, 0].set_xlabel("L"); ax[1, 0].set_ylabel(r"$\ln(\delta Q_f/\delta Q_r)$ (mean ± std)"); ax[1, 0].legend(fontsize=8); ax[1, 0].set_title(r"(d) bond disorder: $\ln r = 2\sum_j\ln|t_{-,j}/t_{+,j}|$ exactly")
for t2, m in zip([0.0, 0.3], "os"):
    rows = [r for r in res["nnn"] if r["t2"] == t2]
    ax[1, 1].errorbar([r["Wd"] for r in rows], [r["ln_mean"] for r in rows], [r["ln_std"] for r in rows], fmt=m + "-", capsize=3, label=rf"$t_2$={t2}")
ax[1, 1].set_xlabel("$W_d$"); ax[1, 1].set_ylabel(r"$\ln(\delta Q_f/\delta Q_r)$"); ax[1, 1].legend(); ax[1, 1].set_title(r"(e) NNN hopping: ratio becomes sample dependent ($\kappa$=0.2)")
for k, m in zip([0.05, 0.1, 0.2], "osd"):
    rows = [r for r in res["onsite"] if r["kappa"] == k]
    ax[1, 2].plot([r["Wd"] for r in rows], [r["R_mean"] for r in rows], m + "-", label=rf"$\kappa$={k}")
ax[1, 2].set_xlabel("$W_d$"); ax[1, 2].set_ylabel(r"$\langle R\rangle$"); ax[1, 2].legend(); ax[1, 2].set_title("(f) rectification coefficient (sample mean)")
plt.tight_layout(); plt.savefig("results/disorder.png", dpi=130)
for r in res["onsite"]: print(f"kappa={r['kappa']} Wd={r['Wd']}: Qf med={r['Qf_median']:.3e} [{r['Qf_p10']:.1e},{r['Qf_p90']:.1e}] ln dev max={r['ln_max_dev']:.1e} R={r['R_mean']:.5f}")
for r in res["Ldep"]: print(f"Ldep kappa={r['kappa']} Wd={r['Wd']} L={r['L']}: <ln Qf>={r['lnQf_mean']:.3f} ± {r['lnQf_std']:.3f}")
for r in res["bond"]: print(f"bond {r['kind']} L={r['L']}: ln={r['ln_mean']:.3f} ± {r['ln_std']:.3f} (clean {r['clean']:.3f}) sum-rule dev {r['max_dev_from_sum_rule']:.1e}")
for r in res["nnn"]: print(f"nnn t2={r['t2']} Wd={r['Wd']}: ln={r['ln_mean']:.3f} ± {r['ln_std']:.3f}  Qf med={r['Qf_median']:.3e}")
