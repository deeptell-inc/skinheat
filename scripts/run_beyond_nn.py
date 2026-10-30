"""Stage 2a (Gaussian): next-nearest-neighbour hopping t2 breaks the tridiagonal identity.
Sweeps t2 in {0,0.1,0.3,0.5} at phi=pi/2, kappa in {0.1,0.2,0.3}, L=40; records energy-resolved
T_RL/T_LR, the non-Bloch prediction 2(L-1) ln|beta_2 beta_3|(eps), and spectral winding numbers.
Writes results/beyond_nn.json (key "nnn") and results/beyond_nn.png (also plots key "interacting"
if scripts/run_beyond_nn_tebd.py has added it)."""
import json, os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from skinheat import chain_heff, delta_heat, rectification

J, L, phi, Th, Tc = 1.0, 40, np.pi / 2, 1.0, 0.5
OUT = "results/beyond_nn.json"


def tpm(kappa):
    return -J - 0.5j * kappa * np.exp(1j * phi), -J - 0.5j * kappa * np.exp(-1j * phi)


def nonbloch_roots(eps, kappa, t2):
    """Roots of t2 b^2 + t+ b + (d-eps) + t- b^-1 + t2 b^-2 = 0 (d=-i kappa), sorted by |b|."""
    tp, tm = tpm(kappa); d = -1j * kappa
    r = np.roots([t2, tp, d - eps, tm, t2] if t2 else [tp, d - eps, tm])
    return r[np.argsort(abs(r))]


def winding(eps, kappa, t2, shift=0.0, nk=8001):
    """Winding of the PBC ring h(k) = d + t+ e^{ik} + t- e^{-ik} + t2(e^{2ik}+e^{-2ik}) around eps - i*shift."""
    tp, tm = tpm(kappa); d = -1j * kappa
    k = np.linspace(-np.pi, np.pi, nk)
    hk = d + tp * np.exp(1j * k) + tm * np.exp(-1j * k) + t2 * (np.exp(2j * k) + np.exp(-2j * k))
    a = np.unwrap(np.angle(hk - (eps - 1j * shift)))
    return float(np.round((a[-1] - a[0]) / (2 * np.pi), 6))


def point(kappa, t2):
    tp, tm = tpm(kappa)
    heff = chain_heff(L, J=J, kappa=kappa, phi=phi, t2=t2)
    dQf, dQr, eps, TRL, TLR = delta_heat(heff, Th, Tc)
    lnT = np.log(TRL / TLR)
    sl = slice(0, len(eps), 25)
    e = eps[sl]
    p = 2 if t2 else 1
    nb = [2 * (L - 1) * np.log(abs(np.prod(nonbloch_roots(x, kappa, t2)[p - 1:p + 1]))) for x in e]
    return dict(kappa=kappa, t2=t2, dQf=dQf, dQr=dQr, R=rectification(dQf, dQr),
                lnratio=float(np.log(dQf / dQr)), pred_nn=float(2 * (L - 1) * np.log(abs(tm / tp))),
                lnT_min=float(lnT.min()), lnT_max=float(lnT.max()),
                eps=e.tolist(), lnT=lnT[sl].tolist(), nonbloch=nb,
                w_real=[winding(x, kappa, t2) for x in e], w_shift=[winding(x, kappa, t2, shift=kappa) for x in e])


res = json.load(open(OUT)) if os.path.exists(OUT) else {}
res["nnn"] = dict(params=dict(J=J, L=L, phi=phi, Th=Th, Tc=Tc, Gamma=1.0, W=4.0),
                  sweep=[point(k, t2) for k in [0.1, 0.2, 0.3] for t2 in [0.0, 0.1, 0.3, 0.5]])
json.dump(res, open(OUT, "w"), indent=1, default=float)

# ---- figure ----
S = res["nnn"]["sweep"]
fig, ax = plt.subplots(3, 3, figsize=(16, 13))
for k, m in zip([0.1, 0.2, 0.3], "osd"):
    rows = [r for r in S if r["kappa"] == k]
    ax[0, 0].plot([r["t2"] for r in rows], [r["lnratio"] for r in rows], m + "-", label=rf"$\kappa$={k}")
    ax[0, 0].plot([r["t2"] for r in rows], [r["pred_nn"] for r in rows], "k:", lw=0.8)
    ax[0, 1].semilogy([r["t2"] for r in rows], [r["dQf"] for r in rows], m + "-", label=rf"$\delta Q_f$, $\kappa$={k}")
    ax[0, 1].semilogy([r["t2"] for r in rows], [r["dQr"] for r in rows], m + "--", mfc="none", label=rf"$\delta Q_r$, $\kappa$={k}")
ax[0, 0].set_xlabel("$t_2$"); ax[0, 0].set_ylabel(r"$\ln(\delta Q_f/\delta Q_r)$"); ax[0, 0].legend(); ax[0, 0].set_title(r"(a) NNN hopping; dotted: NN law $2(L-1)\ln|t_-/t_+|$")
ax[0, 1].set_xlabel("$t_2$"); ax[0, 1].legend(fontsize=7); ax[0, 1].set_title(r"(b) $\Delta T$ response, $\phi=\pi/2$, L=40 (Gauss)")
for r in [r for r in S if r["kappa"] == 0.2]:
    l, = ax[0, 2].plot(r["eps"], r["lnT"], label=rf"$t_2$={r['t2']}")
    ax[0, 2].plot(r["eps"], r["nonbloch"], ":", color=l.get_color())
    ax[1, 0].plot(r["eps"], r["w_shift"], "-", color=l.get_color(), label=rf"$t_2$={r['t2']}")
ax[0, 2].set_ylim(-5, 40); ax[0, 2].set_xlabel(r"$\varepsilon$"); ax[0, 2].set_ylabel(r"$\ln T_{RL}/T_{LR}$"); ax[0, 2].legend(); ax[0, 2].set_title(r"(c) $\kappa$=0.2; dotted: $2(L-1)\ln|\beta_2\beta_3|(\varepsilon)$")
ax[1, 0].set_xlabel(r"$\varepsilon$"); ax[1, 0].set_ylabel(r"winding of $h(k)+i\kappa$ around $\varepsilon$"); ax[1, 0].legend(); ax[1, 0].set_title(r"(d) $w(\varepsilon)$ on the real axis is 0 for all $t_2$")

def rows_of(key):
    return {Lk: [r for r in rows if r.get("dQf") is not None] for Lk, rows in sorted(res[key]["sweep"].items(), key=lambda kv: int(kv[0]))} if key in res else {}

I = rows_of("interacting")      # kappa = 0.2 (reverse response below the MPS floor -> upper bounds)
for Lk, ok in I.items():
    V = [r["V"] for r in ok]
    ax[1, 1].errorbar(V, [r["dQf"] for r in ok], [r["dQf_err"] for r in ok], fmt="o-", capsize=3, label=f"L={Lk}")
    if "gauss" in ok[0]: ax[1, 1].plot([0], [ok[0]["gauss"]["dQf"]], "k*", ms=10, label=f"Gauss V=0, L={Lk}" if Lk == list(I)[0] else None)
    ax[1, 2].errorbar(V, [abs(r["dQr"]) for r in ok], [r["dQr_err"] for r in ok], fmt="s", capsize=3, label=rf"$|\delta Q_r|$ L={Lk} (MPS, $\chi$=64)")
    ax[1, 2].plot(V, [r["dQf"] * np.exp(-res["interacting"]["gauss_lnratio"][Lk]) for r in ok], ":", label=rf"$\delta Q_f e^{{-2(L-1)/\xi}}$ L={Lk}")
ax[1, 1].set_xlabel("V"); ax[1, 1].set_ylabel(r"$\delta Q_f$"); ax[1, 1].legend(fontsize=8); ax[1, 1].set_title(r"(e) interacting, $\kappa$=0.2: forward response (MPS)")
ax[1, 2].set_yscale("log"); ax[1, 2].set_xlabel("V"); ax[1, 2].legend(fontsize=7); ax[1, 2].set_title(r"(f) $\kappa$=0.2: reverse response is below the MPS floor $\sim 10^{-4}$")

# resolvable regime 2(L-1)/xi <~ 4: ln(dQf/dQr) vs L, slope = 2/xi_eff
slopes = {}
for col, kap in enumerate([0.03, 0.05]):
    key = f"interacting_k{kap}"; R = rows_of(key)
    if not R: ax[2, col].set_title(f"(g/h) kappa={kap}: pending"); continue
    byV = {}
    for Lk, ok in R.items():
        for r in ok: byV.setdefault(r["V"], []).append((int(Lk), r["lnratio"], r.get("lnratio_err", 0)))
    tp, tm = -1 - 0.5j * kap * 1j, -1 + 0.5j * kap * 1j
    Ls = np.array(sorted({int(k) for k in R}))
    ax[2, col].plot(Ls, 2 * (Ls - 1) * np.log(abs(tm / tp)), "k--", label=r"Gauss $2(L-1)\ln|t_-/t_+|$")
    for V, pts in sorted(byV.items()):
        pts.sort(); x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts]); e = np.array([p[2] for p in pts])
        ax[2, col].errorbar(x, y, e, fmt="o-", capsize=3, label=f"V={V}")
        if len(x) >= 2:
            w = 1 / np.maximum(e, 1e-3); pf, cov = np.polyfit(x, y, 1, w=w, cov="unscaled")
            slopes.setdefault(kap, []).append((V, pf[0], np.sqrt(cov[0, 0])))
    ax[2, col].set_xlabel("L"); ax[2, col].set_ylabel(r"$\ln(\delta Q_f/\delta Q_r)$"); ax[2, col].legend(fontsize=8); ax[2, col].set_title(rf"({'gh'[col]}) $\kappa$={kap}, MPS $\chi$=128")
for kap, sl in slopes.items():
    tp, tm = -1 - 0.5j * kap * 1j, -1 + 0.5j * kap * 1j
    g2 = 2 * np.log(abs(tm / tp))
    ax[2, 2].errorbar([v for v, _, _ in sl], [s_ / g2 for _, s_, _ in sl], [ds / g2 for _, _, ds in sl], fmt="o-", capsize=3, label=rf"$\kappa$={kap}")
ax[2, 2].axhline(1, color="k", ls="--"); ax[2, 2].set_xlabel("V"); ax[2, 2].set_ylabel(r"$\xi_{Gauss}/\xi_{eff}(V)$  (slope ratio)"); ax[2, 2].legend(); ax[2, 2].set_title("(i) interaction renormalisation of the skin length")
res.setdefault("slopes", {}); res["slopes"] = {str(k): [dict(V=v, slope=float(s_), slope_err=float(ds), gauss_slope=float(2 * np.log(abs((-1 + 0.5j * k * 1j) / (-1 - 0.5j * k * 1j))))) for v, s_, ds in sl] for k, sl in slopes.items()}
json.dump(res, open(OUT, "w"), indent=1, default=float)
plt.tight_layout(); plt.savefig("results/beyond_nn.png", dpi=120)
for r in S:
    print(f"kappa={r['kappa']} t2={r['t2']}: ln(dQf/dQr)={r['lnratio']:.3f} (NN law {r['pred_nn']:.3f})  dQf={r['dQf']:.3e} dQr={r['dQr']:.3e} R={r['R']:.6f}  lnT in [{r['lnT_min']:.2f},{r['lnT_max']:.2f}]  w_real={set(r['w_real'])} w_shift={set(r['w_shift'])}")
