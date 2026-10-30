"""Generate Figs. 1-4 of the skinheat manuscript from results/*.json only (no hand-typed numbers)."""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES, FIG = os.path.join(ROOT, "results"), os.path.join(ROOT, "paper", "figures")
os.makedirs(FIG, exist_ok=True)
load = lambda f: json.load(open(os.path.join(RES, f)))
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"], "font.size": 8,
                     "axes.labelsize": 8, "legend.fontsize": 6.5, "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.6,
                     "lines.linewidth": 1.0, "lines.markersize": 4, "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "none",
                     "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"})
W2 = 7.0  # PRL two-column width (inches)


def tpm(kappa, phi=np.pi / 2, J=1.0):
    return -J - 0.5j * kappa * np.exp(1j * phi), -J - 0.5j * kappa * np.exp(-1j * phi)


def inv_xi(kappa, phi=np.pi / 2):
    tp, tm = tpm(kappa, phi); return np.log(np.abs(tm / tp))


def law(kappa, L):
    return 2 * (np.asarray(L) - 1) * inv_xi(kappa)


def label(ax, s, x=-0.2, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontweight="bold", fontsize=9, va="bottom", ha="left")


def save(fig, name):
    for ext, kw in (("pdf", {}), ("svg", {}), ("png", dict(dpi=300))):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), **kw)
    plt.close(fig)


def fig1():
    d = load("deltaT_rectification.json")
    fig, ax = plt.subplots(1, 3, figsize=(W2, 2.2), gridspec_kw=dict(width_ratios=[1.3, 1, 1], wspace=0.5))
    # (a) schematic
    a = ax[0]; n = 7; a.set_xlim(-1.7, n + 0.7); a.set_ylim(-1.75, 1.85); a.axis("off")
    for j in range(n):
        if j < n - 1:
            a.plot([j + 0.22, j + 0.78], [0, 0], "k-", lw=0.8)
            a.annotate("", xy=(j + 0.5, -1.0), xytext=(j + 0.5, -0.12), arrowprops=dict(arrowstyle="-|>", color="C3", lw=0.8, shrinkA=0, shrinkB=0))
    a.scatter(range(n), [0] * n, s=110, fc="white", ec="k", lw=0.8, zorder=3)
    a.text((n - 1) / 2, -1.45, r"$L_j=\sqrt{\kappa}\,(c_j+e^{i\phi}c_{j+1})$", ha="center", va="center", color="C3", fontsize=7)
    a.add_patch(Rectangle((-1.65, -0.55), 0.95, 1.1, fc="#f6c4b8", ec="C3", lw=0.8)); a.text(-1.18, 0, r"$T_L$", ha="center", va="center")
    a.add_patch(Rectangle((n - 0.3, -0.55), 0.95, 1.1, fc="#bcd3f2", ec="C0", lw=0.8)); a.text(n + 0.18, 0, r"$T_R$", ha="center", va="center")
    a.plot([-0.7, -0.22], [0, 0], "k-", lw=0.8); a.plot([n - 1 + 0.22, n - 0.3], [0, 0], "k-", lw=0.8)
    a.text(0.5, 0.36, r"$-J$", ha="center", fontsize=7)
    a.annotate("", xy=(5.7, 1.05), xytext=(0.3, 1.05), arrowprops=dict(arrowstyle="-|>", color="k", lw=1.0))
    a.text((n - 1) / 2, 1.25, r"forward: $\delta Q_f/\delta Q_r=e^{2(L-1)/\xi}$", ha="center", va="bottom", fontsize=7)
    a.text((n - 1) / 2, 0.62, r"skin length $\xi^{-1}=\ln|t_-/t_+|$", ha="center", fontsize=6.5, color="0.35")
    label(a, "(a)", x=0.0)
    # (b) collapse
    b = ax[1]
    for k, mk in zip(["0.05", "0.1", "0.3"], ["o", "s", "^"]):
        rows = d["L"][k]
        b.plot([law(float(k), r["L"]) for r in rows], [r["lnratio"] for r in rows], mk, mfc="none", label=rf"$\kappa={k}J$, $L=5$–$160$")
    rows = [r for r in d["kappa"] if r["kappa"] <= 0.5]
    b.plot([law(r["kappa"], 40) for r in rows], [r["lnratio"] for r in rows], "x", color="C3", label=r"$L=40$, $\kappa=0.003$–$0.5J$")
    b.plot([0, 50], [0, 50], "k-", lw=0.6, zorder=0)
    b.set_xlabel(r"$2(L-1)/\xi$"); b.set_ylabel(r"$\ln(\delta Q_f/\delta Q_r)$"); b.set_xlim(0, 50); b.set_ylim(0, 50)
    b.legend(loc="upper left", handlelength=1); label(b, "(b)")
    # (c) rectification vs throughput
    c = ax[2]
    rows = [r for r in d["kappa"] if r["kappa"] <= 1.0]; kap = np.array([r["kappa"] for r in rows])
    kk = np.logspace(-2.5, 0, 200); c.semilogx(kk, np.tanh(39 * inv_xi(kk)), "k--", lw=0.6)
    c.semilogx(kap, [r["R"] for r in rows], "o", color="k", ms=3)
    c.set_xlabel(r"$\kappa/J$"); c.set_ylabel(r"$\mathcal{R}$ (black)"); c.set_ylim(0, 1.08)
    c2 = c.twinx(); c2.semilogx(kap, [r["dQf"] * 1e2 for r in rows], "s-", color="C3", ms=3)
    c2.set_ylabel(r"$\delta Q_f\ (10^{-2}J^2/\hbar)$", color="C3"); c2.tick_params(axis="y", colors="C3"); c2.set_ylim(0, 4.6)
    c.axvline(0.07, color="0.6", lw=0.5, ls=":"); c.text(0.078, 0.08, r"$\mathcal{R}>0.99$", fontsize=6.5, rotation=90, va="bottom")
    label(c, "(c)")
    fig.subplots_adjust(left=0.03, right=0.93, top=0.9, bottom=0.2)
    save(fig, "fig1")
    r07 = min(d["kappa"], key=lambda r: abs(r["kappa"] - 0.07)); r0 = d["kappa"][0]
    print(f"fig1: kappa={r07['kappa']:.3f}: R={r07['R']:.4f}, dQf={r07['dQf']:.4f}; kappa={r0['kappa']:.4f}: dQf={r0['dQf']:.4f}")


def fig2():
    b = load("beyond_nn.json")["nnn"]; sw = b["sweep"]; L = b["params"]["L"]
    fig, ax = plt.subplots(1, 3, figsize=(W2, 2.2), gridspec_kw=dict(width_ratios=[1.3, 0.9, 1], wspace=0.5))
    a = ax[0]
    for r, col in zip(sorted([r for r in sw if r["kappa"] == 0.2], key=lambda r: r["t2"]), ["k", "C0", "C2", "C3"]):
        a.plot(r["eps"], r["lnT"], "-", color=col, label=rf"$t_2={r['t2']}J$")
        if r["t2"]:
            a.plot(r["eps"], r["nonbloch"], ":", color=col, lw=1.3)
    a.axhline(0, color="0.7", lw=0.5); a.set_ylim(-5, 40); a.set_xlim(-4, 4)
    a.set_xlabel(r"$\varepsilon/J$"); a.set_ylabel(r"$\ln[T_{RL}(\varepsilon)/T_{LR}(\varepsilon)]$"); a.legend(loc="upper right", bbox_to_anchor=(1.0, 0.86))
    a.text(0.97, 0.96, rf"$\kappa=0.2J$, $L={L}$" + "\n" + r"dotted: $2(L-1)\ln|\beta_2\beta_3|(\varepsilon)$", transform=a.transAxes, va="top", ha="right", fontsize=6.5)
    label(a, "(a)", x=-0.17)
    bx = ax[1]; k = np.linspace(-np.pi, np.pi, 600); kappa = 0.2
    for t2, col in [(0.0, "k"), (0.3, "C2")]:
        E = -2 * np.cos(k) + 2 * t2 * np.cos(2 * k) - 1j * kappa * (1 - np.sin(k))
        bx.plot(E.real, E.imag, "-", color=col, label=rf"$t_2={t2}J$")
    bx.axhline(0, color="C3", lw=0.8)
    bx.text(2.65, 0.012, r"$\varepsilon\in\mathbb{R}$: $w(\varepsilon)=0$", color="C3", fontsize=6.5, ha="right", va="bottom")
    bx.set_xlabel(r"Re$\,E/J$"); bx.set_ylabel(r"Im$\,E/J$"); bx.legend(loc="upper left"); bx.set_ylim(-0.45, 0.12)
    bx.text(0.97, 0.05, r"PBC spectrum, $\kappa=0.2J$", transform=bx.transAxes, ha="right", fontsize=6.5)
    label(bx, "(b)", x=-0.3)
    c = ax[2]
    for kap, mk in zip([0.1, 0.2, 0.3], ["o", "s", "^"]):
        rows = sorted([r for r in sw if r["kappa"] == kap], key=lambda r: r["t2"])
        l, = c.plot([r["t2"] for r in rows], [r["lnratio"] for r in rows], mk + "-", label=rf"$\kappa={kap}J$")
        c.axhline(rows[0]["pred_nn"], color=l.get_color(), ls="--", lw=0.6)
    c.set_xlabel(r"$t_2/J$"); c.set_ylabel(r"$\ln(\delta Q_f/\delta Q_r)$"); c.legend(loc="upper right"); c.set_ylim(0, 27)
    Rmin = min(r["R"] for r in sw)
    c.text(0.03, 0.03, r"dashed: NN law $2(L-1)\ln|t_-/t_+|$" + "\n" + rf"$\mathcal{{R}}>{Rmin:.2f}$ for all points", transform=c.transAxes, fontsize=6.5)
    label(c, "(c)")
    fig.subplots_adjust(left=0.07, right=0.98, top=0.9, bottom=0.2)
    save(fig, "fig2")
    print("fig2: Rmin =", Rmin, "; kappa=0.2 thermal ln ratios:", [(r["t2"], round(r["lnratio"], 2)) for r in sorted([r for r in sw if r["kappa"] == 0.2], key=lambda r: r["t2"])])


def fig3():
    B = load("beyond_nn.json"); cf = B["chi256_final"]; ck = B["chi256_kappa0.1"]; gb = B["gain_bath_V"]
    e05, e1, eg = cf["err_estimate"], 0.04, gb["gain_bath"]["err_V2"]
    fig, ax = plt.subplots(1, 4, figsize=(W2, 2.15), gridspec_kw=dict(wspace=0.6))
    a = ax[0]; Ls = np.array([10, 22])
    for kap, col in [(0.05, "C0"), (0.1, "C3")]:
        a.plot(Ls, law(kap, Ls), "-", color=col, lw=0.8)
    a.errorbar([16], [cf["ln_V0_chi256"]["16"]], yerr=e05, fmt="o", mfc="none", color="C0", label=r"MPS $V=0$")
    a.errorbar([16, 20], [cf["ln_V2"]["16"], cf["ln_V2"]["20"]], yerr=e05, fmt="o", color="C0", label=r"MPS $V=2J$")
    a.errorbar([12, 16], [ck["12"]["0.0"], ck["16"]["0.0"]], yerr=e1, fmt="s", mfc="none", color="C3")
    a.errorbar([12, 16], [ck["12"]["2.0"], ck["16"]["2.0"]], yerr=e1, fmt="s", color="C3")
    a.errorbar([12.4], [gb["gain_bath"]["0.0"]], yerr=eg, fmt="D", mfc="none", color="C2")
    a.errorbar([12.4], [gb["gain_bath"]["2.0"]], yerr=eg, fmt="D", color="C2", label=r"$T_b=\infty$ bath (open: $V=0$)")
    a.set_xlabel(r"$L$"); a.set_ylabel(r"$\ln(\delta Q_f/\delta Q_r)$"); a.set_xlim(10, 22); a.set_ylim(0.5, 3.7)
    a.text(21.7, law(0.05, 21.7) + 0.08, r"$\kappa=0.05J$", color="C0", ha="right", va="bottom", fontsize=6)
    a.text(17.2, 3.02, r"$\kappa=0.1J$", color="C3", ha="left", va="top", fontsize=6)
    a.legend(loc="upper left", fontsize=5.5, handlelength=1.2, labelspacing=0.3); label(a, "(a)")
    b = ax[1]
    b.errorbar([16, 20], [cf["offsets"]["16"], cf["offsets"]["20"]], yerr=e05, fmt="o", color="C0", label=r"$\kappa=0.05J$")
    b.errorbar([12, 16], [ck["12"]["cv"] - ck["12"]["gauss"], ck["16"]["cv"] - ck["16"]["gauss"]], yerr=e1, fmt="s", color="C3", label=r"$\kappa=0.1J$")
    b.errorbar([12.4], [gb["gain_bath"]["offset_V2"]], yerr=eg, fmt="D", color="C2", label=r"$\kappa=0.1J$, $T_b=\infty$")
    b.axhline(0, color="0.6", lw=0.5); b.set_ylim(-0.8, 0.1); b.set_xlim(10, 22)
    b.set_xlabel(r"$L$"); b.set_ylabel(r"$\Delta_V=\ln r(V{=}2J)-\ln r(0)$"); b.legend(loc="lower left", fontsize=5.5); label(b, "(b)")
    c = ax[2]
    for kap, mk, col in [("0.03", "^", "C1"), ("0.05", "o", "C0")]:
        rows = B["slopes"][kap]; g = rows[0]["gauss_slope"]
        c.errorbar([r["V"] for r in rows], [r["slope"] / g for r in rows], yerr=[r["slope_err"] / g for r in rows], fmt=mk, mfc="none", color=col, label=rf"$\kappa={kap}J$, $\chi=128$")
    c.errorbar([2.07], [cf["slope_V2"] / cf["slope_gauss"]], yerr=0.1, fmt="o", color="C0", label=r"$\kappa=0.05J$, $\chi=256$")
    s01 = (ck["16"]["cv"] - ck["12"]["cv"]) / 4; g01 = (ck["16"]["gauss"] - ck["12"]["gauss"]) / 4
    c.errorbar([2.14], [s01 / g01], yerr=0.1, fmt="s", color="C3", label=r"$\kappa=0.1J$, $\chi=256$")
    c.axhline(1, color="k", lw=0.6); c.set_ylim(0.4, 1.7); c.set_xlim(-0.2, 2.4)
    c.set_xlabel(r"$V/J$"); c.set_ylabel(r"$\xi/\xi_{\rm eff}(V)$"); c.legend(loc="upper left", fontsize=5.3, labelspacing=0.25); label(c, "(c)")
    d = ax[3]
    for key, L, mk, col, ls in [("interacting_k0.05", "20", "o", "C0", "-"), ("interacting_k0.05", "12", "o", "C0", "--"),
                                ("interacting_k0.03", "20", "^", "C1", "-"), ("interacting_k0.1", "12", "s", "C3", "--")]:
        rows = sorted(B[key]["sweep"][L], key=lambda r: r["V"]); q0 = rows[0]["dQf"]
        d.plot([r["V"] for r in rows], [r["dQf"] / q0 for r in rows], mk + ls, color=col, mfc=("none" if ls == "--" else col), label=rf"$\kappa={B[key]['params']['kappa']}J$, $L={L}$")
    d.set_xlabel(r"$V/J$"); d.set_ylabel(r"$\delta Q_f(V)/\delta Q_f(0)$"); d.set_ylim(0.5, 1.1); d.legend(fontsize=5.5, loc="lower left"); label(d, "(d)")
    fig.subplots_adjust(left=0.07, right=0.99, top=0.9, bottom=0.21)
    save(fig, "fig3")
    print(f"fig3: slope ratio k0.05 chi256 = {cf['slope_V2']/cf['slope_gauss']:.3f}; k0.1 chi256 = {s01/g01:.3f}; offsets k0.1 CV = {ck['12']['cv']-ck['12']['gauss']:.3f}, {ck['16']['cv']-ck['16']['gauss']:.3f}")


def fig4():
    T = load("thermal_bulk.json"); D = load("disorder.json")
    fig, axs = plt.subplots(2, 2, figsize=(W2, 4.0), gridspec_kw=dict(wspace=0.6, hspace=0.45)); ax = axs.ravel()
    a = ax[0]
    rows = [r for r in T["kappa_Tb"] if r["bulk_bath"]["T"] == 0.5 and r["bulk_bath"]["mu"] == 0.0]
    kap = np.array([r["kappa"] for r in rows]); cost = np.array([r["f"]["sigma"] / r["dQf"] for r in rows])
    a.plot(kap, cost, "o-", color="k"); a.plot(kap, cost[0] / kap[0] * kap, "k:", lw=0.6)
    a.set_xlabel(r"$\kappa/J$"); a.set_ylabel(r"$\sigma_f/\delta Q_f\ \ (1/J)$"); a.set_xlim(0, 0.55); a.set_ylim(0, 20)
    a2 = a.twinx(); a2.plot(kap, [r["dQf"] * 1e2 for r in rows], "s-", color="C3")
    a2.set_ylabel(r"$\delta Q_f\ (10^{-2}J^2/\hbar)$", color="C3"); a2.tick_params(axis="y", colors="C3"); a2.set_ylim(0, 2)
    a.text(0.97, 0.42, r"$T_b=T_c$, $\mu_b=0$, $L=40$" + "\n" + r"dotted: $\propto\kappa\propto1/\xi$", transform=a.transAxes, va="top", ha="right", fontsize=6.5)
    label(a, "(a)")
    b = ax[1]; rows = T["mu_sweep"]; mu = np.array([-r["bulk_bath"]["mu"] for r in rows])
    b.plot(mu, [r["f"]["sigma"] for r in rows], "o-", color="k"); b.set_yscale("log"); b.set_xscale("symlog", linthresh=1.0)
    b.set_xlabel(r"$-\mu_b/J$"); b.set_ylabel(r"$\sigma_f\ (J/\hbar)$")
    b2 = b.twinx(); b2.plot(mu, [r["f"]["n_mean"] for r in rows], "s-", color="C0"); b2.set_ylabel(r"filling $\bar n$", color="C0"); b2.tick_params(axis="y", colors="C0"); b2.set_ylim(0, 0.62)
    label(b, "(b)")
    c = ax[2]
    for kd, col in [(0.05, "C0"), (0.1, "C3"), (0.2, "k")]:
        rows = [r for r in D["onsite"] if r["kappa"] == kd]; W = [r["Wd"] for r in rows]
        c.semilogy(W, [r["Qf_median"] for r in rows], "o-", color=col, ms=3, label=rf"$\kappa={kd}J$")
        c.fill_between(W, [r["Qf_p10"] for r in rows], [r["Qf_p90"] for r in rows], color=col, alpha=0.2, lw=0)
    c.set_xlabel(r"$W_d/J$"); c.set_ylabel(r"$\delta Q_f\ (J^2/\hbar)$, median"); c.legend(loc="lower left", fontsize=5.5); c.set_ylim(1e-7, 6e-2)
    c2 = c.twinx()
    for kd, col in [(0.05, "C0"), (0.1, "C3"), (0.2, "k")]:
        rows = [r for r in D["onsite"] if r["kappa"] == kd]; c2.plot([r["Wd"] for r in rows], [r["R_mean"] for r in rows], "--", color=col, lw=0.8)
    c2.set_ylim(0.95, 1.004); c2.set_ylabel(r"$\mathcal{R}$ (dashed)")
    label(c, "(c)")
    d = ax[3]
    for kind, mk, col, lab in [("kappa", "o", "C0", r"$\kappa_j=\kappa(1\pm0.8u)$"), ("phi", "s", "C2", r"$\phi_j=\pi/2\pm1$"), ("sign", "^", "C3", r"$\phi_j=\pm\pi/2$ (random chirality)")]:
        rows = sorted([r for r in D["bond"] if r["kind"] == kind], key=lambda r: r["L"]); Ls = np.array([r["L"] for r in rows]); s = np.array([r["ln_std"] for r in rows])
        d.plot(Ls, s, mk, color=col, label=lab); d.plot(Ls, s[1] * np.sqrt(Ls / 40), ":", color=col, lw=0.6)
    d.set_xlabel(r"$L$"); d.set_ylabel(r"std of $\ln(\delta Q_f/\delta Q_r)$"); d.legend(loc="upper left", fontsize=5.5); d.set_xlim(10, 90); d.set_ylim(0, 2.3)
    d.text(0.97, 0.97, r"dotted: $\propto\sqrt{L}$, $\kappa=0.1J$", transform=d.transAxes, ha="right", va="top", fontsize=6.5)
    label(d, "(d)")
    fig.subplots_adjust(left=0.08, right=0.92, top=0.95, bottom=0.1)
    save(fig, "fig4")
    print("fig4: cost sigma_f/dQf:", [(k, round(c, 2)) for k, c in zip(kap, cost)])
    print("fig4: bond means vs clean:", [(r["L"], r["kind"], round(r["ln_mean"], 2), round(r["clean"], 2)) for r in D["bond"]])


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4()
