"""chi-convergence + Gauss control variate for the kappa=0.05 subtraction runs (SVD, dt=0.1).
For each (L, V, chi): dQf, dQr, ln ratio.  V=0 is exact (Gauss) -> bias(chi) = MPS - Gauss.
Control variate: dQ_cv(V) = dQ_mps(V) - dQ_mps(0) + dQ_gauss(0); validated by the chi-stability of dQ_mps(V) - dQ_mps(0).
Writes results/beyond_nn.json["chi_convergence"] and results/chi_convergence.png."""
import json, os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from skinheat import build_chain_with_leads, lead_currents

RUNS, KAP = "results/beyond_nn_runs", 0.05
CFG = {"f": (1.0, 0.5), "r": (0.5, 1.0), "iso": (0.5, 0.5)}


def load(L, V, chi, dt=0.1):
    tag = ("" if chi == 128 else f"_chi{chi}") + ("" if dt == 0.1 else f"_dt{dt}")
    r = {}
    for c in CFG:
        fn = f"{RUNS}/L{L}_V{V}_{c}_k{KAP}{tag}.json"
        if not os.path.exists(fn): return None
        d = json.load(open(fn))
        if "final" not in d or d["final"]["t"] < 40: return None   # not yet relaxed (run still in progress)
        r[c] = d["final"]
    dQf = -(r["f"]["ER"] - r["iso"]["ER"]); dQr = -(r["r"]["EL"] - r["iso"]["EL"])
    return dict(dQf=dQf, dQr=dQr, ln=float(np.log(dQf / dQr)))


def gauss(L):
    g = {c: lead_currents(build_chain_with_leads(L, kappa=KAP, phi=np.pi / 2, K=10, TL=TL, TR=TR)) for c, (TL, TR) in CFG.items()}
    dQf = -(g["f"]["R"]["E"] - g["iso"]["R"]["E"]); dQr = -(g["r"]["L"]["E"] - g["iso"]["L"]["E"])
    return dict(dQf=dQf, dQr=dQr, ln=float(np.log(dQf / dQr)))


out = {}
for L in [16, 20]:
    G = gauss(L); out[str(L)] = dict(gauss=G, chi={})
    print(f"L={L}  Gauss: dQf={G['dQf']:.5e} dQr={G['dQr']:.5e} ln={G['ln']:.4f}")
    for chi in [64, 96, 128, 192, 256]:
        m0, m2 = load(L, 0.0, chi), load(L, 2.0, chi)
        if not m0: continue
        row = dict(V0=m0, bias=dict(dQf=m0["dQf"] - G["dQf"], dQr=m0["dQr"] - G["dQr"], ln=m0["ln"] - G["ln"]))
        line = f"  chi={chi:3d}: V=0 bias dQf {row['bias']['dQf']/G['dQf']:+.2%} dQr {row['bias']['dQr']/G['dQr']:+.2%} ln {row['bias']['ln']:+.4f}"
        if m2:
            row["V2"] = m2
            row["cv"] = dict(dQf=m2["dQf"] - m0["dQf"] + G["dQf"], dQr=m2["dQr"] - m0["dQr"] + G["dQr"])
            row["cv"]["ln"] = float(np.log(row["cv"]["dQf"] / row["cv"]["dQr"]))
            row["diff"] = dict(dQf=m2["dQf"] - m0["dQf"], dQr=m2["dQr"] - m0["dQr"], ln=m2["ln"] - m0["ln"])
            line += f" | V=2 raw ln={m2['ln']:.4f} dQr={m2['dQr']:.4e} ; V2-V0: dQf {row['diff']['dQf']:+.4e} dQr {row['diff']['dQr']:+.4e} dln {row['diff']['ln']:+.4f} ; CV ln={row['cv']['ln']:.4f}"
        out[str(L)]["chi"][str(chi)] = row; print(line)
# ---- dt series at chi=128 (+ Richardson extrapolation, error ~ dt^2) ----
dts = {}
for L in [16, 20]:
    G = out[str(L)]["gauss"]; dts[str(L)] = {}
    for V in [0.0, 2.0]:
        ser = {dt: load(L, V, 128, dt) for dt in [0.1, 0.05, 0.025]}
        ser = {dt: v for dt, v in ser.items() if v}
        row = {str(dt): v for dt, v in ser.items()}
        fine = sorted(ser)[:2]                      # two smallest dt available
        if len(fine) == 2:
            d1, d2 = fine                            # d1 < d2
            row["richardson"] = {k: (d2 * ser[d1][k] - d1 * ser[d2][k]) / (d2 - d1) for k in ("dQf", "dQr")}   # linear-in-dt extrapolation
            row["richardson"]["ln"] = float(np.log(row["richardson"]["dQf"] / row["richardson"]["dQr"]))
            row["richardson"]["from"] = fine
        c256 = load(L, V, 256, 0.05)
        if c256:
            row["chi256_dt0.05"] = c256
        dts[str(L)][str(V)] = row
        line = f"L={L} V={V}: " + " ".join(f"dt={dt}: ln={v['ln']:.4f} dQr={v['dQr']:.4e};" for dt, v in sorted(ser.items(), reverse=True))
        if "richardson" in row: line += f"  linear dt->0 from {row['richardson']['from']}: ln={row['richardson']['ln']:.4f} dQr={row['richardson']['dQr']:.4e}"
        if "chi256_dt0.05" in row: line += f"  chi=256,dt=0.05: ln={row['chi256_dt0.05']['ln']:.4f} dQr={row['chi256_dt0.05']['dQr']:.4e}"
        if V == 0: line += f"  | Gauss ln={G['ln']:.4f} dQr={G['dQr']:.4e}"
        print(line)
res = json.load(open("results/beyond_nn.json")); res["chi_convergence"] = out; res["dt_convergence"] = dts
json.dump(res, open("results/beyond_nn.json", "w"), indent=1, default=float)

fig, ax = plt.subplots(1, 4, figsize=(20, 4.5))
for L, mk in zip([16, 20], "os"):
    d = out[str(L)]; chis = sorted(int(c) for c in d["chi"]); G = d["gauss"]
    ax[0].plot(chis, [d["chi"][str(c)]["bias"]["dQr"] / G["dQr"] * 100 for c in chis], mk + "-", label=f"L={L}: $\\delta Q_r$")
    ax[0].plot(chis, [d["chi"][str(c)]["bias"]["dQf"] / G["dQf"] * 100 for c in chis], mk + "--", mfc="none", label=f"L={L}: $\\delta Q_f$")
    c2 = [c for c in chis if "V2" in d["chi"][str(c)]]
    ax[1].plot(c2, [d["chi"][str(c)]["diff"]["ln"] for c in c2], mk + "-", label=f"L={L}: raw $\\Delta\\ln$ (V=2 − V=0)")
    ax[2].plot(c2, [d["chi"][str(c)]["V2"]["ln"] for c in c2], mk + "--", mfc="none", label=f"L={L}: V=2 raw")
    ax[2].plot(c2, [d["chi"][str(c)]["cv"]["ln"] for c in c2], mk + "-", label=f"L={L}: V=2 control-variate")
    ax[2].axhline(G["ln"], color="gray", ls=":", lw=0.8)
ax[0].axhline(0, color="k", lw=0.8); ax[0].set_xlabel(r"$\chi$"); ax[0].set_ylabel("V=0 bias vs Gauss [%]"); ax[0].legend(fontsize=8); ax[0].set_title(r"(a) truncation bias at V=0 ($\kappa$=0.05, SVD, dt=0.1)")
ax[1].set_xlabel(r"$\chi$"); ax[1].set_ylabel(r"$\ln r(V{=}2)-\ln r(V{=}0)$"); ax[1].legend(fontsize=8); ax[1].set_title("(b) interaction shift: chi-stability")
ax[2].set_xlabel(r"$\chi$"); ax[2].set_ylabel(r"$\ln(\delta Q_f/\delta Q_r)$ at V=2"); ax[2].legend(fontsize=8); ax[2].set_title("(c) raw vs bias-corrected (dotted: Gauss V=0)")
for L, mk in zip([16, 20], "os"):
    G = out[str(L)]["gauss"]
    for V, ls in [(0.0, "-"), (2.0, "--")]:
        row = dts[str(L)].get(str(V), {}); pts = sorted((float(k), v["ln"]) for k, v in row.items() if k not in ("richardson", "chi256_dt0.05"))
        if pts:
            ax[3].plot([p[0] for p in pts], [p[1] for p in pts], mk + ls, label=f"L={L}, V={V}")
            if "richardson" in row: ax[3].plot([0], [row["richardson"]["ln"]], mk, mfc="none", ms=9)
            if "chi256_dt0.05" in row: ax[3].plot([0.05], [row["chi256_dt0.05"]["ln"]], "P", color="k", ms=10, label=f"L={L}, V={V}: $\\chi$=256" )
    ax[3].plot([0], [G["ln"]], "k*", ms=10)
ax[3].set_xlabel("dt"); ax[3].set_ylabel(r"$\ln(\delta Q_f/\delta Q_r)$"); ax[3].legend(fontsize=8); ax[3].set_title(r"(d) dt series, $\chi$=128; open: linear $dt \to 0$; star: Gauss V=0")
plt.tight_layout(); plt.savefig("results/chi_convergence.png", dpi=130)
