"""Stage 2b (interacting): NESS heat currents of the bond-loss chain + V n_j n_{j+1} with mesoscopic leads,
via the superfermion MPS (skinheat.lindblad_mps).  Forward (hot left), reverse (hot right), isothermal.

  python scripts/run_beyond_nn_tebd.py run --L 20 --V 1.0 --cfg f [--K 10 --chi 64 --dt 0.1 --tmax 80]
  python scripts/run_beyond_nn_tebd.py merge        # -> results/beyond_nn.json["interacting"]

Each run writes results/beyond_nn_runs/L{L}_V{V}_{cfg}.json after every time chunk (intermediate results first)."""
import argparse, json, os, time, warnings, numpy as np
warnings.filterwarnings("ignore")
from skinheat import build_chain_with_leads, lead_currents
from skinheat.lindblad_mps import LindbladMPS, layout

KAPPA, PHI, TH, TC, W = 0.2, np.pi / 2, 1.0, 0.5, 4.0
CFG = {"f": (TH, TC), "r": (TC, TH), "iso": (TC, TC)}
RUNS = "results/beyond_nn_runs"


def fname(L, V, cfg, tag="", kappa=KAPPA):
    k = "" if kappa == KAPPA else f"_k{kappa}"
    return f"{RUNS}/L{L}_V{V}_{cfg}{k}{tag}.json"


def run(a):
    TL, TR = CFG[a.cfg]
    m = build_chain_with_leads(a.L, kappa=a.kappa, phi=PHI, gain_ratio=a.gain_ratio, K=a.K, W=W, TL=TL, TR=TR)
    S = LindbladMPS(m, pos=layout(a.L, a.K), V=a.V, int_bonds=[(j, j + 1) for j in range(a.L - 1)])
    out = dict(L=a.L, V=a.V, cfg=a.cfg, K=a.K, chi=a.chi, dt=a.dt, comp=a.comp, kappa=a.kappa, gain_ratio=a.gain_ratio, phi=PHI, TL=TL, TR=TR,
               gauss=lead_currents(m) if a.V == 0 else None, trace=[], converged=False)
    t0 = time.time(); prev = None
    while S.t < a.tmax - 1e-9:
        S.evolve(a.chunk, dt=a.dt, chi=a.chi, max_trunc_err=np.inf, compression_method=a.comp)
        c = S.currents()
        row = dict(t=S.t, EL=c["L"]["E"], ER=c["R"]["E"], NL=c["L"]["N"], NR=c["R"]["N"],
                   tr=float(np.real(S.trace())), chi=int(max(S.psi.chi)), wall=time.time() - t0)
        out["trace"].append(row)
        if prev is not None:
            d = max(abs(row["EL"] - prev["EL"]), abs(row["ER"] - prev["ER"]))
            row["dE"] = d
            if d < a.tol and S.t >= a.tmin:
                out["converged"] = True
        prev = row
        out["final"] = dict(EL=row["EL"], ER=row["ER"], NL=row["NL"], NR=row["NR"], t=row["t"])
        json.dump(out, open(fname(a.L, a.V, a.cfg, a.tag, a.kappa), "w"), indent=1)
        print(f"L={a.L} V={a.V} {a.cfg}: t={S.t:.1f} EL={row['EL']:+.6f} ER={row['ER']:+.6f} chi={row['chi']} dE={row.get('dE', np.nan):.1e} wall={row['wall']:.0f}s", flush=True)
        if out["converged"]:
            break


def merge(a):
    res = json.load(open("results/beyond_nn.json")) if os.path.exists("results/beyond_nn.json") else {}
    sweep, gl = {}, {}
    for L in a.Ls:
        rows = []
        for V in a.Vs:
            r = {}
            for cfg in CFG:
                if os.path.exists(fname(L, V, cfg, a.tag, a.kappa)):
                    r[cfg] = json.load(open(fname(L, V, cfg, a.tag, a.kappa)))
                    tr = r[cfg]["trace"][-5:]
                    r[cfg]["err"] = {k: float(np.std([x[k] for x in tr])) for k in ("EL", "ER")}
            row = dict(V=V, have=sorted(r), converged=all(x["converged"] for x in r.values()) if r else False,
                       t_final={k: x["final"]["t"] for k, x in r.items()}, chi={k: x["trace"][-1]["chi"] for k, x in r.items()})
            if len(r) == 3:
                dQf = -(r["f"]["final"]["ER"] - r["iso"]["final"]["ER"])
                dQr = -(r["r"]["final"]["EL"] - r["iso"]["final"]["EL"])
                ef = np.hypot(r["f"]["err"]["ER"], r["iso"]["err"]["ER"]); er = np.hypot(r["r"]["err"]["EL"], r["iso"]["err"]["EL"])
                row.update(dQf=dQf, dQr=dQr, dQf_err=ef, dQr_err=er, lnratio=float(np.log(abs(dQf / dQr))), R=(abs(dQf) - abs(dQr)) / (abs(dQf) + abs(dQr)),
                           R_lower=(abs(dQf) - ef - abs(dQr) - er) / (abs(dQf) + abs(dQr)), chi_run=r["f"]["chi"], dt=r["f"]["dt"],
                           EL={k: x["final"]["EL"] for k, x in r.items()}, ER={k: x["final"]["ER"] for k, x in r.items()})
                if V == 0 and r["f"]["gauss"]:
                    g = {k: x["gauss"] for k, x in r.items()}
                    gf = -(g["f"]["R"]["E"] - g["iso"]["R"]["E"]); gr = -(g["r"]["L"]["E"] - g["iso"]["L"]["E"])
                    row["gauss"] = dict(dQf=gf, dQr=gr, lnratio=float(np.log(abs(gf / gr))),
                                        ER={k: g[k]["R"]["E"] for k in g}, EL={k: g[k]["L"]["E"] for k in g})
                    gl[str(L)] = row["gauss"]["lnratio"]
            rows.append(row)
        # systematic error from the V=0 Gauss comparison (chi truncation bias), applied to all V at this L
        g0 = next((r for r in rows if r["V"] == 0 and "gauss" in r), None)
        for row in rows:
            if "dQf" in row and g0:
                row["dQf_err"] = max(row["dQf_err"], abs(g0["dQf"] - g0["gauss"]["dQf"]))
                row["dQr_err"] = max(row["dQr_err"], abs(g0["dQr"] - g0["gauss"]["dQr"]))
                row["lnratio_err"] = float(np.hypot(row["dQf_err"] / abs(row["dQf"]), row["dQr_err"] / abs(row["dQr"])))
                row["R_lower"] = (abs(row["dQf"]) - row["dQf_err"] - abs(row["dQr"]) - row["dQr_err"]) / (abs(row["dQf"]) + abs(row["dQr"]))
        sweep[str(L)] = rows
    key = "interacting" if a.kappa == KAPPA and not a.tag else f"interacting_k{a.kappa}{a.tag}"
    res[key] = dict(params=dict(kappa=a.kappa, phi=PHI, Th=TH, Tc=TC, W=W, K=a.K, method="superfermion MPS, W^II order 2, variational compression"),
                    sweep=sweep, gauss_lnratio=gl)
    json.dump(res, open("results/beyond_nn.json", "w"), indent=1, default=float)
    for L, rows in sweep.items():
        for row in rows:
            if "dQf" in row:
                print(f"L={L} V={row['V']}: dQf={row['dQf']:+.5e}({row['dQf_err']:.0e}) dQr={row['dQr']:+.5e}({row['dQr_err']:.0e}) ln ratio={row['lnratio']:.3f} R={row['R']:.5f} (>={row['R_lower']:.3f}) conv={row['converged']} t={row['t_final']} chi={row['chi']}"
                      + (f"  | Gauss: dQf={row['gauss']['dQf']:+.5e} dQr={row['gauss']['dQr']:+.5e} ln={row['gauss']['lnratio']:.3f}" if "gauss" in row else ""))
            else:
                print(f"L={L} V={row['V']}: incomplete {row['have']}")


p = argparse.ArgumentParser(); sub = p.add_subparsers(dest="cmd")
q = sub.add_parser("run"); q.add_argument("--L", type=int, default=20); q.add_argument("--V", type=float, default=0.0); q.add_argument("--cfg", default="f")
q.add_argument("--K", type=int, default=10); q.add_argument("--chi", type=int, default=64); q.add_argument("--dt", type=float, default=0.1)
q.add_argument("--tmax", type=float, default=80.0); q.add_argument("--tmin", type=float, default=20.0); q.add_argument("--chunk", type=float, default=2.0); q.add_argument("--tol", type=float, default=1e-5); q.add_argument("--tag", default=""); q.add_argument("--kappa", type=float, default=KAPPA); q.add_argument("--comp", default="SVD"); q.add_argument("--gain_ratio", type=float, default=0.0)
q = sub.add_parser("merge"); q.add_argument("--Ls", type=int, nargs="+", default=[20, 40]); q.add_argument("--Vs", type=float, nargs="+", default=[0.0, 0.5, 1.0, 2.0]); q.add_argument("--K", type=int, default=10); q.add_argument("--kappa", type=float, default=KAPPA); q.add_argument("--tag", default="")
a = p.parse_args()
{"run": run, "merge": merge}[a.cmd](a)
