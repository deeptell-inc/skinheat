"""Stage 2c: Delta-T responses by the exact perturbative identity (no subtraction of two NESS):
    rho_f - rho_iso = int_0^inf exp(L_f s) (L_f - L_iso) rho_iso ds   (and likewise for r),
so dQ_f = -int E_R[chi_f(s)] ds, dQ_r = -int E_L[chi_r(s)] ds with chi the decaying perturbation (skinheat.lindblad_mps.LindbladResponse).
One base NESS (isothermal) + two perturbation evolutions per (kappa, L, V).

  python scripts/run_response.py run --kappa 0.05 --L 20 --V 1.0 [--K 10 --chi 128 --dt 0.1 --tbase 100 --tresp 100]
  python scripts/run_response.py merge --kappa 0.05 --Ls 12 16 20     -> results/beyond_nn.json["response_k<kappa>"]
Writes results/beyond_nn_runs/resp_L{L}_V{V}_k{kappa}.json incrementally (base trace, then f, then r)."""
import argparse, json, os, pickle, time, warnings, numpy as np
warnings.filterwarnings("ignore")
from skinheat import build_chain_with_leads, lead_currents
from skinheat.lindblad_mps import LindbladMPS, LindbladResponse, layout

PHI, TH, TC, W = np.pi / 2, 1.0, 0.5, 4.0
RUNS = "results/beyond_nn_runs"


def fname(L, V, kappa, tag=""):
    return f"{RUNS}/resp_L{L}_V{V}_k{kappa}{tag}.json"


def run(a):
    mk = lambda TL, TR: build_chain_with_leads(a.L, kappa=a.kappa, phi=PHI, K=a.K, W=W, TL=TL, TR=TR)
    m_iso, m_f, m_r = mk(TC, TC), mk(TH, TC), mk(TC, TH)
    g = {c: lead_currents(m) for c, m in [("iso", m_iso), ("f", m_f), ("r", m_r)]} if a.V == 0 else None
    out = dict(L=a.L, V=a.V, kappa=a.kappa, K=a.K, chi=a.chi, dt=a.dt, method="response",
               gauss=dict(dQf=-(g["f"]["R"]["E"] - g["iso"]["R"]["E"]), dQr=-(g["r"]["L"]["E"] - g["iso"]["L"]["E"])) if g else None,
               base_trace=[], resp={})
    fn = fname(a.L, a.V, a.kappa, a.tag); t0 = time.time()
    base = LindbladMPS(m_iso, pos=layout(a.L, a.K), V=a.V, int_bonds=[(j, j + 1) for j in range(a.L - 1)])
    prev = None
    while base.t < a.tbase - 1e-9:
        base.evolve(2.0, dt=a.dt, chi=a.chi, compression_method="SVD", max_trunc_err=np.inf)
        c = base.currents(); row = dict(t=base.t, EL=c["L"]["E"], ER=c["R"]["E"], chi=int(max(base.psi.chi)), wall=time.time() - t0)
        if prev: row["dE"] = max(abs(row["EL"] - prev["EL"]), abs(row["ER"] - prev["ER"]))
        out["base_trace"].append(row); prev = row
        json.dump(out, open(fn, "w"), indent=1)
        if base.t >= a.tmin and row.get("dE", 1) < a.tol:
            break
    out["base"] = dict(EL=prev["EL"], ER=prev["ER"], t=base.t)
    pickle.dump(base.psi, open(fn.replace(".json", "_base.pkl"), "wb"))
    print(f"base done t={base.t} EL={prev['EL']:+.6f} ER={prev['ER']:+.6f} ({time.time()-t0:.0f}s)", flush=True)
    for cfg, m_t, lab in [("f", m_f, "R"), ("r", m_r, "L")]:
        t1 = time.time(); R = LindbladResponse(base, m_t, chi=a.chi)
        while R.T.t < a.tresp - 1e-9:
            R.run(min(R.T.t + 10.0, a.tresp), dt=a.dt)
            out["resp"][cfg] = dict(ts=R.ts, E=[e[lab] for e in R.E], E_other=[e["L" if lab == "R" else "R"] for e in R.E],
                                    integral=-R.integral(lab), integral_notail=-R.integral(lab, tail=False), wall=time.time() - t1)
            json.dump(out, open(fn, "w"), indent=1)
            print(f"{cfg}: t={R.T.t:.0f} E_{lab}={R.E[-1][lab]:+.3e} dQ={out['resp'][cfg]['integral']:+.6e} (notail {out['resp'][cfg]['integral_notail']:+.6e}) ({time.time()-t1:.0f}s)", flush=True)
            if abs(R.E[-1][lab]) < a.etol * abs(out["resp"][cfg]["integral"]) and R.T.t >= 20:
                break
    dQf, dQr = out["resp"]["f"]["integral"], out["resp"]["r"]["integral"]
    out["final"] = dict(dQf=dQf, dQr=dQr, lnratio=float(np.log(abs(dQf / dQr))), R=(abs(dQf) - abs(dQr)) / (abs(dQf) + abs(dQr)))
    json.dump(out, open(fn, "w"), indent=1)
    print("final", out["final"], "gauss", out["gauss"], flush=True)


def merge(a):
    res = json.load(open("results/beyond_nn.json"))
    sweep = {}
    for L in a.Ls:
        rows = []
        for V in a.Vs:
            fn = fname(L, V, a.kappa, a.tag)
            if not os.path.exists(fn):
                continue
            d = json.load(open(fn))
            if "final" not in d:
                rows.append(dict(V=V, incomplete=True)); continue
            row = dict(V=V, **d["final"], base=d["base"], chi=d["chi"], dt=d["dt"], gauss=d["gauss"],
                       tail_frac={c: abs(d["resp"][c]["integral"] - d["resp"][c]["integral_notail"]) / abs(d["resp"][c]["integral"]) for c in "fr"})
            rows.append(row)
        # systematic error from the V=0 Gauss comparison (same L, kappa)
        g0 = next((r for r in rows if r.get("V") == 0 and r.get("gauss")), None)
        for r in rows:
            if "dQf" in r:
                ef = abs(g0["dQf"] - g0["gauss"]["dQf"]) / abs(g0["gauss"]["dQf"]) if g0 else 0.05
                er = abs(g0["dQr"] - g0["gauss"]["dQr"]) / abs(g0["gauss"]["dQr"]) if g0 else 0.05
                r["dQf_err"], r["dQr_err"] = ef * abs(r["dQf"]), er * abs(r["dQr"])
                r["lnratio_err"] = float(np.hypot(ef, er))
        sweep[str(L)] = rows
    res[f"response_k{a.kappa}{a.tag}"] = dict(params=dict(kappa=a.kappa, phi=PHI, Th=TH, Tc=TC, W=W, K=a.K, method="exact perturbative response, SVD-compressed W^II"), sweep=sweep)
    json.dump(res, open("results/beyond_nn.json", "w"), indent=1, default=float)
    for L, rows in sweep.items():
        for r in rows:
            if "dQf" in r:
                print(f"L={L} V={r['V']}: dQf={r['dQf']:+.5e} dQr={r['dQr']:+.5e} ln={r['lnratio']:.4f} ± {r['lnratio_err']:.3f} R={r['R']:.4f} tail={r['tail_frac']}" + (f" | Gauss dQf={r['gauss']['dQf']:+.5e} dQr={r['gauss']['dQr']:+.5e} ln={np.log(r['gauss']['dQf']/r['gauss']['dQr']):.4f}" if r["gauss"] else ""))
            else:
                print(f"L={L} V={r['V']}: incomplete")


p = argparse.ArgumentParser(); sub = p.add_subparsers(dest="cmd")
q = sub.add_parser("run")
for n, t, d in [("L", int, 20), ("V", float, 0.0), ("kappa", float, 0.05), ("K", int, 10), ("chi", int, 128), ("dt", float, 0.1), ("tbase", float, 100.0), ("tmin", float, 40.0), ("tol", float, 1e-5), ("tresp", float, 120.0), ("etol", float, 1e-4)]:
    q.add_argument(f"--{n}", type=t, default=d)
q.add_argument("--tag", default="")
q = sub.add_parser("merge"); q.add_argument("--Ls", type=int, nargs="+", default=[12, 16, 20]); q.add_argument("--Vs", type=float, nargs="+", default=[0.0, 0.5, 1.0, 2.0]); q.add_argument("--kappa", type=float, default=0.05); q.add_argument("--K", type=int, default=10); q.add_argument("--tag", default="")
a = p.parse_args()
{"run": run, "merge": merge}[a.cmd](a)
