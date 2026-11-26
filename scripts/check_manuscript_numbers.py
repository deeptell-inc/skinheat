"""Recompute every number quoted in paper/prb/main.tex from results/*.json (or directly from the package)
and compare with the quoted value.  Prints a Markdown table; exits 1 if any row fails.

    python scripts/check_manuscript_numbers.py > paper/prb/number-check.md
"""
import json, os, sys, glob
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from skinheat import chain_heff, delta_heat, build_chain_with_leads, lead_currents, ness_covariance, entropy_production

R = lambda f: json.load(open(os.path.join(ROOT, "results", f)))
DT, BN, TB, DS = R("deltaT_rectification.json"), R("beyond_nn.json"), R("thermal_bulk.json"), R("disorder.json")
rows, fails = [], 0


def chk(sec, claim, quoted, value, tol=None, src="", op="eq"):
    """op: 'eq' |value-quoted|<=tol ; 'le' value<=quoted ; 'ge' value>=quoted ; 'info' no test."""
    global fails
    q = np.atleast_1d(np.array(quoted, float)); v = np.atleast_1d(np.array(value, float))
    if op == "eq":
        ok = v.shape == q.shape and bool(np.all(np.abs(v - q) <= tol))
    elif op == "le":
        ok = bool(np.all(v <= q))
    elif op == "ge":
        ok = bool(np.all(v >= q))
    else:
        ok = None
    fails += ok is False
    fmt = lambda a: ", ".join(f"{x:.5g}" for x in a)
    rows.append(f"| {sec} | {claim} | {fmt(q)} | {fmt(v)} | {src} | {'OK' if ok else ('—' if ok is None else '**FAIL**')} |")


def inv_xi(k):
    return np.log(abs((-1 - 0.5j * k * np.exp(-1j * np.pi / 2)) / (-1 - 0.5j * k * np.exp(1j * np.pi / 2))))


# ---------------- Sec. III, Fig. 1 ----------------
chk("III B", "threshold (L-1)/xi for R>0.99", 2.65, np.arctanh(0.99), 0.005, "atanh(0.99)")
r07 = min(DT["kappa"], key=lambda r: abs(r["kappa"] - 0.07))
q0 = delta_heat(chain_heff(40, kappa=0.0), 1.0, 0.5)[0]
chk("III B", "R at kappa=0.07J, L=40", 0.992, r07["R"], 0.0005, "deltaT_rectification.json kappa[]")
chk("III B", "dQf(kappa=0.07)/dQf(kappa=0) in %", 45, 100 * r07["dQf"] / q0, 1.0, "json + chain_heff(kappa=0)")
dev = max(abs(r["lnratio"] - 2 * (r["L"] - 1) * inv_xi(float(k))) for k, v in DT["L"].items() for r in v if r["lnratio"] < 30)
chk("III B", "Fig. 1(b): max deviation from the law (ln ratio < 30)", 1e-8, dev, op="le", src="deltaT_rectification.json L{}")
chk("III B", "L range in Fig. 1(b)", [5, 160], [min(r["L"] for r in DT["L"]["0.1"]), max(r["L"] for r in DT["L"]["0.1"])], 0, "json")
kk = [r["kappa"] for r in DT["kappa"] if r["kappa"] <= 0.5]
chk("III B", "kappa range in Fig. 1(b)", [0.003, 0.5], [min(kk), max(kk)], 0.012, "json")

# ---------------- Sec. II C: two lead implementations ----------------
def meso(L, k, K):
    def E(TL, TR):
        m = build_chain_with_leads(L, kappa=k, phi=np.pi / 2, K=K, W=4.0, TL=TL, TR=TR); c = lead_currents(m, ness_covariance(m)); return c
    f, r, i = E(1.0, 0.5), E(0.5, 1.0), E(0.5, 0.5)
    return -(f["R"]["E"] - i["R"]["E"]), -(r["L"]["E"] - i["L"]["E"])
CASES = ((12, 0.05), (20, 0.05), (20, 0.1), (20, 0.2), (40, 0.07), (40, 0.2))
for K in (80, 10):
    d = []
    for L, k in CASES:
        mf, mr = meso(L, k, K); wf, wr, *_ = delta_heat(chain_heff(L, kappa=k), 1.0, 0.5); d.append(100 * (mf / wf - 1))
        if mr > 1e-7:
            chk("III B", f"mesoscopic leads obey the law (K={K}, L={L}, kappa={k}): |ln ratio - law|", 1e-6, abs(np.log(mf / mr) - 2 * (L - 1) * inv_xi(k)), op="le", src="package")
    if K == 80:
        chk("II C", "K=80: largest |dQf(mesoscopic)/dQf(wide band) - 1| in % over six (L, kappa)", 6, max(abs(x) for x in d), op="le", src="package")
    else:
        chk("II C", "K=10: dQf below the wide-band value by (min, max) % over six (L, kappa)", [28, 38], [-max(d), -min(d)], 0.6, "package")

# ---------------- Sec. IV, Fig. 2 ----------------
nn = {(r["kappa"], r["t2"]): r for r in BN["nnn"]["sweep"]}
chk("IV", "thermal ln ratio kappa=0.2: t2=0 and t2=0.5", [15.65, 8.11], [nn[(0.2, 0.0)]["lnratio"], nn[(0.2, 0.5)]["lnratio"]], 0.005, "beyond_nn.json nnn")
chk("IV", "R stays above 0.97 (all t2, kappa)", 0.97, min(r["R"] for r in nn.values()), op="ge", src="beyond_nn.json nnn")
chk("IV", "real-axis winding number w(eps) (max |w| over all runs)", 0, max(abs(x) for r in nn.values() for x in r["w_real"]), 1e-6, "beyond_nn.json nnn w_real")
for k in (0.1, 0.2, 0.3):
    dv = max(np.max(np.abs(np.array(nn[(k, t)]["lnT"]) - np.array(nn[(k, t)]["nonbloch"]))) for t in (0.1, 0.3, 0.5))
    chk("IV", f"non-Bloch form: max |exact - 2(L-1)ln|b2 b3|| over energy, kappa={k} (the O(1) term)", 0, dv, op="info", src="beyond_nn.json nnn")

# ---------------- Sec. V ----------------
dc = BN["dt_convergence"]["16"]["0.0"]; g16 = BN["chi_convergence"]["16"]["gauss"]
chk("V B / Table I", "exact ln ratio, dQr (kappa=0.05, L=16)", [1.5003, 4.918e-3], [g16["ln"], g16["dQr"]], [5e-5, 5e-7], "beyond_nn.json chi_convergence")
chk("V B / Table I", "ln ratio: chi=128 dt=0.1/0.05/0.025, chi=256", [1.347, 1.412, 1.462, 1.5033], [dc["0.1"]["ln"], dc["0.05"]["ln"], dc["0.025"]["ln"], dc["chi256_dt0.05"]["ln"]], [5e-4, 5e-4, 5e-4, 5e-5], "beyond_nn.json dt_convergence")
chk("V B / Table I", "dQr bias in %", [16.8, 7.7, 2.6, 0.03], [100 * (dc[k]["dQr"] / g16["dQr"] - 1) for k in ("0.1", "0.05", "0.025", "chi256_dt0.05")], [0.05, 0.05, 0.05, 0.005], "beyond_nn.json")
chk("V B", "chi=256 ln-ratio error at kappa=0.05", 0.003, abs(dc["chi256_dt0.05"]["ln"] - g16["ln"]), 0.0005, "beyond_nn.json")
ck = BN["chi256_kappa0.1"]; cf = BN["chi256_final"]; gb = BN["gain_bath_V"]
chk("V B", "V=0 bias at kappa=0.1, L=12 and 16", [0.045, -0.038], [ck["12"]["0.0"] - ck["12"]["gauss"], ck["16"]["0.0"] - ck["16"]["gauss"]], 0.0006, "beyond_nn.json chi256_kappa0.1")


def run(L, V, cfg, kappa, tag):
    return json.load(open(os.path.join(ROOT, "results", "beyond_nn_runs", f"L{L}_V{V}_{cfg}_k{kappa}{tag}.json")))


def lnfluct(L, V, kappa, tag):
    f, r, i = (run(L, V, c, kappa, tag) for c in "f r iso".split())
    s = lambda x, k: np.std([t[k] for t in x["trace"][-5:]])
    dQf = -(f["final"]["ER"] - i["final"]["ER"]); dQr = -(r["final"]["EL"] - i["final"]["EL"])
    return float(np.hypot(np.hypot(s(f, "ER"), s(i, "ER")) / abs(dQf), np.hypot(s(r, "EL"), s(i, "EL")) / abs(dQr))), f, r, i


fl = {}
for (L, V, k, tag) in [(16, "0.0", 0.05, "_chi256_dt0.05"), (16, "2.0", 0.05, "_chi256_dt0.05"), (20, "2.0", 0.05, "_chi256_dt0.05"),
                       (12, "0.0", 0.1, "_chi256_dt0.05"), (12, "2.0", 0.1, "_chi256_dt0.05"), (16, "0.0", 0.1, "_chi256_dt0.05"), (16, "2.0", 0.1, "_chi256_dt0.05"),
                       (12, "0.0", 0.05, "_gain1_chi256_dt0.05"), (12, "2.0", 0.05, "_gain1_chi256_dt0.05")]:
    try:
        fl[(L, V, k, tag)] = lnfluct(L, V, k, tag)
    except FileNotFoundError as e:
        rows.append(f"| V B | run file missing: {os.path.basename(e.filename)} | | | | — |")
pure = lambda kk: max(v[0] for key, v in fl.items() if "gain" not in key[3] and key[2] == kk)
chk("V B", "time fluctuation of ln ratio over the last five windows, kappa=0.05 (max over chi=256 runs)", 0.011, pure(0.05), 0.0005, "beyond_nn_runs/*chi256*.json")
chk("V B", "same, kappa=0.1", 0.03, pure(0.1), 0.005, "beyond_nn_runs/*chi256*.json")
chk("V B / Fig. 3", "same, infinite-temperature-bath run at V=2J (also the stored error bar)", [0.08, 0.08], [fl[(12, "2.0", 0.05, "_gain1_chi256_dt0.05")][0], gb["gain_bath"]["err_V2"]], 0.005, "beyond_nn_runs/*gain1*.json, gain_bath_V")
chk("V C", "slope uncertainty at kappa=0.05 in % (time fluctuation of the two V=2J runs)", 4, 100 * np.hypot(fl[(16, "2.0", 0.05, "_chi256_dt0.05")][0], fl[(20, "2.0", 0.05, "_chi256_dt0.05")][0]) / 4 / cf["slope_gauss"], 0.7, "beyond_nn_runs")
chk("Fig. 3 caption", "error bar at kappa=0.1: hypot(V=0 bias, time fluctuation)", 0.06, np.hypot(0.045, pure(0.1)), 0.006, "see rows above")
# absolute-current error and dQr error at chi=128 and chi=256 (kappa=0.05, L=16, V=0)
for tag, lab in (("_dt0.05", "chi=128, dt=0.05"), ("_chi256_dt0.05", "chi=256, dt=0.05")):
    try:
        f, r, i = (run(16, "0.0", c, 0.05, tag) for c in "f r iso".split())
        ae = max(abs(x["final"][k] - x["gauss"][s]["E"]) for x in (f, r, i) for k, s in (("EL", "L"), ("ER", "R")))
        dqr = -(r["final"]["EL"] - i["final"]["EL"])
        qa, qd = ((7e-4, 3.8e-4) if "128" in lab else (6e-4, 1.5e-6))
        chk("V B", f"max absolute lead-current error, {lab} (kappa=0.05, L=16, V=0)", qa, ae, 0.5e-4, "beyond_nn_runs")
        chk("V B", f"|dQr - exact|, {lab}", qd, abs(dqr - g16["dQr"]), 0.05 * qd, "beyond_nn_runs")
    except FileNotFoundError as e:
        rows.append(f"| V B | run file missing: {os.path.basename(e.filename)} | | | | — |")
it = BN["interacting"]["sweep"]
chk("V B", "kappa=0.2: lower bounds on R (min, max over L=20,40 and V)", [0.68, 0.94], [min(r["R_lower"] for L in it for r in it[L]), max(r["R_lower"] for L in it for r in it[L])], 0.006, "beyond_nn.json interacting R_lower")
chk("V B", "kappa=0.2: dQf (1e-2) at V=0 and 2J, L=20 then L=40", [1.10, 0.63, 0.70, 0.33], [100 * it["20"][0]["dQf"], 100 * it["20"][-1]["dQf"], 100 * it["40"][0]["dQf"], 100 * it["40"][-1]["dQf"]], 0.006, "beyond_nn.json interacting")
chk("V C / Table II", "kappa=0.05: exact V=0 (L=16, 20)", [1.5003, 1.9004], [cf["ln_V0_gauss"]["16"], cf["ln_V0_gauss"]["20"]], 5e-5, "beyond_nn.json chi256_final")
chk("V C / Table II", "kappa=0.05: MPS V=2J (L=16, 20)", [1.083, 1.473], [cf["ln_V2"]["16"], cf["ln_V2"]["20"]], 5e-4, "chi256_final")
chk("V C / Table II", "kappa=0.05: offsets Delta_V (L=16, 20)", [-0.42, -0.43], [cf["offsets"]["16"], cf["offsets"]["20"]], 0.005, "chi256_final")
chk("V C / Table II", "kappa=0.05: slope V=2J and exact", [0.0975, 0.1000], [cf["slope_V2"], cf["slope_gauss"]], 5e-5, "chi256_final")
chk("V C", "slope deviation in %", 2.5, 100 * (1 - cf["slope_V2"] / cf["slope_gauss"]), 0.05, "chi256_final")
chk("V C / Table II", "kappa=0.1: exact, MPS V=0, MPS V=2J at L=12", [2.2018, 2.247, 1.704], [ck["12"]["gauss"], ck["12"]["0.0"], ck["12"]["2.0"]], [5e-5, 5e-4, 5e-4], "chi256_kappa0.1")
chk("V C / Table II", "kappa=0.1: exact, MPS V=0, MPS V=2J at L=16", [3.0025, 2.964, 2.448], [ck["16"]["gauss"], ck["16"]["0.0"], ck["16"]["2.0"]], [5e-5, 5e-4, 5e-4], "chi256_kappa0.1")
chk("V C / Table II", "kappa=0.1: Delta_V (L=12, 16)", [-0.54, -0.52], [ck["12"]["2.0"] - ck["12"]["0.0"], ck["16"]["2.0"] - ck["16"]["0.0"]], 0.005, "chi256_kappa0.1")
chk("V C / Table II", "kappa=0.1: slope (control variate) and exact", [0.207, 0.200], [(ck["16"]["cv"] - ck["12"]["cv"]) / 4, (ck["16"]["gauss"] - ck["12"]["gauss"]) / 4], 5e-4, "chi256_kappa0.1")
for k, q, e, gq in (("0.03", [0.055, 0.066, 0.046, 0.043], 0.013, 0.060), ("0.05", [0.096, 0.101, 0.125, 0.094], 0.019, 0.100)):
    s = BN["slopes"][k]
    chk("V C", f"chi=128 slopes at kappa={k} (V=0, 0.5, 1, 2)", q, [x["slope"] for x in s], 6e-4, "beyond_nn.json slopes")
    chk("V C", f"their uncertainty (max) and the exact slope, kappa={k}", [e, gq], [max(x["slope_err"] for x in s[:3]) if k == "0.03" else max(x["slope_err"] for x in s if x["V"] != 1.0), s[0]["gauss_slope"]], [0.0015, 5e-4], "slopes")
red = []
for key, L in (("interacting_k0.05", "20"), ("interacting_k0.05", "12"), ("interacting_k0.03", "20"), ("interacting_k0.1", "12")):
    rr = sorted(BN[key]["sweep"][L], key=lambda r: r["V"]); red.append(100 * (1 - rr[-1]["dQf"] / rr[0]["dQf"]))
chk("V C", "forward reduction at V=2J in % (min, max of the four series of Fig. 3d)", [25, 35], [min(red), max(red)], 1.0, "beyond_nn.json interacting_k*")
chk("V C / Table II", "bath run: MPS V=0, V=2J, Delta_V", [2.203, 1.662, -0.54], [gb["gain_bath"]["0.0"], gb["gain_bath"]["2.0"], gb["gain_bath"]["offset_V2"]], [5e-4, 5e-4, 0.005], "beyond_nn.json gain_bath_V")
chk("V C", "difference of Delta_V between bath and pure loss", 0.002, abs(gb["gain_bath"]["offset_V2"] - gb["pure_loss_k0p1"]["offset_V2"]), 0.0006, "gain_bath_V")
chk("V C", "V=0 MPS bias: bath and pure loss", [0.001, 0.045], [gb["gain_bath"]["0.0"] - ck["12"]["gauss"], ck["12"]["0.0"] - ck["12"]["gauss"]], 0.0006, "gain_bath_V, chi256_kappa0.1")

# ---------------- Sec. VI ----------------
kt = TB["kappa_Tb"]
sp = max(max(np.ptp([r[q] for r in kt if r["kappa"] == k]) for q in ("dQf", "dQr")) for k in (0.05, 0.1, 0.2, 0.3, 0.5))
chk("VI B", "spread of dQf and dQr over T_b = 0.25-2 (all kappa)", 2e-13, sp, op="le", src="thermal_bulk.json kappa_Tb")
chk("VI B", "spread of dQf and dQr over mu_b = 0..-30 (kappa=0.2)", 2e-13, max(np.ptp([r[q] for r in TB["mu_sweep"]]) for q in ("dQf", "dQr")), op="le", src="thermal_bulk.json mu_sweep")
chk("VI B", "T_b values of the sweep", [0.25, 0.5, 1.0, 2.0], sorted(set(r["bulk_bath"]["T"] for r in kt)), 0, "kappa_Tb")
chk("VI B", "mu_b values of the sweep", [0, -0.5, -1, -2, -4, -30], [r["bulk_bath"]["mu"] for r in TB["mu_sweep"]], 0, "mu_sweep")
dd = {(r["bulk_bath"]["K"], r["bulk_bath"]["W"]): r["lnratio"] for r in TB["discretisation"]}
chk("VI B", "ln ratio for (K_b, W_b) = (10,4), (10,6), (20,4), (20,6), (40,6); wide band", [14.23, 14.48, 14.68, 14.68, 14.95, 15.65], [dd[(10, 4.0)], dd[(10, 6.0)], dd[(20, 4.0)], dd[(20, 6.0)], dd[(40, 6.0)], TB["discretisation"][0]["nn_law"]], 0.006, "thermal_bulk.json discretisation")
c5 = [r for r in kt if r["bulk_bath"]["T"] == 0.5]
chk("VI C", "sigma_f for kappa = 0.05..0.5 (T_b=T_c)", [0.042, 0.044, 0.047, 0.049, 0.052], [r["f"]["sigma"] for r in c5], 6e-4, "kappa_Tb")
chk("VI C", "dQf (1e-2)", [1.64, 1.11, 0.66, 0.47, 0.29], [100 * r["dQf"] for r in c5], 0.006, "kappa_Tb")
chk("VI C", "sigma_f/dQf", [2.6, 4.0, 7.0, 10.4, 18.1], [r["f"]["sigma"] / r["dQf"] for r in c5], 0.06, "kappa_Tb")
chk("VI C", "heat into the bulk bath: kappa=0.05 and 0.5", [0.026, 0.049], [-c5[0]["f"]["E"]["bulk"], -c5[-1]["f"]["E"]["bulk"]], 6e-4, "kappa_Tb f.E.bulk")
chk("VI A", "isothermal sigma (max |sigma_iso| at T_b=T_c)", 1e-8, max(abs(r["iso"]["sigma"]) for r in c5), op="le", src="kappa_Tb")
chk("VI A", "sigma >= 0 in all configurations (min)", -1e-9, min(r[c]["sigma"] for r in kt + TB["mu_sweep"] for c in "f r iso".split()), op="ge", src="thermal_bulk.json")
ms = TB["mu_sweep"]
chk("VI C", "sigma_f over the mu_b sweep", [0.047, 0.12, 0.34, 1.15, 3.2, 25.5], [r["f"]["sigma"] for r in ms], [6e-4, 6e-3, 6e-3, 6e-3, 0.06, 0.06], "mu_sweep")
chk("VI C", "mean filling over the mu_b sweep", [0.50, 0.44, 0.38, 0.25, 0.12, 0.10], [r["f"]["n_mean"] for r in ms], 0.006, "mu_sweep")
chk("Fig. 4(b) caption", "sigma_f(mu_b=-30)/sigma_f(mu_b=0), quoted as a factor", 500, ms[-1]["f"]["sigma"] / ms[0]["f"]["sigma"], 60, "mu_sweep")
chk("Fig. 4(b) caption", "ln ratio quoted in the caption", 14.23, ms[0]["lnratio"], 0.006, "mu_sweep")


def sig(onsite, Tb, L=16, K=6, phi=np.pi / 3):
    m = build_chain_with_leads(L, kappa=0.2, phi=phi, K=K, TL=1.0, TR=0.5, onsite=onsite, bulk_bath=dict(K=6, W=4.0, T=Tb, mu=0.0))
    return entropy_production(m, lead_currents(m))


rng = np.random.default_rng(3); dis = rng.uniform(-1, 1, 16); sym = dis.copy(); sym[8:] = dis[:8][::-1]
chk("VI D", "mirror identity, clean chain: |sigma(T_b=T_c) - sigma(T_b=T_h)|", 1e-12, abs(sig(0.0, 0.5) - sig(0.0, 1.0)), op="le", src="package (as tests/test_thermal_bulk.py)")
chk("VI D", "same with on-site disorder of width W_d = J (16 sites)", 3e-4, abs(sig(dis / 2, 0.5) - sig(dis / 2, 1.0)), 0.5e-4, "package")
chk("VI D", "same with mirror-symmetric disorder", 1e-12, abs(sig(sym, 0.5) - sig(sym, 1.0)), op="le", src="package")

# ---------------- Sec. VII ----------------
on = DS["onsite"]
chk("VII A", "R for kappa = 0.05, 0.1, 0.2 (every W_d)", [0.9604, 0.9992, 1.0000], [np.mean([r["R_mean"] for r in on if r["kappa"] == k]) for k in (0.05, 0.1, 0.2)], 5e-5, "disorder.json onsite")
chk("VII A", "sample independence of the ln ratio (max deviation from the law)", 1e-13, max(r["ln_max_dev"] for r in on), op="le", src="onsite ln_max_dev")
m1 = {r["Wd"]: r["Qf_median"] for r in on if r["kappa"] == 0.1}
chk("VII A", "median dQf (1e-2) at kappa=0.1 for W_d = 0, 0.5, 1, 2, 3, 4", [1.62, 1.26, 0.61, 0.051, 0.0014, 0.000016], [100 * m1[w] for w in (0, 0.5, 1.0, 2.0, 3.0, 4.0)], [0.006, 0.006, 0.006, 0.0006, 0.00006, 0.000001], "onsite Qf_median")
chk("VII A", "largest W_d with median dQf >= 1e-3 for all three kappa", 1.5, max(w for w in sorted(set(r["Wd"] for r in on)) if all(r["Qf_median"] >= 1e-3 for r in on if r["Wd"] == w)), 0, "onsite Qf_median")
ld = {(r["kappa"], r["Wd"], r["L"]): r["lnQf_mean"] for r in DS["Ldep"]}
chk("VII A / Table III", "forward decay rates (rows kappa = 0, 0.05, 0.2; columns W_d = 0.5, 1, 2)", [0.006, 0.020, 0.055, 0.017, 0.034, 0.090, 0.027, 0.044, 0.107],
    [(ld[(k, w, 10)] - ld[(k, w, 80)]) / 70 for k in (0.0, 0.05, 0.2) for w in (0.5, 1.0, 2.0)], 6e-4, "disorder.json Ldep")
bd = {(r["kind"], r["L"]): r for r in DS["bond"]}
chk("VII B", "sum rule: max deviation over samples", 1e-13, max(r["max_dev_from_sum_rule"] for r in DS["bond"]), op="le", src="disorder.json bond")
chk("VII B", "kappa spread: means (L=20,40,80)", [3.75, 7.85, 15.83], [bd[("kappa", L)]["ln_mean"] for L in (20, 40, 80)], 0.006, "bond")
chk("VII B", "clean values", [3.80, 7.81, 15.81], [bd[("kappa", L)]["clean"] for L in (20, 40, 80)], 0.006, "bond")
chk("VII B", "kappa spread: standard deviations", [0.43, 0.55, 0.88], [bd[("kappa", L)]["ln_std"] for L in (20, 40, 80)], 0.006, "bond")
chk("VII B", "phase spread: mean at L=40, clean x sin(1), sin(1)", [6.56, 6.57, 0.84], [bd[("phi", 40)]["ln_mean"], bd[("phi", 40)]["clean"] * np.sin(1), np.sin(1)], 0.006, "bond")
chk("VII B", "phase spread: standard deviations (min, max)", [0.13, 0.25], [min(bd[("phi", L)]["ln_std"] for L in (20, 40, 80)), max(bd[("phi", L)]["ln_std"] for L in (20, 40, 80))], 0.006, "bond")
chk("VII B", "random chirality: means", [-0.07, -0.04, 0.16], [bd[("sign", L)]["ln_mean"] for L in (20, 40, 80)], 0.006, "bond")
chk("VII B", "random chirality: standard deviations", [0.91, 1.34, 1.77], [bd[("sign", L)]["ln_std"] for L in (20, 40, 80)], 0.006, "bond")
nd = {(r["t2"], r["Wd"]): r for r in DS["nnn"]}
chk("VII C", "t2=0.3, kappa=0.2: means for W_d = 0.5, 1, 2 and clean", [10.12, 9.71, 9.42, 10.31], [nd[(0.3, w)]["ln_mean"] for w in (0.5, 1.0, 2.0, 0)], 0.006, "disorder.json nnn")
chk("VII C", "their standard deviations", [0.11, 0.19, 0.30], [nd[(0.3, w)]["ln_std"] for w in (0.5, 1.0, 2.0)], 0.006, "nnn")

# ---------------- Sec. VIII ----------------
need = lambda k: int(np.ceil(np.arctanh(0.99) / inv_xi(k) + 1))
chk("VIII / Table IV", "L needed for R>0.99 at kappa/J = 0.2 and 0.05", [15, 54], [need(0.2), need(0.05)], 0, "law")
chk("VIII", "kappa/2pi (Hz) for J/2pi = 0.5 kHz and kappa/J = 0.05, 0.2", [25, 100], [500 * 0.05, 500 * 0.2], 0, "arithmetic")
chk("VIII / Table IV", "kappa/2pi = 4 g^2/Gamma_a (MHz) for g = 1.1, 2.2 MHz, Gamma_a = 10 MHz", [0.5, 2.0], [4 * 1.1 ** 2 / 10, 4 * 2.2 ** 2 / 10], 0.07, "arithmetic")
chk("VIII / Table IV", "kappa/J for kappa/2pi = 0.5, 2 MHz and J/2pi = 10 MHz", [0.05, 0.2], [0.5 / 10, 2 / 10], 0, "arithmetic")

# ---------------- the quoted values above must be the ones printed in the manuscript ----------------
TEX = open(os.path.join(ROOT, "paper", "prb", "main.tex")).read()
SNIPPETS = [r"(L-1)/\xi>2.65", r"$\R=0.992$", r"$45\%$ of its lossless value", r"agree within $6\%$", r"$28$--$38\%$", r"$15.65\to8.11$", r"stays above $0.97$",
            r"1.347 & 1.412 & 1.462 & 1.5033", r"$+16.8\%$ & $+7.7\%$ & $+2.6\%$ & $+0.03\%$", r"$\ln(\dQ_f/\dQ_r)=1.5003$", r"4.918\times10^{-3}J^2$",
            r"$+0.045$ ($L=12$) and $-0.038$ ($L=16$)", r"at most $0.011$", r"$0.03$ at $\kappa=0.1J$, and $0.08$", r"$7\times10^{-4}J^2$", r"$6\times10^{-4}J^2$",
            r"$3.8\times10^{-4}J^2$ to $1.5\times10^{-6}J^2$", r"$\R\ge0.68$--$0.94$", r"from $1.10$ to $0.63\times10^{-2}J^2$", r"from $0.70$ to $0.33\times10^{-2}J^2$",
            r"slope is $0.0975$ against the free-fermion $0.1000$", r"is $-0.42$ at both", r"slope is $0.207$ against $0.200$", r"$\Delta_V=-0.54$ and $-0.52$",
            r"0.055$, $0.066$, $0.046$, $0.043$ ($\pm0.013$; exact $0.060$)", r"$0.096$, $0.101$, $0.125$, $0.094$ ($\pm0.019$; exact $0.100$)", r"$25$--$35\%$",
            r"$\Delta_V=-0.54\pm0.08$", r"($+0.001$ against $+0.045$)", r"0.05 & 16 & 1.5003 & 1.5033 & 1.083 & $-0.42$", r"0.05 & 20 & 1.9004 & --- & 1.473 & $-0.43$",
            r"0.1 & 12 & 2.2018 & 2.247 & 1.704 & $-0.54$", r"0.1 & 16 & 3.0025 & 2.964 & 2.448 & $-0.52$", r"0.1 (bath) & 12 & 2.2018 & 2.203 & 1.662 & $-0.54$",
            r"$2\times10^{-13}J^2$", r"14.23$, $14.48$, $14.68$ and $14.95$", r"wide-band value $15.65$", r"0.042$, $0.044$, $0.047$, $0.049$ and $0.052\,J$",
            r"1.64$, $1.11$, $0.66$, $0.47$ and $0.29\times10^{-2}J^2$", r"$2.6$, $4.0$, $7.0$, $10.4$ and $18.1\,J^{-1}$", r"from $0.026$ to $0.049\,J^2$",
            r"0.047$, $0.12$, $0.34$, $1.15$, $3.2$ and $25.5\,J$", r"$0.50$, $0.44$, $0.38$, $0.25$, $0.12$ and $0.10$", r"about $500$ times", r"($\ln=14.23$)",
            r"$3\times10^{-4}J$ for $W_d=J$", r"$\R=0.9604$, $0.9992$ and $1.0000$", r"$1.62$, $1.26$, $0.61$, $0.051$, $0.0014$ and $0.000016\times10^{-2}J^2$", r"$W_d\lesssim1.5J$",
            r"$\kappa=0$ & 0.006 & 0.020 & 0.055", r"$\kappa=0.05J$ & 0.017 & 0.034 & 0.090", r"$\kappa=0.2J$ & 0.027 & 0.044 & 0.107",
            r"($3.75$, $7.85$ and $15.83$ against the clean $3.80$, $7.81$ and $15.81$", r"$0.43$, $0.55$ and $0.88$", r"($6.56$ against $7.81\times0.84=6.57$", r"$0.13$--$0.25$",
            r"$-0.07$, $-0.04$ and $0.16$", r"$0.91$, $1.34$ and $1.77$", r"10.12\pm0.11$, $9.71\pm0.19$ and $9.42\pm0.30$", r"clean value $10.31$",
            r"15--54 & 15--54", r"$g=2\pi\times1.1$--$2.2$~MHz", r"$\kappa=2\pi\times25$--$100$~Hz", r"about $4\%$ of the slope"]
missing = [x for x in SNIPPETS if x not in TEX]
for x in missing:
    rows.append(f"| text | quoted value not found verbatim in main.tex: `{x}` | | | paper/prb/main.tex | **FAIL** |")
fails += len(missing)
rows.append(f"| text | {len(SNIPPETS) - len(missing)} of {len(SNIPPETS)} quoted strings found verbatim in main.tex | | | paper/prb/main.tex | {'OK' if not missing else '**FAIL**'} |")

print("| Section | Quantity | Quoted in the manuscript | Recomputed | Source | Result |\n|---|---|---|---|---|---|")
print("\n".join(rows))
print(f"\n{len(rows)} rows, {fails} failed.")
sys.exit(1 if fails else 0)
