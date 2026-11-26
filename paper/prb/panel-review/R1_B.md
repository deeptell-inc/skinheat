# Round-1 record of Panelist B (reproduction from raw run files and the project package)
B1 UNSUPPORTED: "xi not renormalised within about 4% at kappa=0.05". 4% is only the time scatter of two V=2 runs. Omitted: no chi=256 V=0 run at L=20; bias at kappa=0.1 flips sign (+0.045 at L=12, -0.038 at L=16 = 10% of slope); L=20 V=2 runs and f/iso at L=16 V=2 stopped at tmax=60 unconverged (median dE 2-6e-5 vs criterion 1e-5); L=16 V=2 ln ratio still drifting by 0.032 over the last 20/J; slope 0.0975 (mixed finals), 0.0966 (last-5 means), 0.1017 (time-aligned finals); V=2 L=16 value 1.083 (chi=256), 1.097/1.101/1.148 (chi=128, dt=0.025/0.05/0.1), 1.173 (chi=96). Defensible: about +-10%. kappa=0.1: raw V=2 slope 0.186 (-7%), V=0 MPS slope 0.179 (-10% from exact); control-variate slope 0.207 assumes V-independent bias. Fig.3(c) error bar 0.04 is hard-coded.
B2 ERROR: "truncation error accumulates with the number of time steps, so it is not reduced by a smaller dt; it is reduced only by chi" contradicts Table I: at chi=128 the dQr bias falls 16.8% -> 7.7% -> 2.6% for dt 0.1 -> 0.05 -> 0.025; absolute current error 7.3e-4 (chi=128) vs 6.3e-4 (chi=256), non-monotonic in dt (7.8, 7.3, 11.6 e-4).
B3 IMPRECISE: "dQr within 0.03%, error 1.5e-6" is a single snapshot; last-5-window scatter of dQr is 1.6e-5 (0.3%); ln ranges 1.5025-1.5096 over the last 10/J.
B4 Table II reproduces; but Delta_V is defined against exact V=0 at kappa=0.05 and against MPS V=0 at kappa=0.1; text says "-0.42 at both L=16 and L=20" while table/data give -0.43 at L=20; table mixes runs stopped at different times.
B5 ERROR (minor): bath V=0 bias is +0.0017 -> +0.002, not +0.001.
B6 ERROR: chi=128 slope uncertainties +-0.013 and +-0.019 are not the maxima (0.0142, 0.0217); fits weight by a systematic (V=0 bias) treated as random; three of eight slopes deviate by more than 20% from exact (0.125 +25%, 0.046 -23%, 0.043 -28%); all 24 chi=128 groups unconverged; these series constrain the slope to +-25-30%.
B7 IMPRECISE: "forward throughput reduced by 25-35% at V=2J": plotted series 24.8-34.4%; unplotted: 19.4% (kappa=0.03, L=12), 37.6% (chi=256, kappa=0.1, L=16), 42-53% (kappa=0.2), 76% (bath run).
B8: kappa=0.2 runs used chi=64, dt=0.1, variational compression (not disclosed); V=0 throughput 4% (L=20) and 14% (L=40) above exact; measured dQr negative in three cases, so "R >= 0.68-0.94" is not a bound there; mismatched stop times.
B9 REPRODUCED: time fluctuations 0.011/0.03/0.08.
B10: Fig. 3 error bar 0.02 at kappa=0.05 is not what the stated recipe gives (hypot(0.003,0.011)=0.011); 0.06, 0.04 hard-coded; Fig. 3(b) reuses bars although the bias should cancel there.
B11: convergence flag set on a single quiet window; 59 of 209 run files flagged.
B12: "agrees with the Lyapunov solution to 2e-3": the test tolerance is 2e-3*max(1,|E|/0.05) = 4.2-4.5e-3.
B13, B14 REPRODUCED: 45%, R=0.992; law precision 2e-11 / 2e-6.
B15 IMPRECISE: K=80 within 6% fails for L=12, kappa=0.03 (6.26%); K=160 gives +1.8 to +9.1% (worse than K=80), so mesoscopic leads do not converge to wide band with K at W=4J; with K=10, L=40, kappa=0.2 the law deviates by 1e-4.
B16 IMPRECISE: "decays as 1/L for any kappa" tested for a single kappa (R = -0.070, -0.039, -0.019 at L=10, 20, 40).
B17 ERROR: "Eq.(13) reproduces the exact transmissions for all t2" but the highlighted sign change is never reproduced: non-Bloch curve minimum +0.070; exact negative at several (kappa,t2); e.g. kappa=0.2, t2=0.3, eps=-3.4: exact -0.30 vs non-Bloch +0.51; max mismatch 1.0-3.8. Negative window at kappa=0.1, t2=0.3 spans eps in [-4,-2.95] (not narrow); no sign change at kappa=0.3, t2=0.3; for t2=0.1 the peak lies 0.6-0.75J below the Hermitian band bottom.
B19 ERROR: "sigma_f/dQf ∝ kappa ∝ 1/xi": sigma_f/(dQf kappa) = 51.8, 39.7, 35.2, 34.8, 36.2; affine fit 34.6 kappa + 0.47; dotted line overshoots the kappa=0.5 point by 43%; kappa->0 limit is 1/T_c-1/T_h = 1/J.
B20: bath independence 2e-13 absolute is meaningless at kappa>=0.3 where dQr ~1e-12-1e-14 (ln ratios scatter 21.213-21.233; nan at 0.5).
B22: "sigma=0 when all reservoirs share one temperature" also requires equal chemical potentials (0.329 at mu_b=-1).
B24 IMPRECISE: mirror-identity disorder number 3e-4 is one seed; seeds 0-11: median 9e-4, range 8.8e-5 to 2.8e-3. "holds for any t2" not tested numerically (builder has no t2 argument).
B25 IMPRECISE: "nearly additive": excess over clean+Anderson is +-2% (W_d=0.5J), 7-10% (J), 37-40% (2J).
B26, B28 REPRODUCED: on-site invariance, sum rule, bond statistics, Table IV arithmetic. "10^2 samples" is 900 in total; tests assert only 1e-8-1e-9, finer precision rests on data.
B27: Fig. 4(c) uses 100 samples; "60-100" applies to Table III and VII C.
Checker defects: filters hide failing entries (slope-uncertainty maxima; mesoscopic-law check skipped for small dQr); case selection; hand-typed JSON values; loose tolerances; many claims never checked.
Top three: B1, B6, B17+B19 (then B2).
