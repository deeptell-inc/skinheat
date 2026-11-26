| Section | Quantity | Quoted in the manuscript | Recomputed | Source | Result |
|---|---|---|---|---|---|
| III B | threshold (L-1)/xi for R>0.99 | 2.65 | 2.6467 | atanh(0.99) | OK |
| III B | R at kappa=0.07J, L=40 | 0.992 | 0.99179 | deltaT_rectification.json kappa[] | OK |
| III B | dQf(kappa=0.07)/dQf(kappa=0) in % | 45 | 44.676 | json + chain_heff(kappa=0) | OK |
| III B | Fig. 1(b): max |ln ratio - law| for predicted ln ratio <= 48 | 2e-11 | 1.0481e-11 | deltaT_rectification.json | OK |
| III B | Fig. 1(b): |ln ratio - law| at the point with predicted ln ratio 72 | 2e-06 | 1.625e-06 | deltaT_rectification.json | OK |
| III B | Fig. 1(b): number of stored points beyond 75 (omitted from the figure) and its predicted ln ratio | 1, 96 | 1, 96.125 | deltaT_rectification.json | OK |
| III B | L range in Fig. 1(b) | 5, 160 | 5, 160 | json | OK |
| III B | kappa range in Fig. 1(b) | 0.003, 0.5 | 0.0031623, 0.48919 | json | OK |
| III B | mesoscopic leads obey the law (K=80, L=12, kappa=0.05): |ln ratio - law| | 1e-06 | 4.8901e-12 | package | OK |
| III B | mesoscopic leads obey the law (K=80, L=20, kappa=0.05): |ln ratio - law| | 1e-06 | 2.651e-11 | package | OK |
| III B | mesoscopic leads obey the law (K=80, L=20, kappa=0.1): |ln ratio - law| | 1e-06 | 6.1788e-11 | package | OK |
| III B | mesoscopic leads obey the law (K=80, L=20, kappa=0.2): |ln ratio - law| | 1e-06 | 8.7149e-09 | package | OK |
| III B | mesoscopic leads obey the law (K=80, L=40, kappa=0.07): |ln ratio - law| | 1e-06 | 1.0691e-09 | package | OK |
| II C | K=80: largest |dQf(mesoscopic)/dQf(wide band) - 1| in % over six (L, kappa) | 6 | 5.826 | package | OK |
| III B | mesoscopic leads obey the law (K=10, L=12, kappa=0.05): |ln ratio - law| | 1e-06 | 2.4836e-12 | package | OK |
| III B | mesoscopic leads obey the law (K=10, L=20, kappa=0.05): |ln ratio - law| | 1e-06 | 1.9385e-11 | package | OK |
| III B | mesoscopic leads obey the law (K=10, L=20, kappa=0.1): |ln ratio - law| | 1e-06 | 2.0151e-10 | package | OK |
| III B | mesoscopic leads obey the law (K=10, L=20, kappa=0.2): |ln ratio - law| | 1e-06 | 1.4064e-08 | package | OK |
| III B | mesoscopic leads obey the law (K=10, L=40, kappa=0.07): |ln ratio - law| | 1e-06 | 5.3342e-10 | package | OK |
| II C | K=10: dQf below the wide-band value by (min, max) % over six (L, kappa) | 28, 38 | 28.216, 37.947 | package | OK |
| IV | thermal ln ratio kappa=0.2: t2=0 and t2=0.5 | 15.65, 8.11 | 15.652, 8.1062 | beyond_nn.json nnn | OK |
| IV | R stays above 0.97 (all t2, kappa) | 0.97 | 0.97132 | beyond_nn.json nnn | OK |
| IV | real-axis winding number w(eps) (max |w| over all runs) | 0 | 0 | beyond_nn.json nnn w_real | OK |
| IV | non-Bloch form: max |exact - 2(L-1)ln|b2 b3|| over energy, kappa=0.1 (the O(1) term) | 0 | 3.8265 | beyond_nn.json nnn | — |
| IV | non-Bloch form: max |exact - 2(L-1)ln|b2 b3|| over energy, kappa=0.2 (the O(1) term) | 0 | 2.1166 | beyond_nn.json nnn | — |
| IV | non-Bloch form: max |exact - 2(L-1)ln|b2 b3|| over energy, kappa=0.3 (the O(1) term) | 0 | 3.0516 | beyond_nn.json nnn | — |
| V B / Table I | exact ln ratio, dQr (kappa=0.05, L=16) | 1.5003, 0.004918 | 1.5003, 0.0049178 | beyond_nn.json chi_convergence | OK |
| V B / Table I | ln ratio: chi=128 dt=0.1/0.05/0.025, chi=256 | 1.347, 1.412, 1.462, 1.5033 | 1.347, 1.4122, 1.4616, 1.5033 | beyond_nn.json dt_convergence | OK |
| V B / Table I | dQr bias in % | 16.8, 7.7, 2.6, 0.03 | 16.802, 7.7012, 2.636, 0.030374 | beyond_nn.json | OK |
| V B | chi=256 ln-ratio error at kappa=0.05 | 0.003 | 0.0029498 | beyond_nn.json | OK |
| V B | V=0 bias at kappa=0.1, L=12 and 16 | 0.045, -0.038 | 0.045441, -0.038162 | beyond_nn.json chi256_kappa0.1 | OK |
| V B | time fluctuation of ln ratio over the last five windows, kappa=0.05 (max over chi=256 runs) | 0.011 | 0.010917 | beyond_nn_runs/*chi256*.json | OK |
| V B | same, kappa=0.1 | 0.03 | 0.031223 | beyond_nn_runs/*chi256*.json | OK |
| V B / Fig. 3 | same, infinite-temperature-bath run at V=2J (also the stored error bar) | 0.08, 0.08 | 0.081996, 0.08 | beyond_nn_runs/*gain1*.json, gain_bath_V | OK |
| V C | slope uncertainty at kappa=0.05 in % (time fluctuation of the two V=2J runs) | 4 | 3.4347 | beyond_nn_runs | OK |
| Fig. 3 caption | error bar at kappa=0.1: hypot(V=0 bias, time fluctuation) | 0.06 | 0.054771 | see rows above | OK |
| V B | max absolute lead-current error, chi=128, dt=0.05 (kappa=0.05, L=16, V=0) | 0.0007 | 0.00072636 | beyond_nn_runs | OK |
| V B | |dQr - exact|, chi=128, dt=0.05 | 0.00038 | 0.00037873 | beyond_nn_runs | OK |
| V B | max absolute lead-current error, chi=256, dt=0.05 (kappa=0.05, L=16, V=0) | 0.0006 | 0.00063307 | beyond_nn_runs | OK |
| V B | |dQr - exact|, chi=256, dt=0.05 | 1.5e-06 | 1.4938e-06 | beyond_nn_runs | OK |
| V B | kappa=0.2: lower bounds on R (min, max over L=20,40 and V) | 0.68, 0.94 | 0.68152, 0.93865 | beyond_nn.json interacting R_lower | OK |
| V B | kappa=0.2: dQf (1e-2) at V=0 and 2J, L=20 then L=40 | 1.1, 0.63, 0.7, 0.33 | 1.0958, 0.634, 0.69942, 0.33071 | beyond_nn.json interacting | OK |
| V C / Table II | kappa=0.05: exact V=0 (L=16, 20) | 1.5003, 1.9004 | 1.5003, 1.9004 | beyond_nn.json chi256_final | OK |
| V C / Table II | kappa=0.05: MPS V=2J (L=16, 20) | 1.083, 1.473 | 1.0832, 1.4733 | chi256_final | OK |
| V C / Table II | kappa=0.05: offsets Delta_V (L=16, 20) | -0.42, -0.43 | -0.41708, -0.42714 | chi256_final | OK |
| V C / Table II | kappa=0.05: slope V=2J and exact | 0.0975, 0.1 | 0.097506, 0.1 | chi256_final | OK |
| V C | slope deviation in % | 2.5 | 2.4939 | chi256_final | OK |
| V C / Table II | kappa=0.1: exact, MPS V=0, MPS V=2J at L=12 | 2.2018, 2.247, 1.704 | 2.2018, 2.2473, 1.704 | chi256_kappa0.1 | OK |
| V C / Table II | kappa=0.1: exact, MPS V=0, MPS V=2J at L=16 | 3.0025, 2.964, 2.448 | 3.0025, 2.9643, 2.4477 | chi256_kappa0.1 | OK |
| V C / Table II | kappa=0.1: Delta_V (L=12, 16) | -0.54, -0.52 | -0.54323, -0.51668 | chi256_kappa0.1 | OK |
| V C / Table II | kappa=0.1: slope (control variate) and exact | 0.207, 0.2 | 0.2068, 0.20017 | chi256_kappa0.1 | OK |
| V C | chi=128 slopes at kappa=0.03 (V=0, 0.5, 1, 2) | 0.055, 0.066, 0.046, 0.043 | 0.055494, 0.065612, 0.04643, 0.04327 | beyond_nn.json slopes | OK |
| V C | their uncertainty (max) and the exact slope, kappa=0.03 | 0.013, 0.06 | 0.013427, 0.060005 | slopes | OK |
| V C | chi=128 slopes at kappa=0.05 (V=0, 0.5, 1, 2) | 0.096, 0.101, 0.125, 0.094 | 0.096204, 0.10109, 0.12538, 0.093809 | beyond_nn.json slopes | OK |
| V C | their uncertainty (max) and the exact slope, kappa=0.05 | 0.019, 0.1 | 0.018876, 0.10002 | slopes | OK |
| V C | forward reduction at V=2J in % (min, max of the four series of Fig. 3d) | 25, 35 | 24.817, 34.383 | beyond_nn.json interacting_k* | OK |
| V C / Table II | bath run: MPS V=0, V=2J, Delta_V | 2.203, 1.662, -0.54 | 2.203, 1.662, -0.541 | beyond_nn.json gain_bath_V | OK |
| V C | difference of Delta_V between bath and pure loss | 0.002 | 0.002 | gain_bath_V | OK |
| V C | V=0 MPS bias: bath and pure loss | 0.001, 0.045 | 0.0011639, 0.045441 | gain_bath_V, chi256_kappa0.1 | OK |
| VI B | spread of dQf and dQr over T_b = 0.25-2 (all kappa) | 2e-13 | 8.412e-14 | thermal_bulk.json kappa_Tb | OK |
| VI B | spread of dQf and dQr over mu_b = 0..-30 (kappa=0.2) | 2e-13 | 1.8902e-13 | thermal_bulk.json mu_sweep | OK |
| VI B | T_b values of the sweep | 0.25, 0.5, 1, 2 | 0.25, 0.5, 1, 2 | kappa_Tb | OK |
| VI B | mu_b values of the sweep | 0, -0.5, -1, -2, -4, -30 | 0, -0.5, -1, -2, -4, -30 | mu_sweep | OK |
| VI B | ln ratio for (K_b, W_b) = (10,4), (10,6), (20,4), (20,6), (40,6); wide band | 14.23, 14.48, 14.68, 14.68, 14.95, 15.65 | 14.232, 14.475, 14.675, 14.682, 14.946, 15.652 | thermal_bulk.json discretisation | OK |
| VI C | sigma_f for kappa = 0.05..0.5 (T_b=T_c) | 0.042, 0.044, 0.047, 0.049, 0.052 | 0.042427, 0.044012, 0.046594, 0.048716, 0.052032 | kappa_Tb | OK |
| VI C | dQf (1e-2) | 1.64, 1.11, 0.66, 0.47, 0.29 | 1.6372, 1.1082, 0.6619, 0.46705, 0.28727 | kappa_Tb | OK |
| VI C | sigma_f/dQf | 2.6, 4, 7, 10.4, 18.1 | 2.5914, 3.9713, 7.0395, 10.431, 18.113 | kappa_Tb | OK |
| VI C | heat into the bulk bath: kappa=0.05 and 0.5 | 0.026, 0.049 | 0.026055, 0.049159 | kappa_Tb f.E.bulk | OK |
| VI A | isothermal sigma (max |sigma_iso| at T_b=T_c) | 1e-08 | 1.7523e-12 | kappa_Tb | OK |
| VI A | sigma >= 0 in all configurations (min) | -1e-09 | -1.7523e-12 | thermal_bulk.json | OK |
| VI C | sigma_f over the mu_b sweep | 0.047, 0.12, 0.34, 1.15, 3.2, 25.5 | 0.046594, 0.12141, 0.3431, 1.1498, 3.1786, 25.54 | mu_sweep | OK |
| VI C | mean filling over the mu_b sweep | 0.5, 0.44, 0.38, 0.25, 0.12, 0.1 | 0.5, 0.43966, 0.37706, 0.25178, 0.1184, 0.10087 | mu_sweep | OK |
| Fig. 4(b) caption | sigma_f(mu_b=-30)/sigma_f(mu_b=0), quoted as a factor | 500 | 548.14 | mu_sweep | OK |
| Fig. 4(b) caption | ln ratio quoted in the caption | 14.23 | 14.232 | mu_sweep | OK |
| VI D | mirror identity, clean chain: |sigma(T_b=T_c) - sigma(T_b=T_h)| | 1e-12 | 1.435e-14 | package (as tests/test_thermal_bulk.py) | OK |
| VI D | same with on-site disorder of width W_d = J (16 sites) | 0.0003 | 0.00031534 | package | OK |
| VI D | same with mirror-symmetric disorder | 1e-12 | 1.2209e-13 | package | OK |
| VII A | R for kappa = 0.05, 0.1, 0.2 (every W_d) | 0.9604, 0.9992, 1 | 0.96035, 0.99919, 1 | disorder.json onsite | OK |
| VII A | sample independence of the ln ratio (max deviation from the law) | 1e-14 | 7.1054e-15 | onsite ln_max_dev | OK |
| VII A | median dQf (1e-2) at kappa=0.1 for W_d = 0, 0.5, 1, 2, 3, 4 | 1.62, 1.26, 0.61, 0.051, 0.0014, 1.6e-05 | 1.6212, 1.2562, 0.6088, 0.051479, 0.0013976, 1.6441e-05 | onsite Qf_median | OK |
| VII A | largest W_d with median dQf >= 1e-3 for all three kappa | 1.5 | 1.5 | onsite Qf_median | OK |
| VII A / Table III | forward decay rates (rows kappa = 0, 0.05, 0.2; columns W_d = 0.5, 1, 2) | 0.006, 0.02, 0.055, 0.017, 0.034, 0.09, 0.027, 0.044, 0.107 | 0.0061991, 0.019803, 0.055182, 0.017073, 0.033572, 0.090052, 0.027028, 0.044027, 0.10718 | disorder.json Ldep | OK |
| VII B | sum rule: max deviation over samples | 3e-14 | 2.1316e-14 | disorder.json bond | OK |
| VII B | kappa spread: means (L=20,40,80) | 3.75, 7.85, 15.83 | 3.7535, 7.8514, 15.826 | bond | OK |
| VII B | clean values | 3.8, 7.81, 15.81 | 3.8032, 7.8065, 15.813 | bond | OK |
| VII B | kappa spread: standard deviations | 0.43, 0.55, 0.88 | 0.42545, 0.54899, 0.87564 | bond | OK |
| VII B | phase spread: mean at L=40, clean x sin(1), sin(1) | 6.56, 6.57, 0.84 | 6.5606, 6.569, 0.84147 | bond | OK |
| VII B | phase spread: standard deviations (min, max) | 0.13, 0.25 | 0.12762, 0.24557 | bond | OK |
| VII B | random chirality: means | -0.07, -0.04, 0.16 | -0.07206, -0.03603, 0.16414 | bond | OK |
| VII B | random chirality: standard deviations | 0.91, 1.34, 1.77 | 0.91093, 1.3435, 1.7661 | bond | OK |
| VII C | t2=0.3, kappa=0.2: means for W_d = 0.5, 1, 2 and clean | 10.12, 9.71, 9.42, 10.31 | 10.124, 9.7065, 9.4235, 10.313 | disorder.json nnn | OK |
| VII C | their standard deviations | 0.11, 0.19, 0.3 | 0.11009, 0.18705, 0.29606 | nnn | OK |
| VIII / Table IV | L needed for R>0.99 at kappa/J = 0.2 and 0.05 | 15, 54 | 15, 54 | law | OK |
| VIII | kappa/2pi (Hz) for J/2pi = 0.5 kHz and kappa/J = 0.05, 0.2 | 25, 100 | 25, 100 | arithmetic | OK |
| VIII / Table IV | kappa/2pi = 4 g^2/Gamma_a (MHz) for g = 1.1, 2.2 MHz, Gamma_a = 10 MHz | 0.5, 2 | 0.484, 1.936 | arithmetic | OK |
| VIII / Table IV | kappa/J for kappa/2pi = 0.5, 2 MHz and J/2pi = 10 MHz | 0.05, 0.2 | 0.05, 0.2 | arithmetic | OK |
| text | 70 of 70 quoted strings found verbatim in main.tex | | | paper/prb/main.tex | OK |

97 rows, 0 failed.
