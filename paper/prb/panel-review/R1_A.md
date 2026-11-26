# Round-1 record of Panelist A (independent first-principles code; did not use the project package)
A1 VERIFIED: covariance equation Eq.(2) incl. transposes (brute force 2-3 modes, error 1e-15).
A2 VERIFIED: h_eff elements, t_+ superdiagonal, t_- subdiagonal; phi=pi/2: t_+=-0.9, t_-=-1.1 (kappa=0.2); right eigenvectors localise to the right; G_{L1} is left-to-right; directions consistent.
A3 IMPRECISE: "right eigenvectors localised with inverse localisation length xi^-1 = ln|t_-/t_+|": amplitude envelope is |t_-/t_+|^{j/2}, so xi is the intensity length; amplitude length is 2 xi. Law unaffected.
A4 VERIFIED: corner identity, law, R=tanh, atanh(0.99)=2.6467, L=15 (kappa=0.2) and 54 (kappa=0.05).
A5 VERIFIED: Landauer form Eq.(8) holds for the background-subtracted response (not for total current) with bulk loss; mu=0 so heat = energy.
A6 IMPRECISE: bosonic extension. With the band at negative energies and mu=0 the Bose function is negative (n_B(-J;T=J)=-1.58): no stable bosonic equilibrium; needs positive carrier frequency or mu below the band; heat weight is omega_0+eps for photons. Constant ratio survives.
A7 VERIFIED: characteristic equation Eq.(12) placement of t_+/t_-; root product 1; t2=0 limit product t_-/t_+.
A8 ERROR: text before Eq.(13) says left-to-right amplitude ~|beta_3|^L and right-to-left ~|beta_2|^-L. Numerics (t2=0.3, eps=0.7, kappa=0.2; roots 0.2296, 0.9904, 1.1435, 3.8449): dln|G_L1|/dL=-0.009653=ln|beta_2|, dln|G_1L|/dL=-0.134124=-ln|beta_3|. So |G_L1|~|beta_2|^L, |G_1L|~|beta_3|^-L; the ratio exponent Eq.(13) is still correct.
A9 VERIFIED with caveat: Eq.(13) residual D_L=ln(T_RL/T_LR)-2(L-1)ln|b2 b3| max over eps in [-4,4]: t2=0.1: 1.86 (L=40), 1.83 (80), 1.83 (160); t2=0.3: 1.76, 0.49, 0.46; t2=0.5: 2.12, 3.05, 2.38. Non-extensive. But finite-L sign reversals are NOT predicted by the leading exponent (t2=0.3, L=40, eps=-3.493: exact -0.4895 vs leading +0.4791).
A10 IMPRECISE: "thermally averaged ratio is the Fermi-window average of xi(eps)^-1": actually ln of a transmission-weighted average of exp. Example L=40,t2=0.5,kappa=0.2: actual 8.106; T_LR-weighted average of the log ratio 6.826; of the root exponent 6.648; Fermi-weight-only average 5.198.
A11 IMPRECISE: winding. E(pi/2)=-2 t2 is real, so w is undefined (not 0) at eps=-2 t2; w(-i kappa)=-1 only for t2<J (0 for t2>J).
A12 VERIFIED: same-chirality fermionic gain gives kappa->kappa(1+g).
A13 VERIFIED: mirror identity P h_eff P = h_eff^T, G_11=G_LL; sum rules hold only when all reservoir channels are included; sigma_f equality to 1e-16 (L=8 example).
A14 VERIFIED: sign of Eq.(15); but the symbol Edot_alpha is used with the opposite direction in Eq.(6) (absorbed by lead) and Eq.(15) (from reservoir into system).
A15 ERROR: Sec. VIII coupling g(c_j^dag a + e^{i phi} c_{j+1}^dag a + h.c.) eliminates to L_eff ∝ c_j + e^{-i phi} c_{j+1}, the conjugate phase of Eq.(1); at phi=pi/2 the diode direction is reversed. Factor 4g^2/Gamma_a correct. Exact Schur complement check: discrepancy 2.6e-4 with the manuscript's phase at g=0.01, Gamma_a=1, phi=0.7, zero with the conjugate.
A16 ERROR: circuit parameters (J=Gamma_a=10 MHz, g=2.2 MHz) do not justify the constant-kappa model: kappa(eps)/kappa(0)=1/(1+(2eps/Gamma_a)^2) = 1/5 at eps=J. Explicit-auxiliary chain L=15 (phase corrected): ln(T_RL/T_LR)=5.438 at eps=0, 1.083 at eps=J, 0.319 at eps=2J, versus Markov 5.438.
Top three: A16, A15, A8.
