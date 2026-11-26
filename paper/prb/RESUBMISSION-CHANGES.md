# Resubmission text for the APS form (previous manuscript LX20531)

Paste the two blocks below into the text box. Plain text, no title change requested.

---

RESPONSE TO THE REVIEWER COMMENTS

The previous version of this work (LX20531, submitted to Physical Review Letters as a Letter) was not sent for external review, so there are no reviewer comments to answer. The editors stated that they had not judged the validity of the manuscript and concluded that it was unlikely to meet the criteria of that journal for impact, innovation and interest. In response we have rewritten the work as a full-length regular article for Physical Review B. The derivations, the precision analysis of the numerics and the supporting data, which were in the Supplemental Material of the Letter, are now part of the main text, so that every statement can be checked without a separate file. The scientific results are the same as in LX20531. While preparing this version we rechecked every quoted number against the archived result files; items 8 and 16 below list the statements and numbers that were corrected. The title is unchanged.

LIST OF CHANGES (relative to LX20531)

1. Format. The Letter and its Supplemental Material are merged into one regular article with nine numbered sections, four figures, four tables and 50 references. There is no Supplemental Material.

2. Abstract. Rewritten for the regular-article format (228 words). It now states the origin of the law (an identity of tridiagonal Green functions), the sum rule for bond disorder and the validity for bosonic carriers.

3. Introduction. A paragraph on related work was added (diffusive skin effect in thermal metamaterials, nonreciprocal routing of thermal noise in cascaded networks, transport theory of non-Hermitian nonreciprocal systems, boundary-driven open quantum systems); these references were previously cited only in the Supplemental Material. A paragraph outlining the sections and stating the units and default parameters was added.

4. Section II (model and formalism). The Lyapunov equation for the Gaussian steady state, the two implementations of the leads and the definition of the heat-current responses as explicit equations were moved from the Supplemental Material into the main text.

5. Section III (exact law). The proof of the tridiagonal identity is now given in full in the main text. The consequences of the law are listed explicitly (independence of the lead distributions and of the size of the temperature difference, validity for bosons, independence of on-site potentials, validity for mesoscopic leads). The treatment of same-chirality gain was moved from the Supplemental Material.

6. Sections IV to VIII. The derivation of the non-Bloch form of the exponent and the periodic-boundary spectrum (Sec. IV), the superfermion matrix-product-state method, its calibration table, the resolvable regime and the table of interacting results (Sec. V), the construction of the finite-temperature bath, the dependence on its discretisation and the proof of the mirror identity for the entropy production (Sec. VI), the disorder data with the table of decay rates (Sec. VII) and the table of experimental parameters (Sec. VIII) were moved from the Supplemental Material into the main text.

7. Section V. One sentence was added noting that the free-fermion bias of the matrix-product-state result is an order of magnitude smaller for the infinite-temperature bath than for pure loss, with the two values (0.001 and 0.045) that were already listed in the tables of the previous version.

8. Correction in the Conclusion. The previous version stated that interactions, next-nearest-neighbour hopping and Anderson localisation all lower the throughput and never change the slope. This was inaccurate for next-nearest-neighbour hopping, which changes the exponent, as the results section of the previous version already showed. The Conclusion now states that interactions and Anderson localisation act on the prefactor and leave the slope unchanged, and that next-nearest-neighbour hopping replaces the constant exponent by its energy-resolved non-Bloch form.

9. Conclusion. The comparison with nonlinear thermal rectifiers no longer quotes a factor ("two orders of magnitude"); it now states only that their rectification ratios are independent of size.

10. The condition for a rectification above 0.99 is quoted as (L-1)/xi > 2.65 throughout; the previous version rounded it to 2.6 in one place.

11. A remark in the Supplemental Material about an alternative linear-response evaluation that was not used for any reported result has been removed.

12. Wording. Long sentences were split and several formulations were made more direct. Figure captions no longer use bold lead sentences. These edits do not change any claim or number.

13. Figures. The error bars of Fig. 3 now include the time fluctuation of the steady state (0.02 at kappa = 0.05J, 0.06 at kappa = 0.1J, 0.08 for the bath run; 4% for the slopes in panel c), and the caption states this. [DELETE THE NEXT SENTENCE IF THE OTHER FIGURES ARE THE SAME FILES AS IN LX20531.] The axis range of Fig. 4(c) was extended so that all data are visible, the range of Fig. 1(b) was set to 0-75 (see item 16), and the legend of Fig. 3(a) was moved off the data. The data are unchanged.

14. Cover letter. The list of suggested referees was revised for the readership of Physical Review B.

15. An "AI use disclosure" section was added before the Data availability statement, in accordance with the APS policy on AI tools.

16. Numerical corrections. A scripted comparison of every number quoted in the text with the archived result files (now part of the public repository) found the following errors in the previous version, which are corrected. None of them affects the exact law, the figures' data or the conclusions.
(a) Sec. III and the caption of Fig. 1: at kappa = 0.07J the forward throughput is 45% of its lossless value, not 49%.
(a2) Sec. III and Fig. 1(b): the agreement with the exact law was described as holding "to all digits" for all points. It holds to 2e-11 in the logarithm up to a ratio of e^48 and to 2e-6 at e^72; one parameter set beyond this (kappa = 0.3J, L = 160), whose reverse response is below double precision, is now omitted from the figure and this is stated in the text.
(b) Sec. II: the two lead implementations agree within 6% for K = 80; for K = 10 the throughput is 28-38% below the wide-band value (previously "within 5-15%"). The ratio is the same in all cases.
(c) Sec. V: the exact check of the vectorised Lindbladian holds to better than 1e-10 (previously quoted as 1e-15); the time fluctuation of the logarithmic ratio is at most 0.011 at kappa = 0.05J, 0.03 at kappa = 0.1J and 0.08 for the bath run (previously "less than 0.01" and 0.05); the statement on the error floor now distinguishes the absolute lead currents (about 7e-4 and 6e-4 at bond dimension 128 and 256) from the reverse response (3.8e-4 and 1.5e-6); the lower bounds on the rectification at kappa = 0.2J are 0.68-0.94 (previously 0.85-0.98).
(c2) Sec. V and the caption of Fig. 3: the text and caption now state which results use bond dimension 256 with time step 0.05/J (Table II, Figs. 3(a), 3(b) and the filled symbols of Fig. 3(c)) and which use bond dimension 128 with time step 0.1/J (open symbols of Fig. 3(c) and Fig. 3(d)); the previous version attributed the former setting to all runs.
(d) Sec. V: the statement on the renormalisation of the skin length now gives the measured deviation (2.5%) and the uncertainty (about 4%) separately; the interaction offset for the infinite-temperature bath is quoted with its uncertainty, -0.54 +/- 0.08, and described as consistent with the pure-loss value.
(e) Sec. VI: the two responses are independent of the bath state to 2e-13 (previously the ratio was said to agree to 1e-9, which does not hold at kappa = 0.2J where the reverse response is of order 4e-9); the mirror identity is broken by 3e-4 for disorder of width J (previously 3e-3); the caption of Fig. 4 now says that a particle sink costs about 500 times more entropy (previously "three orders of magnitude").
(f) Sec. VII: the sum rule for bond disorder is verified to 3e-14 (previously 1e-14); the median throughput at disorder width 4J is 1.6e-7 (previously misprinted as 1.6e-6); a throughput above 1e-3 requires a disorder width below 1.5J (previously 2J).
(g) Sec. VIII and Table IV: the circuit parameters are now mutually consistent (J = 10 MHz, g = 1.1-2.2 MHz, giving kappa = 0.5-2 MHz and kappa/J = 0.05-0.2); the chain length needed for a rectification above 0.99 at kappa/J = 0.05 is 54 (previously rounded to 60).
