# AI use disclosure for the APS submission form

APS policy (journals.aps.org/authors/appropriate-use-ai-tools): substantive AI use must be declared; the statement is printed in the paper under "AI Use Disclosure"; authors are fully responsible for accuracy; undisclosed substantive use is treated as misconduct. The use in this work is substantive in every category the policy lists (code, simulations, analysis, derivations, figures, drafting).

Confirmed by the authors on 2026-11-26: the Gaussian solver and the formulation of the exact law were also produced by Claude Fable 5.1. The two remaining items in [square brackets] are statements about what the authors did. Only the authors can confirm them. Edit or delete each one so that it is literally true, and do not submit the text with brackets left in.

---

AI tool name and version.
Claude Code (Anthropic), an agentic coding assistant, running the models Claude Fable 5.1 (October to early November 2026) and Claude Fable 5.2 (from 11 November 2026).

How the AI assisted.
(1) Software. The AI wrote the skinheat package: the Gaussian (Lyapunov) solver and transmission routines, the superfermion matrix-product-state Lindblad solver, the finite-temperature bath and disorder routines, the reproduction scripts and the test suite.
(2) Simulations and data analysis. The AI launched and monitored the numerical runs, carried out the convergence studies of the matrix-product-state calculations (bond dimension and time step) and analysed the results.
(3) Derivations. The AI worked out the formulation of the exact law through the tridiagonal Green-function identity, the non-Bloch form of the rectification exponent, the sum rule for bond disorder and the mirror identity for the entropy production.
(4) Figures. All figures were produced by scripts written by the AI from the archived data files.
(5) Literature. The AI searched for related work and checked the bibliographic data of every reference against Crossref or arXiv.
(6) Writing. The AI produced the first draft of the full text and restructured and edited it in later versions.

How the authors directed and verified the AI output.
The authors defined the research question and set each task, decided which extensions to study and which journal to target, and reviewed the results and the text at each stage, requesting corrections where needed. The numerical methods are validated by automated tests that are archived with the code: the exact law is reproduced to machine precision, the Gaussian solver agrees with the full many-body Lindbladian on small systems, and the matrix-product-state solver agrees with the exact solution in the non-interacting limit, which is the calibration reported in the paper. Every figure is regenerated from the archived data by a script. [The authors independently re-derived the exact law, Eq. (11), the non-Bloch form, Eq. (13), the mirror identity, Eq. (17), and the sum rule, Eq. (19).] [The authors read the complete manuscript, checked every numerical value quoted in the text against the archived result files, and verified the cited references.] The authors take full responsibility for the content of the article.
