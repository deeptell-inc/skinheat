# PRB submission package — "Exact exponential heat rectification from the Liouvillian skin effect"

Regular article for Physical Review B, built 2026-11-11 from the PRL Letter and its Supplemental Material (both integrated into one main text; no Supplemental Material). Nothing here has been submitted.

| File | Purpose |
|---|---|
| `main.tex`, `main.pdf`, `main.bbl` | Article: REVTeX 4.2 `prb` reprint, 10 pages, 9 sections, 4 figures, 4 tables, 50 references; abstract 227 words; main text about 4603 words |
| `refs.bib`, `figures/fig1–4.pdf` | Bibliography (all entries verified) and figures (regenerate with `../make_figures.py`) |
| `cover_letter.tex`, `cover_letter.pdf` | Summary, why PRB, **disclosure of the PRL submission LX20531**, five suggested referees, data statement |

## Data Availability Statement (submission form)

> The data that support the findings of this article, including the matrix-product-state time series behind Fig. 3, are publicly available on Zenodo [46]. The code used to generate and analyse the data and the scripts that regenerate every figure from the stored data files are available from Refs. [46,47].

## Checklist
- [x] Compiles with pdflatex + bibtex, 0 errors, 0 undefined references, 0 overfull boxes
- [x] Hygiene grep (`corrected|earlier (revision|draft)|previously reported|referee|phase-1|Letter|Supplemental`) on `main.tex`: no hits
- [x] Previous APS submission disclosed in the cover letter (LX20531)
- [ ] Zenodo deposit 23801927 published, with a tarball rebuilt by `../submission-bundle/zenodo/make_deposit.sh` (the draft uploaded on 2026-10-30 predates the figure fixes and contains the LaTeX sources)
- [ ] arXiv: after publication (submission limit)
- [ ] Confirm the current affiliations of the suggested referees (Yoshida, Yokomizo in particular) in the submission form
