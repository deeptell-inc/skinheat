# Submission bundle — "Exact exponential heat rectification from the Liouvillian skin effect"

Target: Physical Review Letters. Built 2026-10-29 from `skinheat/paper/` (sources) and `skinheat/results/` (data). All PDFs compile inside this directory with pdflatex + bibtex; nothing here has been submitted or uploaded.

| File | Purpose |
|---|---|
| `main.tex`, `main.pdf`, `main.bbl` | Letter (REVTeX 4.2 `prl` reprint, 4 figures, 6 pages incl. references) |
| `sm.tex`, `sm.pdf`, `sm.bbl` | Supplemental Material S1–S7 (7 pages, 4 tables) |
| `refs.bib` | Bibliography source, all entries verified |
| `figures/fig1–4.pdf` | Figures (regenerate with `../make_figures.py`) |
| `cover_letter.tex`, `cover_letter.pdf` | Cover letter: claim, why PRL, 5 suggested referees, length statement, prior-submission statement |
| `DECLARATIONS.md` | Competing interests, funding, contributions, data availability |
| `submission-requirements.md` | PRL requirement × deliverable checklist and the items the authors still fill in |
| `arxiv/` | Single-file arXiv version (`main_arxiv.tex`, Letter + SM, 10 pages), `main_arxiv.bbl`, figures, `abstract.txt`, `SUBMISSION_NOTES.md` (form fields), `skinheat_arxiv.tar.gz` (300 KB, compiles with two pdflatex passes) |
| `zenodo/` | `README.md`, `CITATION.cff`, `submission.json` (zenodo-deposit config), `make_deposit.sh` (rebuilds `skinheat-v0.1.0.tar.gz`, 4.5 MB, and runs an offline dry-run only) |

Order of operations at submission time: (1) `zenodo-deposit --production` (draft) → publish from the Zenodo UI → copy the DOI into `main.tex` (Data availability), `cover_letter.tex`, `DECLARATIONS.md`, `arxiv/SUBMISSION_NOTES.md`; (2) submit `main.tex` + figures + `sm.pdf` + `cover_letter.pdf` + the statement in `data-availability-statement.md` to PRL; (3) after publication, post `arxiv/skinheat_arxiv.tar.gz` to arXiv with the journal reference in the Comments line.
