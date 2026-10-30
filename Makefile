PY ?= python
TEX = pdflatex -interaction=nonstopmode

.PHONY: install test test-all figures paper bundle zenodo clean

install:            ## editable install with test deps (add [mps] for TeNPy)
	$(PY) -m pip install -e ".[dev]"

test:               ## fast tests (~4 min)
	$(PY) -m pytest tests -q -m "not slow"

test-all:           ## including the slow MPS-vs-exact comparison
	$(PY) -m pytest tests -q

figures:            ## Figs. 1-4 from results/*.json
	$(PY) paper/make_figures.py

paper: figures      ## Letter + Supplemental Material
	cd paper && $(TEX) main >/dev/null && bibtex main >/dev/null && $(TEX) main >/dev/null && $(TEX) main >/dev/null
	cd paper && $(TEX) sm >/dev/null && bibtex sm >/dev/null && $(TEX) sm >/dev/null && $(TEX) sm >/dev/null
	@grep -c '^!' paper/main.log paper/sm.log || true

bundle:             ## arXiv tarball (single-file Letter + SM)
	cd paper/submission-bundle/arxiv && $(TEX) main_arxiv >/dev/null && bibtex main_arxiv >/dev/null && $(TEX) main_arxiv >/dev/null && $(TEX) main_arxiv >/dev/null \
	  && tar czf skinheat_arxiv.tar.gz main_arxiv.tex main_arxiv.bbl refs.bib figures abstract.txt && ls -la skinheat_arxiv.tar.gz

zenodo:             ## Zenodo tarball + offline dry-run (no upload)
	bash paper/submission-bundle/zenodo/make_deposit.sh

clean:
	rm -rf .pytest_cache skinheat.egg-info skinheat/__pycache__ tests/__pycache__
	cd paper && rm -f *.aux *.log *.out *.blg *.fls *.fdb_latexmk *.synctex.gz *Notes.bib *.compile.log
