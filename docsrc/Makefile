# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -aE -b $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

github:
	@echo "---- Compiling latex file of diagram ------"
	@pdflatex -output-directory=./texsource/diagram/ ./texsource/diagram/diagram_models.tex > ./texsource/diagram/pdflatex.log
	@echo "---- Converting pdf to png ------"
	@convert -density 400 ./texsource/diagram/diagram_models.pdf ./_static/diagram_models.png
	@echo "---- Building html files ------"
	@make html
	@echo "---- Copying html files to ../docs/ ------"
	@cp -a _build/. ../docs/
