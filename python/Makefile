.PHONY: lib figures notebook clean

PYTHON ?= $(shell if [ -x /opt/homebrew/bin/python3 ]; then printf /opt/homebrew/bin/python3; else command -v python3; fi)
JUPYTER ?= $(shell if [ -x /opt/homebrew/bin/jupyter ]; then printf /opt/homebrew/bin/jupyter; else command -v jupyter; fi)
MPLCONFIGDIR ?= $(CURDIR)/.matplotlib

lib:
	$(PYTHON) -m cspack.build --quiet

figures: lib
	mkdir -p "$(MPLCONFIGDIR)"
	MPLCONFIGDIR="$(MPLCONFIGDIR)" $(PYTHON) notebooks/render_figures.py --write-notebook

notebook: figures
	mkdir -p "$(MPLCONFIGDIR)"
	MPLCONFIGDIR="$(MPLCONFIGDIR)" $(JUPYTER) lab notebooks/CSpack_examples.ipynb

clean:
	rm -rf build dist *.egg-info cspack/lib .matplotlib
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	find . -type d -name .ipynb_checkpoints -prune -exec rm -rf {} +
