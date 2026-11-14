LATEX_IMAGE ?= pandoc/latex:latest
UID_GID     := $(shell id -u):$(shell id -g)

.PHONY: all pdf paper test analysis clean

all: pdf

pdf paper: docs/paper.pdf

docs/paper.pdf: docs/paper.tex
	docker run --rm --user "$(UID_GID)" -e HOME=/tmp --entrypoint sh \
		-v "$(CURDIR)/docs:/work" -w /tmp $(LATEX_IMAGE) \
		-c "cp /work/paper.tex /tmp/ && pdflatex -interaction=nonstopmode paper.tex && pdflatex -interaction=nonstopmode paper.tex && pdflatex -interaction=nonstopmode paper.tex && cat /tmp/paper.pdf > /work/paper.pdf"
	rm -f docs/paper.aux docs/paper.log docs/paper.out docs/paper.toc

test:
	python3 run_falsification.py

analysis:
	python3 run_analysis.py

clean:
	rm -f docs/*.aux docs/*.log docs/*.out docs/*.toc
