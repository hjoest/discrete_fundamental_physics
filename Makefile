LATEX_IMAGE ?= pandoc/latex:latest
UID_GID     := $(shell id -u):$(shell id -g)

.PHONY: all pdf paper paper2 simulation test scrutiny analysis clean

all: pdf

pdf: docs/paper.pdf docs/paper2.pdf

paper: docs/paper.pdf

paper2: docs/paper2.pdf

docs/paper.pdf: docs/paper.tex
	docker run --rm --user "$(UID_GID)" -e HOME=/tmp --entrypoint sh \
		-v "$(CURDIR)/docs:/work" -w /tmp $(LATEX_IMAGE) \
		-c "cp /work/paper.tex /tmp/ && pdflatex -interaction=nonstopmode paper.tex && pdflatex -interaction=nonstopmode paper.tex && pdflatex -interaction=nonstopmode paper.tex && cat /tmp/paper.pdf > /work/paper.pdf"
	rm -f docs/paper.aux docs/paper.log docs/paper.out docs/paper.toc

docs/paper2.pdf: docs/paper2.tex
	docker run --rm --user "$(UID_GID)" -e HOME=/tmp --entrypoint sh \
		-v "$(CURDIR)/docs:/work" -w /tmp $(LATEX_IMAGE) \
		-c "cp /work/paper2.tex /tmp/ && pdflatex -interaction=nonstopmode paper2.tex && pdflatex -interaction=nonstopmode paper2.tex && pdflatex -interaction=nonstopmode paper2.tex && cat /tmp/paper2.pdf > /work/paper2.pdf"
	rm -f docs/paper2.aux docs/paper2.log docs/paper2.out docs/paper2.toc

simulation:
	python3 run_simulation.py

test:
	python3 run_falsification.py

scrutiny:
	python3 run_scrutiny.py

analysis:
	python3 run_analysis.py

clean:
	rm -f docs/*.aux docs/*.log docs/*.out docs/*.toc
