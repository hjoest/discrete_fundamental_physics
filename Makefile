LATEX_IMAGE ?= pandoc/latex:latest
UID_GID     := $(shell id -u):$(shell id -g)

.PHONY: all pdf paper paper2 paper3 anonymous simulation test scrutiny analysis clean

all: pdf

pdf: docs/paper.pdf docs/paper2.pdf docs/paper3.pdf

paper: docs/paper.pdf

paper2: docs/paper2.pdf

paper3: docs/paper3.pdf

anonymous: docs/paper2_anonymous.pdf docs/paper3_anonymous.pdf

docs/paper.pdf: docs/paper.tex
	docker run --rm --user "$(UID_GID)" -e HOME=/tmp --entrypoint sh \
		-v "$(CURDIR)/docs:/work" -w /tmp $(LATEX_IMAGE) \
		-c "cp /work/paper.tex /tmp/ && pdflatex -interaction=nonstopmode paper.tex && pdflatex -interaction=nonstopmode paper.tex && pdflatex -interaction=nonstopmode paper.tex && cat /tmp/paper.pdf > /work/paper.pdf"
	rm -f docs/paper.aux docs/paper.log docs/paper.out docs/paper.toc

docs/paper2.pdf: docs/paper2.tex
	docker run --rm --user "$(UID_GID)" -e HOME=/tmp --entrypoint sh \
		-v "$(CURDIR)/docs:/work" -w /tmp $(LATEX_IMAGE) \
		-c "cp /work/paper2.tex /tmp/ && pdflatex -interaction=nonstopmode paper2.tex && pdflatex -interaction=nonstopmode paper2.tex && pdflatex -interaction=nonstopmode paper2.tex && cat /tmp/paper2.pdf > /work/paper2.pdf"
	rm -f docs/paper2.aux docs/paper2.log docs/paper.out docs/paper2.toc

docs/paper3.pdf: docs/paper3.tex
	docker run --rm --user "$(UID_GID)" -e HOME=/tmp --entrypoint sh \
		-v "$(CURDIR)/docs:/work" -w /tmp $(LATEX_IMAGE) \
		-c "cp /work/paper3.tex /tmp/ && pdflatex -interaction=nonstopmode paper3.tex && pdflatex -interaction=nonstopmode paper3.tex && pdflatex -interaction=nonstopmode paper3.tex && cat /tmp/paper3.pdf > /work/paper3.pdf"
	rm -f docs/paper3.aux docs/paper3.log docs/paper3.out docs/paper3.toc

docs/paper2_anonymous.pdf: docs/paper2.tex
	docker run --rm --user "$(UID_GID)" -e HOME=/tmp --entrypoint sh \
		-v "$(CURDIR)/docs:/work" -w /tmp $(LATEX_IMAGE) \
		-c "sed -e 's/pdfauthor={Holger Joest}/pdfauthor={}/g' \
		        -e 's|\\\\author{Holger Joest\\\\thanks{.*}|\\\\author{\\\\textit{Author name(s) and repository URL anonymized for double-anonymous peer review}}|g' \
		        -e 's|available at \\\\url{https://github.com/hjoest/discrete_fundamental_physics}\\.|provided as supplementary material (repository URL anonymized for double-anonymous peer review).|g' \
		        -e 's|H.~Joest, \`\`Discrete|[Author(s) anonymized], \`\`Discrete|g' \
		        -e 's|, \\\\href{https://doi.org/10.5281/zenodo.[0-9]*}{doi:10.5281/zenodo.[0-9]*}| [DOI anonymized for double-anonymous review]|g' \
		        /work/paper2.tex > /tmp/paper2_anonymous.tex && \
		    pdflatex -interaction=nonstopmode paper2_anonymous.tex && \
		    pdflatex -interaction=nonstopmode paper2_anonymous.tex && \
		    pdflatex -interaction=nonstopmode paper2_anonymous.tex && \
		    cat /tmp/paper2_anonymous.tex > /work/paper2_anonymous.tex && \
		    cat /tmp/paper2_anonymous.pdf > /work/paper2_anonymous.pdf"

docs/paper3_anonymous.pdf: docs/paper3.tex
	docker run --rm --user "$(UID_GID)" -e HOME=/tmp --entrypoint sh \
		-v "$(CURDIR)/docs:/work" -w /tmp $(LATEX_IMAGE) \
		-c "sed -e 's/pdfauthor={Holger Joest}/pdfauthor={}/g' \
		        -e 's|\\\\author{Holger Joest\\\\thanks{.*}|\\\\author{\\\\textit{Author name(s) and repository URL anonymized for double-anonymous peer review}}|g' \
		        -e 's|available at \\\\url{https://github.com/hjoest/discrete_fundamental_physics}\\.|provided as supplementary material (repository URL anonymized for double-anonymous peer review).|g' \
		        -e 's|H.~Joest, \`\`Discrete|[Author(s) anonymized], \`\`Discrete|g' \
		        -e 's|H.~Joest, \`\`Numerical|[Author(s) anonymized], \`\`Numerical|g' \
		        -e 's|, \\\\href{https://doi.org/10.5281/zenodo.[0-9]*}{doi:10.5281/zenodo.[0-9]*}| [DOI anonymized for double-anonymous review]|g' \
		        /work/paper3.tex > /tmp/paper3_anonymous.tex && \
		    pdflatex -interaction=nonstopmode paper3_anonymous.tex && \
		    pdflatex -interaction=nonstopmode paper3_anonymous.tex && \
		    pdflatex -interaction=nonstopmode paper3_anonymous.tex && \
		    cat /tmp/paper3_anonymous.tex > /work/paper3_anonymous.tex && \
		    cat /tmp/paper3_anonymous.pdf > /work/paper3_anonymous.pdf"

simulation:
	python3 run_simulation.py

test:
	python3 run_falsification.py

scrutiny:
	python3 run_scrutiny.py

analysis:
	python3 run_analysis.py

clean:
	rm -f docs/*.aux docs/*.log docs/*.out docs/*.toc docs/*_anonymous.pdf docs/*_anonymous.tex
