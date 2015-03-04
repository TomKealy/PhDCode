.phony: all dvi pdf  ps view clean

BASE = by-example

all: pdf


pdf: $(BASE).pdf


$(BASE).pdf: $(BASE).tex
	pdflatex $(BASE)

view: $(BASE).pdf
	open $(BASE).pdf


clean:
	rm -f $(BASE).log $(BASE).aux $(BASE).dvi $(BASE).out $(BASE).ps