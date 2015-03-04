.phony: all dvi pdf  ps view clean

BASE = matgraph

all: dvi pdf

dvi: $(BASE).dvi

pdf: $(BASE).pdf

$(BASE).dvi: $(BASE).tex
	latex $(BASE)
	latex $(BASE)

$(BASE).pdf: $(BASE).tex
	pdflatex $(BASE)

view: $(BASE).dvi
	xdvik $(BASE) &

ps: $(BASE).dvi
	dvips $(BASE) -t letter -o $(BASE).ps

clean:
	rm -f $(BASE).ps $(BASE).log $(BASE).aux $(BASE).dvi $(BASE).out