all: main.pdf

main.pdf: main.tex abstract.tex introduction.tex methods.tex results.tex discussion.tex refs.bib evolution.sty
	pdflatex main.tex
	bibtex main
	pdflatex main.tex
	pdflatex main.tex

clean:
	rm -f *.bbl *.blg *.log *.aux *.toc *.fff *.ent *.ttt *.synctex.gz
