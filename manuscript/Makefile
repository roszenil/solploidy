all: main.pdf

main.pdf: main.tex abstract.tex introduction.tex methods.tex results.tex discussion.tex refs.bib evolution.sty
	pdflatex main.tex
	bibtex main
	pdflatex main.tex
	pdflatex main.tex
	mv main.pdf temp.pdf
	pdftk temp.pdf ../finalfigs/allSSEmodels.pdf output main.pdf
	rm temp.pdf

clean:
	rm -f *.bbl *.blg *.log *.aux *.toc *.fff *.ent *.ttt *.out *.efsuppfigure *.efsupptable *.synctex.gz
