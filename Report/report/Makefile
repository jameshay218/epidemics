all: litreview

litreview: thesis.tex
	make clean
	pdflatex thesis 
	pdflatex thesis
	pdflatex thesis
	pdflatex thesis
	make clean
	evince thesis.pdf


clean:
	rm -f *.aux *.bbl *.log *.xml *.toc *.blg *.dvi *.out *~* \#*
