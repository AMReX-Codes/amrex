EPStoPDF = epstopdf

ALL: paper.ps paper.pdf 

eps_source = $(wildcard *.eps)

pdf_source = $(eps_source:.eps=.pdf)

paper.dvi: paper.tex architecture.tex $(eps_source) ../maestro_doc.bib
	latex paper.tex < /dev/null
	bibtex paper
	latex paper.tex < /dev/null
	latex paper.tex < /dev/null
	latex paper.tex < /dev/null

paper.pdf: paper.tex architecture.tex $(pdf_source) ../maestro_doc.bib
	pdflatex paper.tex < /dev/null
	bibtex paper
	pdflatex paper.tex < /dev/null
	pdflatex paper.tex < /dev/null
	pdflatex paper.tex < /dev/null

pdf:	paper.pdf 

%.ps: %.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) paper.dvi 
	$(RM) paper.blg paper.log
	$(RM) paper.aux paper.ps paper.bbl
	$(RM) *~

.PHONY: clean
