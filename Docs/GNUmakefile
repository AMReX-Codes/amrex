EPStoPDF = epstopdf

ALL: paper.ps paper.pdf 

eps_source = $(wildcard BoxLib/*.eps) 

pdf_source = $(eps_source:.eps=.pdf)

section_sources = BoxLib/boxlib.tex 

paper.dvi: paper.tex $(eps_source) $(section_sources) $(eps_source) 
	latex paper.tex < /dev/null
	bibtex paper
	latex paper.tex < /dev/null
	latex paper.tex < /dev/null
	latex paper.tex < /dev/null

paper.pdf: paper.tex $(pdf_source) $(section_sources) $(pdf_source)
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
