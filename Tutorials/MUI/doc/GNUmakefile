
EPStoPDF = epstopdf

ALL: MUIcouplingNotes.ps MUIcouplingNotes.pdf

eps_source = $(wildcard *.eps)

pdf_source = $(eps_source:.eps=.pdf)

MUIcouplingNotes.dvi: MUIcouplingNotes.tex $(eps_source)
	latex MUIcouplingNotes.tex < /dev/null
	latex MUIcouplingNotes.tex < /dev/null
	latex MUIcouplingNotes.tex < /dev/null
	latex MUIcouplingNotes.tex < /dev/null

MUIcouplingNotes.pdf: MUIcouplingNotes.tex $(pdf_source)
	pdflatex MUIcouplingNotes.tex < /dev/null
	pdflatex MUIcouplingNotes.tex < /dev/null
	pdflatex MUIcouplingNotes.tex < /dev/null
	pdflatex MUIcouplingNotes.tex < /dev/null

pdf:	MUIcouplingNotes.pdf

MUIcouplingNotes.ps: MUIcouplingNotes.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) MUIcouplingNotes.dvi 
	$(RM) MUIcouplingNotes.blg MUIcouplingNotes.log MUIcouplingNotes.pdf
	$(RM) MUIcouplingNotes.aux MUIcouplingNotes.ps MUIcouplingNotes.bbl
	$(RM) MUIcouplingNotes.out *~

.PHONY: clean
