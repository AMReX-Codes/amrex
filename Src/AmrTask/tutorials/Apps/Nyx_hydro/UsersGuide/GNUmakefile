EPStoPDF = epstopdf

ALL: NyxUserGuide.ps NyxUserGuide.pdf 

eps_source = $(wildcard *.eps) 

pdf_source = $(eps_source:.eps=.pdf)

section_sources = Introduction/NyxIntroduction.tex \
                  GettingStarted/NyxGettingStarted.tex \
	          Inputs/NyxInputs.tex \
	          UnitsAndConstants/NyxUnits.tex \
	          ComovingHydro/Equations.tex \
	          Gravity/NyxGravity.tex \
                  Particles/Particles.tex \
                  Visualization/NyxVisualization.tex 

NyxUserGuide.dvi: NyxUserGuide.tex $(eps_source) $(section_sources) $(eps_source)
	latex NyxUserGuide.tex < /dev/null
	bibtex NyxUserGuide
	latex NyxUserGuide.tex < /dev/null
	latex NyxUserGuide.tex < /dev/null
	latex NyxUserGuide.tex < /dev/null

NyxUserGuide.pdf: NyxUserGuide.tex $(pdf_source) $(section_sources) $(pdf_source)
	pdflatex NyxUserGuide.tex < /dev/null
	bibtex NyxUserGuide
	pdflatex NyxUserGuide.tex < /dev/null
	pdflatex NyxUserGuide.tex < /dev/null
	pdflatex NyxUserGuide.tex < /dev/null

pdf:	NyxUserGuide.pdf 

%.ps: %.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) NyxUserGuide.dvi
	$(RM) NyxUserGuide.blg NyxUserGuide.log NyxUserGuide.toc
	$(RM) NyxUserGuide.aux NyxUserGuide.ps NyxUserGuide.bbl NyxUserGuide.pdf
	$(RM) *~

.PHONY: clean
