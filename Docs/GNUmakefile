EPStoPDF = epstopdf

ALL: UsersGuide.ps UsersGuide.pdf 

eps_source = $(wildcard Preface/*.eps) \
	     $(wildcard Introduction/*.eps) \
	     $(wildcard GettingStarted/*.eps) \
	     $(wildcard F_AdvancedTopics/*.eps)

pdf_source = $(eps_source:.eps=.pdf)

section_sources = Preface/Preface.tex \
                  Introduction/Introduction.tex \
                  GettingStarted/GettingStarted.tex \
                  F_AdvancedTopics/F_AdvancedTopics.tex \
	          C_AdvancedTopics/C_AdvancedTopics.tex \
                  Regression/test_suite.tex

UsersGuide.dvi: UsersGuide.tex $(eps_source) $(section_sources) $(eps_source) 
	latex UsersGuide.tex < /dev/null
	latex UsersGuide.tex < /dev/null

UsersGuide.pdf: UsersGuide.tex $(pdf_source) $(section_sources) $(pdf_source)
	pdflatex UsersGuide.tex < /dev/null
	pdflatex UsersGuide.tex < /dev/null

pdf: UsersGuide.pdf 

%.ps: %.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) UsersGuide.dvi
	$(RM) UsersGuide.blg UsersGuide.log UsersGuide.lof
	$(RM) UsersGuide.aux UsersGuide.ps UsersGuide.bbl
	$(RM) UsersGuide.toc UsersGuide.pdf *~

.PHONY: clean
