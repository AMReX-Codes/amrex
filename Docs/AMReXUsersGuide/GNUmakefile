EPStoPDF = epstopdf

ALL: AMReXUsersGuide.pdf 

eps_source = $(wildcard *.eps) $(wildcard */*.eps)

pdf_source = $(eps_source:.eps=.pdf)

section_sources = Preface/Preface.tex \
                  Introduction/Introduction.tex \
                  GettingStarted/GettingStarted.tex \
                  BuildingAMReX/BuildingAMReX.tex \
                  Basics/Basics.tex \
                  Boundary/Boundary.tex \
                  AmrCore/AmrCore.tex \
                  Particle/Particle.tex \
                  EB/EB.tex \
                  AmrLevel/AmrLevel.tex \
                  Fortran/Fortran.tex \
                  Visualization/Visualization.tex \
                  Profiling/Profiling.tex \
                  Debugging/Debugging.tex \
                  CVODE/CVODE.tex

AMReXUsersGuide.pdf: AMReXUsersGuide.tex $(pdf_source) $(section_sources) $(pdf_source)
	pdflatex AMReXUsersGuide.tex < /dev/null
#	bibtex AMReXUsersGuide
	pdflatex AMReXUsersGuide.tex < /dev/null
	makeindex AMReXUsersGuide < /dev/null
	pdflatex AMReXUsersGuide.tex < /dev/null
	pdflatex AMReXUsersGuide.tex < /dev/null

pdf:	AMReXUsersGuide.pdf 

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source)
	$(RM) AMReXUsersGuide.blg AMReXUsersGuide.log AMReXUsersGuide.toc
	$(RM) AMReXUsersGuide.aux AMReXUsersGuide.bbl
	$(RM) AMReXUsersGuide.lot AMReXUsersGuide.lof AMReXUsersGuide.out
	$(RM) AMReXUsersGuide.idx AMReXUsersGuide.ind AMReXUsersGuide.ilg
	$(RM) *~

realclean: clean
	$(RM) AMReXUsersGuide.pdf

.PHONY: clean realclean
