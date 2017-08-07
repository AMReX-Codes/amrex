SRC_FILES = theory.tex \
            AMR/AMR.tex \
			PML/PML.tex

all: $(SRC_FILES)
	rm -f allbibs.bib;cat *.bib > allbibs.bib
	rm -f theory.rst
	pandoc theory.tex --mathjax -S --wrap=preserve --bibliography allbibs.bib -o theorycore.rst
	cat theory_header.rst theorycore.rst > theory.rst
