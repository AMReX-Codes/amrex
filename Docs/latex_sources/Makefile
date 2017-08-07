SRC_FILES = theory.tex \
            AMR/AMR.tex \
			PML/PML.tex

all: $(SRC_FILES)
	rm -f allbibs.bib;cat *.bib > allbibs.bib
	pandoc theory.tex --mathjax -S --wrap=preserve --bibliography allbibs.bib -o warpx_theory.rst
	pandoc ../../../picsar/Doxygen/latex_theory/theory.tex --mathjax -S --wrap=preserve --bibliography allbibs.bib -o picsar_theory.rst
	cp *.rst ../source/theory/
