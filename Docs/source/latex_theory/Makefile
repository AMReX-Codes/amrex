SRC_FILES = theory.tex \
            AMR/AMR.tex \
	    PML/PML.tex

all: $(SRC_FILES) clean
	pandoc theory.tex --mathjax -S --wrap=preserve --bibliography allbibs.bib -o warpx_theory.rst
	cp warpx_theory.rst ../theory
	cd ../../../../picsar/Doxygen/pages/latex_theory/; pandoc theory.tex --mathjax -S --wrap=preserve --bibliography allbibs.bib -o picsar_theory.rst
	cp ../../../../picsar/Doxygen/pages/latex_theory/picsar_theory.rst ../theory

clean: 
	rm -f warpx_theory.rst
	rm -f ../../../../picsar/Doxygen/pages/latex_theory/picsar_theory.rst
