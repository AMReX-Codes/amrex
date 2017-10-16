SRC_FILES = theory.tex \
            AMR/AMR.tex \
	    PML/PML.tex

all: $(SRC_FILES) clean
	pandoc intro.tex --mathjax -S --wrap=preserve --bibliography allbibs.bib -o intro.rst
	mv intro.rst ../theory
	pandoc theory.tex --mathjax -S --wrap=preserve --bibliography allbibs.bib -o warpx_theory.rst
	mv warpx_theory.rst ../theory
	cd ../../../../picsar/Doxygen/pages/latex_theory/; pandoc theory.tex --mathjax -S --wrap=preserve --bibliography allbibs.bib -o picsar_theory.rst
	mv ../../../../picsar/Doxygen/pages/latex_theory/picsar_theory.rst ../theory
	cp ../../../../picsar/Doxygen/images/PIC.png ../theory
	cp ../../../../picsar/Doxygen/images/Yee_grid.png ../theory

clean:
	rm -f ../theory/intro.rst
	rm -f ../theory/warpx_theory.rst
	rm -f ../theory/picsar_theory.rst
	rm -f ../theory/PIC.png
	rm -f ../theory/Yee_grid.png
