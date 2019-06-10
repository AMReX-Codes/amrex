SRC_FILES = theory.tex \
            AMR/AMR.tex \
	    PML/PML.tex

all: $(SRC_FILES) clean
	pandoc intro.tex --mathjax --wrap=preserve --bibliography allbibs.bib -o intro.rst
	pandoc AMR/AMR.tex --mathjax --wrap=preserve --bibliography allbibs.bib -o amr.rst
	pandoc Boosted_frame/Boosted_frame.tex --mathjax --wrap=preserve --bibliography allbibs.bib -o boosted_frame.rst
	pandoc input_output/input_output.tex --mathjax --wrap=preserve --bibliography allbibs.bib -o input_output.rst
	pandoc PML/PML.tex --mathjax --wrap=preserve --bibliography allbibs.bib -o PML.rst
	mv *.rst ../theory
	cd ../../../../picsar/Doxygen/pages/latex_theory/; pandoc theory.tex --mathjax --wrap=preserve --bibliography allbibs.bib -o picsar_theory.rst
	mv ../../../../picsar/Doxygen/pages/latex_theory/picsar_theory.rst ../theory
	cp ../../../../picsar/Doxygen/images/PIC.png ../theory
	cp ../../../../picsar/Doxygen/images/Yee_grid.png ../theory

clean:
	rm -f ../theory/intro.rst
	rm -f ../theory/warpx_theory.rst
	rm -f ../theory/picsar_theory.rst
	rm -f ../theory/PIC.png
	rm -f ../theory/Yee_grid.png
