#
# Makefile for C++ versions of PyBoxLib
#

all:
	make -f GNUmakefile_CXX.mak DIM=1
	make -f GNUmakefile_CXX.mak DIM=2
	make -f GNUmakefile_CXX.mak DIM=3

wrappers:
	make -f GNUmakefile_CXX.mak DIM=1 wrapper
	make -f GNUmakefile_CXX.mak DIM=2 wrapper
	make -f GNUmakefile_CXX.mak DIM=3 wrapper

clean:
	rm -rf boxlib/*.so boxlib/*.pyc *.pyc

realclean:
	rm -rf boxlib/*.so boxlib/*.pyc *.pyc d f o src/* boxlib/bl*.py


.PHONY: clean
