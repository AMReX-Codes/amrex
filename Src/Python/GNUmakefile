#
# Makefile for C++ versions of PyBoxLib
#

OUT ?= .


all:
	make -f GNUmakefile_CXX.mak DIM=1 OUT=$(OUT)
	make -f GNUmakefile_CXX.mak DIM=2 OUT=$(OUT)
	make -f GNUmakefile_CXX.mak DIM=3 OUT=$(OUT)

wrappers:
	make -f GNUmakefile_CXX.mak DIM=1 wrapper
	make -f GNUmakefile_CXX.mak DIM=2 wrapper
	make -f GNUmakefile_CXX.mak DIM=3 wrapper

clean:
	rm -rf boxlib/*.so boxlib/*.pyc *.pyc

realclean:
	rm -rf build/ boxlib/*.so boxlib/*.pyc *.pyc d f o


.PHONY: clean
