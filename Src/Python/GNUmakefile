#
# Makefile for C++ versions of PyBoxLib
#

OUT ?= .

COMP  ?= gcc
FCOMP ?= gfortran
CXX   ?= mpic++
CC    ?= mpicc

all:
	make -f GNUmakefile_CXX.mak DIM=1 COMP=$(COMP) FCOMP=$(FCOMP) CC=$(CC) CXX=$(CXX) OUT=$(OUT)
	make -f GNUmakefile_CXX.mak DIM=2 COMP=$(COMP) FCOMP=$(FCOMP) CC=$(CC) CXX=$(CXX) OUT=$(OUT)
	make -f GNUmakefile_CXX.mak DIM=3 COMP=$(COMP) FCOMP=$(FCOMP) CC=$(CC) CXX=$(CXX) OUT=$(OUT)

wrappers:
	make -f GNUmakefile_CXX.mak DIM=1 wrapper
	make -f GNUmakefile_CXX.mak DIM=2 wrapper
	make -f GNUmakefile_CXX.mak DIM=3 wrapper

clean:
	rm -rf boxlib/*.so boxlib/*.pyc *.pyc

realclean:
	rm -rf build/ boxlib/*.so boxlib/*.pyc *.pyc d f o

.PHONY: clean
