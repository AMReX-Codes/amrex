#
# Makefile for PyBoxLib
#

BOXLIB_HOME ?= ../../..

COMP   := gfortran
NDEBUG := t
OMP    :=
MPI    := t

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

F90      := mpif90
FC       := mpif90
CC       := mpicc
F90FLAGS += -fPIC
CFLAGS   += -fPIC

PYINCLUDE := $(shell python -c 'import distutils.sysconfig; print distutils.sysconfig.get_python_inc()')
PYINCLUDE += $(shell python -c 'import numpy; print numpy.get_include()')

all: $(OUT)/fcboxlib.so

VPATH_LOCATIONS += src
INCLUDE_LOCATIONS += $(PYINCLUDE)
csources += fboxlib_c.c

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
include $(BOXLIB_HOME)/Src/Python/GPackage.mak

VPATH_LOCATIONS += src
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak

$(OUT)/fcboxlib.so: $(objects)
	mkdir -p $(OUT)
	$(F90) -shared -o $@ $^
