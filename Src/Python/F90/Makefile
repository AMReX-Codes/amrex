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

PYINCLUDE := -I$(shell python -c 'import distutils.sysconfig; print distutils.sysconfig.get_python_inc()')
PYINCLUDE += -I$(shell python -c 'import numpy; print numpy.get_include()')
PYLDFLAGS  =   $(shell python-config --ldflags)

#all: $(OUT)/libboxlib.so
all: $(OUT)/fcboxlib.so

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
include $(BOXLIB_HOME)/Src/Python/GPackage.mak

VPATH_LOCATIONS += src
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak

$(OUT)/libboxlib.so: $(objects)
	mkdir -p $(OUT)
	$(F90) $(F90FLAGS) -shared -o $@ $^ $(PYLDFLAGS)

$(OUT)/fcboxlib.so: src/fboxlib_c.c $(objects)
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -shared -o $@ $(PYINCLUDE) $^ $(PYLDFLAGS)
