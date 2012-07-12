#
# Makefile for PyBoxLib
#

BOXLIB_HOME ?= $(HOME)/Development/BoxLib

COMP   := gfortran
NDEBUG := t
OMP    := t
MPI    := t

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

F90      := mpif90
FC       := mpif90
CC       := mpicc
F90FLAGS += -fPIC
CFLAGS   += -fPIC

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib
VPATH_LOCATIONS += src

f90sources += blobjects.f90
f90sources += fboxlib.f90
f90sources += boxlib_numpy_f.f90

PYFBOXLIB = pyboxlib/libpyfboxlib.so
PYCBOXLIB = pyboxlib/libpycboxlib.so
PYBOXLIBS = $(PYFBOXLIB) $(PYCBOXLIB)

all: $(PYBOXLIBS)

NPINCLUDE:=-I$(shell python -c 'import distutils.sysconfig; print distutils.sysconfig.get_python_inc()')
NPINCLUDE+=-I$(shell python -c 'import numpy; print numpy.get_include()')/numpy

$(PYFBOXLIB): $(objects)
	$(F90) $(F90FLAGS) -shared -o $@ $^

$(PYCBOXLIB): boxlib_numpy_c.c
	cd pyboxlib && $(CC) $(CFLAGS) $(NPINCLUDE) -shared -o libpycboxlib.so ../src/boxlib_numpy_c.c -L. -lpyfboxlib

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
