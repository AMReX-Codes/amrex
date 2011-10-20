MPI  := t
COMP := gfortran

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

CC       := mpicc
FC       := mpif90
F90      := mpif90
F90FLAGS += -fPIC		# needed for dynamic loading with Python
CFLAGS   += -fPIC

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

#include $(BOXLIB_HOME)/Src/LinearSolvers/F_MG/GPackage.mak
#VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/F_MG

PyBL := src/blobjects.f90 src/pyfboxlib.pyf src/fboxlib.pyf

pyboxlib/pyfboxlib.so: $(objects) $(PyBL)
	@echo Running f2py...
	@f2py --quiet --fcompiler=gnu95 --f90exec=$(FC) --f90flags="-I $(mdir)" \
		-c src/pyfboxlib.pyf src/blobjects.f90 src/boxlib_numpy.f90 src/fboxlib.f90 $(objects)
	@mv pyfboxlib.so pyboxlib

# PyBoxLib
src/blobjects.f90: src/blobjects.py
	cd src && python blobjects.py

src/fboxlib.pyf: src/fboxlib.f90
	cd src && f2py --overwrite-signature -h fboxlib.pyf fboxlib.f90

src/pyfboxlib.pyf: src/pyfboxlib.m4 src/boxlib_numpy.c src/fboxlib.pyf
	cd src && m4 pyfboxlib.m4 > pyfboxlib.pyf

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak

