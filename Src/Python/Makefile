MPI  := t
COMP := gfortran

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

CC       := mpicc
FC       := mpif90
F90      := mpif90
F90FLAGS += -fPIC
CFLAGS   += -fPIC

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

f90sources += blobjects.f90
f90sources += fboxlib.f90
f90sources += boxlib_numpy.f90
VPATH_LOCATIONS += src

all: libpyfboxlib.so libpycboxlib.so

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
include $(BOXLIB_HOME)/Src/Python/GMakerules.mak
