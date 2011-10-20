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

include $(BOXLIB_HOME)/Src/Python/GPackage.mak

# add extra modules
#pybl_sources +=

all: pyfboxlib.so

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
include $(BOXLIB_HOME)/Src/Python/GMakerules.mak
