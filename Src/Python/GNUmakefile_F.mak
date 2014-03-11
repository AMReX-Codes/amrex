#
# Makefile for FORTRAN version of PyBoxLib
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

include $(BOXLIB_HOME)/Src/Python/GPackage.mak

all: $(PYFBOXLIB)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
include $(BOXLIB_HOME)/Src/Python/GMakerules.mak



