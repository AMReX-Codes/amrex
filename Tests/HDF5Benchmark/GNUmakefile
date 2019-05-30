AMREX_HOME ?= ../../

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP    = gnu 

PRECISION = DOUBLE

USE_MPI   = TRUE
USE_OMP   = FALSE

TINY_PROFILE = TRUE

###################################################

EBASE     = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

HDF5_HOME   = /home/khl7265/.local
USE_HDF5 = TRUE
ifeq ($(USE_HDF5), TRUE)
DEFINES += -DBL_HDF5 
INCLUDE_LOCATIONS += $(HDF5_HOME)/include
LIBRARIES         += -lhdf5 -lz -ldl -L$(HDF5_HOME)/lib
endif
