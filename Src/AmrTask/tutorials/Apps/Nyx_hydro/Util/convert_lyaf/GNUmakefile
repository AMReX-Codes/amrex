# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../amrex
HDF5_DIR   ?=

DEBUG     = FALSE
USE_MPI   = TRUE
USE_OMP   = FALSE
PROFILE   = TRUE

COMP      = intel

PRECISION = DOUBLE
DEFINES += -DBL_NOLINEVALUES

DIM       = 3
EBASE     = convert

include ./Make.package
include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Extern/amrdata/Make.package

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Boundary
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Amr
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata

# hdf5 stuff
INCLUDE_LOCATIONS += $(HDF5_DIR)/include
LIBRARY_LOCATIONS += $(HDF5_DIR)/lib
LIBRARIES += -lhdf5 #$(GPERFTOOLS_TCMALLOC)

vpathdir += $(AMREX_HOME)/Src/Boundary
vpathdir += $(AMREX_HOME)/Src/Base
vpathdir += $(AMREX_HOME)/Src/Amr
vpathdir += $(AMREX_HOME)/Src/Extern/amrdata

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)
	@echo SUCCESS

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
