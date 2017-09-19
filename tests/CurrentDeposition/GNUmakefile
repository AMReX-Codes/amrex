AMREX_HOME ?= ../../../amrex
PICSAR_HOME ?= ../../../picsar

USE_PARTICLES = TRUE

DEBUG     = FALSE
USE_MPI   = TRUE
USE_OMP   = TRUE
PROFILE   = FALSE
COMP      = gnu
DIM       = 3
PRECISION = DOUBLE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(PICSAR_HOME)/src/Make.package

default: $(executable)
	@echo SUCCESS

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
