
AMREX_HOME ?= ../..

# compilation options
COMP    = g++

USE_MPI   = TRUE
USE_UPCXX = TRUE

USE_OMP   = FALSE

PROFILE   = FALSE

PRECISION = DOUBLE

DEBUG     = FALSE

DIM      = 3

Bpack := ./Make.package

EBASE = SMC

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

all: $(executable)
	@echo SUCCESS

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules


#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)

