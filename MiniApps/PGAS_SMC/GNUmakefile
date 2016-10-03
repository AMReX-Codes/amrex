
BOXLIB_HOME ?= ../..

# compilation options
COMP    = g++
FCOMP   = gfortran

USE_MPI   = TRUE
USE_UPCXX = TRUE

USE_OMP   = FALSE

PROFILE   = FALSE

PRECISION = DOUBLE

DEBUG     = FALSE

DIM      = 3

Bpack := ./Make.package

EBASE = SMC

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

all: $(executable)
	@echo SUCCESS

include ./Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules


#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)

