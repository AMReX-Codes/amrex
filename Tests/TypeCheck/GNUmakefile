AMREX_HOME = ../../

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 3

COMP    = gcc

PRECISION = DOUBLE

USE_MPI   = TRUE
USE_OMP   = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

default: typecheck
	@echo "\n7 errors should be found.  No errros should be found in functions starting with pass_."

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
