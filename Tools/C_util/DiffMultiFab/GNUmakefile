AMREX_HOME ?= ../../..

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 2
DIM	= 3

COMP    = gcc

PRECISION = DOUBLE

BL_NO_FORT = TRUE

USE_MPI   = FALSE
USE_OMP   = FALSE

TEST=TRUE
USE_ASSERTION=TRUE

###################################################

EBASE     = diffmultifab

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += ${EBASE}.cpp

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
