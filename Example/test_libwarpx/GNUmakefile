AMREX_HOME ?= ../../../amrex

DEBUG	= FALSE
#DEBUG	= TRUE

DIM	= 3

COMP    = gcc

USE_MPI   = TRUE
USE_OMP   = FALSE

LIBRARIES += -lwarpx

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += main.cpp
CEXE_headers += WarpXWrappers.h

INCLUDE_LOCATIONS += ../../Source

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
