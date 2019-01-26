AMREX_HOME ?= ../../../amrex

DEBUG	= TRUE

DIM	= 3

#COMP    = intel  # gcc
COMP    = gcc

TINY_PROFILE = FALSE
USE_PARTICLES = FALSE

PRECISION = DOUBLE

USE_MPI   = FALSE
USE_OMP   = FALSE

###################################################

EBASE     = sliceutils

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
