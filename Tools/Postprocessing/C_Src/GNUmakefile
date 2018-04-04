AMREX_HOME ?= ../../..

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP    = gcc

PRECISION = DOUBLE

USE_MPI   = TRUE
USE_OMP   = FALSE

###################################################

EBASE     = particle_compare
#EBASE     = PlotfileToMatLab

# If NEEDS_F90_SRC=TRUE, look for ${EBASE}_F.F90
NEEDS_F90_SRC = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

ifeq ($(EBASE), particle_compare)
  executable = $(EBASE).exe
endif

include ./Make.package
include $(AMREX_HOME)/Tools/GNUMake/Make.rules
