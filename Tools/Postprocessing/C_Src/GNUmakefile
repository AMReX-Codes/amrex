AMREX_HOME ?= ../../..

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3
DIM	= 2

COMP    = gcc

PRECISION = DOUBLE

USE_MPI   = TRUE
USE_OMP   = FALSE

###################################################

EBASE     = particle_compare
#EBASE     = PlotfileToMatLab
#EBASE     = MultiFabToMatLab
#EBASE     = PtwisePltTransform
#EBASE     = PlotfileGradient
#EBASE     = cube_extract

# If NEEDS_f90_SRC=TRUE, look for ${EBASE}_nd.f90
NEEDS_f90_SRC = FALSE
#NEEDS_f90_SRC = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

ifeq ($(EBASE), particle_compare)
  executable = $(EBASE).exe
endif

include ./Make.package

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
