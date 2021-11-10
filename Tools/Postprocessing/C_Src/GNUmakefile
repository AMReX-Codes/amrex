AMREX_HOME ?= ../../..

DEBUG	= TRUE
#DEBUG	= FALSE

DIM	= 3
DIM	= 2

COMP    = gcc

PRECISION = DOUBLE

USE_MPI   = FALSE
USE_OMP   = FALSE

###################################################

EBASE     = particle_compare
#EBASE     = PlotfileToMatLab
#EBASE     = PlotfileToTurb
#EBASE     = WritePlotfileToASCII
#EBASE     = MultiFabToMatLab
#EBASE     = PtwisePltTransform
#EBASE     = cube_extract
#EBASE     = HorizontalAvg
#EBASE     = IntegrateComp

# If NEEDS_f90_SRC=TRUE, look for ${EBASE}_nd.f90
ifeq ($(EBASE), PtwisePltTransform)
  NEEDS_f90_SRC = TRUE
else ifeq ($(EBASE), PlotfileToTurb)
  NEEDS_f90_SRC = TRUE
else
  NEEDS_f90_SRC = FALSE
endif

ifeq ($(EBASE),particle_compare)
  DIM = 3
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

ifeq ($(EBASE), particle_compare)
  executable = $(EBASE).exe
endif

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

clean::
	$(SILENT) $(RM) particle_compare.exe
