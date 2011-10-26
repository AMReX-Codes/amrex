BOXLIB_HOME ?= ../../..

TOP = $(BOXLIB_HOME)
#
# Variables for the user to set ...
#
PRECISION     = DOUBLE
DEBUG	      = TRUE
DEBUG	      = FALSE
DIM	      = 2
COMP          = g++
FCOMP         = gfortran
USE_MPI       = TRUE
#
# Base name of the executable.
#
EBASE = DiffUniform
#EBASE = PltFileNorm
#EBASE = PltFileScalConvRate
#EBASE = DiffFab
#EBASE = DiffSameGrid
#EBASE = PltFileNorm
#EBASE = PltFileNormB
#EBASE = PltFileScalConvRate
#EBASE = DiffSameGrid
#EBASE = DiffSameGrid2
#EBASE = DiffUniform
#EBASE = DiffSameGridRefined
#EBASE = DiffFab

CEXE_sources += $(EBASE).cpp

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

DEFINES += -DBL_NOLINEVALUES -DBL_PARALLEL_IO

include Make.package
include ../Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package
include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
#include $(BOXLIB_HOME)/Src/C_AMRLib/Make.package
include $(BOXLIB_HOME)/Src/Extern/amrdata/Make.package

INCLUDE_LOCATIONS += .
INCLUDE_LOCATIONS += ..
vpathdir += ..

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
vpathdir += $(BOXLIB_HOME)/Src/C_BaseLib

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BoundaryLib
vpathdir += $(BOXLIB_HOME)/Src/C_BoundaryLib

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib
vpathdir += $(BOXLIB_HOME)/Src/C_AMRLib

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Extern/amrdata
vpathdir += $(BOXLIB_HOME)/Src/Extern/amrdata

CEXE_sources += BCRec.cpp
CEXE_headers += BCRec.H BC_TYPES.H

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules

#include $(TOP)/Tools/C_mk/Make.defs ./Make.package

#ifeq ($(USE_ARRAYVIEW),TRUE)
#  DEFINES += -DBL_USE_ARRAYVIEW
#  DEFINES += -DBL_ARRAYVIEW_TAGBOX
#endif





