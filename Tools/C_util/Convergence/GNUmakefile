
TOP = $(AMREX_HOME)
USE_EB =FALSE
#
# Variables for the user to set ...
#
PRECISION     = DOUBLE
DEBUG	      = TRUE
DEBUG	      = FALSE
DIM	      = 3
#DIM	       = 2
#DIM           = 1
COMP          = g++
#USE_MPI       = TRUE 
USE_MPI       = FALSE
#
# Base name of the executable.
#
#EBASE = Add
### this has no main EBASE = ComputeAmrDataNorms
EBASE = ComparePlotfiles
#EBASE = DiffFab
#EBASE = DiffSameDomainRefined
#EBASE = DiffSameDomainRefinedComposite
#EBASE = DiffSameDomainRefinedStag
#EBASE = DiffSameGrid2
#EBASE = DiffSameGrid
#EBASE = DiffSameGridRefined
#EBASE = DiffUniform
#EBASE = PltFileNormB
#EBASE = PltFileNorm
#EBASE = PltFileScalConvRate
#EBASE  = RichardsonConvergenceTest

CEXE_sources += $(EBASE).cpp

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

DEFINES += -DBL_NOLINEVALUES -DBL_PARALLEL_IO

include Make.package
#include ../Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
#include $(AMREX_HOME)/Src/Amr/Make.package
include $(AMREX_HOME)/Src/Extern/amrdata/Make.package

INCLUDE_LOCATIONS += .
INCLUDE_LOCATIONS += ..
vpathdir += ..

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
vpathdir += $(AMREX_HOME)/Src/Base

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Boundary
vpathdir += $(AMREX_HOME)/Src/Boundary

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Amr
vpathdir += $(AMREX_HOME)/Src/Amr

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata
vpathdir += $(AMREX_HOME)/Src/Extern/amrdata

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/ProfParser
vpathdir += $(AMREX_HOME)/Src/Extern/ProfParser

CEXE_sources += AMReX_BCRec.cpp
CEXE_headers += AMReX_BCRec.H AMReX_BC_TYPES.H

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

#include $(TOP)/Tools/GNUMake/Make.defs ./Make.package

#ifeq ($(USE_ARRAYVIEW),TRUE)
#  DEFINES += -DBL_USE_ARRAYVIEW
#  DEFINES += -DBL_ARRAYVIEW_TAGBOX
#endif





