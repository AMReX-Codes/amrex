PBOXLIB_HOME := ../..
FBOXLIB_HOME =
FBOXLIB_HOME = ../../../fParallel

PRECISION = DOUBLE
PROFILE = TRUE
PROFILE = FALSE
DEBUG	= FALSE
DEBUG	= TRUE
DIM	= 3
DIM	= 2
COMP = KCC
COMP = g++

COMP=Intel
FCOMP=Intel

USE_MPI=TRUE
USE_MPI=FALSE
USE_THREADS=TRUE
USE_THREADS=FALSE

STRICTLY=TRUE

USE_ARRAYVIEW = TRUE
USE_ARRAYVIEW = FALSE

EBASE = main

include $(PBOXLIB_HOME)/mk/Make.defs

# multigrid test sources
CEXE_sources += $(EBASE).cpp

include Make.package
include $(PBOXLIB_HOME)/bndrylib/Make.package
include $(PBOXLIB_HOME)/mglib/Make.package
include $(PBOXLIB_HOME)/BoxLib/Make.package

#include $(PBOXLIB_HOME)/util/Make.package

INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/bndrylib
vpathdir += $(PBOXLIB_HOME)/bndrylib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/BoxLib
vpathdir += $(PBOXLIB_HOME)/BoxLib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/iamrlib
vpathdir += $(PBOXLIB_HOME)/iamrlib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/mglib
vpathdir += $(PBOXLIB_HOME)/mglib

#INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/util
#vpathdir += $(PBOXLIB_HOME)/util

ifdef FBOXLIB_HOME
  DEFINES += -DMG_USE_FBOXLIB

  include $(PBOXLIB_HOME)/MGT_Solver/Make.package
  INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/MGT_Solver
  vpathdir += $(PBOXLIB_HOME)/MGT_Solver

  include $(FBOXLIB_HOME)/mg/FParallelMG.mak
  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/mg/
  Fdirs   := boxlib mg extern/SPARSKIT extern/LAPACK
  Flocs   := $(foreach dir, $(Fdirs), $(FBOXLIB_HOME)/$(dir))
  vpathdir += $(Flocs)
endif



# We need these from amrlib 

CEXE_sources += BCRec.cpp 
CEXE_headers += BCRec.H BC_Types.H
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/amrlib
vpathdir += $(PBOXLIB_HOME)/amrlib

# Needed from AmrVis
CEXE_headers += AmrData.H
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/pAmrvis
vpathdir += $(PBOXLIB_HOME)/pAmrvis

ifeq ($(DEBUG),TRUE)
  ifeq ($(USE_ARRAYVIEW),TRUE)
    ARRAYVIEWDIR  = $(PBOXLIB_HOME)/util
    INCLUDE_LOCATIONS += $(ARRAYVIEWDIR)
    CEXE_headers += DatasetClient.H ArrayView.H
    CEXE_sources += DatasetClient.cpp
    DEFINES += -DBL_USE_ARRAYVIEW
    #DEFINES += -DBL_ARRAYVIEW_TAGBOX
  endif
  CEXE_headers += TV_TempWrite.H
endif

CXXFLAGS +=
CXXOPTF +=
CXXDEBF +=

CFLAGS +=
COPTF +=
CDEBF +=

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(PBOXLIB_HOME)/mk/Make.rules
