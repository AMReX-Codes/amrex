PBOXLIB_HOME := ../..
FBOXLIB_HOME = ../../../fParallel
FBOXLIB_HOME =

PRECISION = DOUBLE
PROFILE = TRUE
PROFILE = FALSE
DEBUG	= TRUE
DEBUG	= FALSE
DIM	= 2
DIM	= 3

COMP=g++
FCOMP=Intel

USE_MPI=FALSE
USE_MPI=TRUE
USE_THREADS=TRUE
USE_THREADS=FALSE

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

INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/bndrylib
vpathdir += $(PBOXLIB_HOME)/bndrylib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/BoxLib
vpathdir += $(PBOXLIB_HOME)/BoxLib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/iamrlib
vpathdir += $(PBOXLIB_HOME)/iamrlib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/mglib
vpathdir += $(PBOXLIB_HOME)/mglib

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
CEXE_headers += BCRec.H BC_TYPES.H
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/amrlib
vpathdir += $(PBOXLIB_HOME)/amrlib

# Needed from AmrVis
#CEXE_headers += AmrData.H
#INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/pAmrvis
#vpathdir += $(PBOXLIB_HOME)/pAmrvis

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(PBOXLIB_HOME)/mk/Make.rules
