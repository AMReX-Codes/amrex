BOXLIB_HOME := ../..
FBOXLIB_HOME =
FBOXLIB_HOME = ../../../fParallel

PRECISION = DOUBLE
PROFILE = FALSE
DEBUG	= FALSE
DEBUG	= TRUE
DIM	= 3
DIM	= 2

COMP=g++
FCOMP=gfortran

USE_MPI=FALSE
USE_ARRAYVIEW = FALSE

EBASE = main

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

# multigrid test sources
CEXE_sources += $(EBASE).cpp

include Make.package
include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
include $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/bndrylib
vpathdir += $(BOXLIB_HOME)/bndrylib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/BoxLib
vpathdir += $(BOXLIB_HOME)/BoxLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/iamrlib
vpathdir += $(BOXLIB_HOME)/iamrlib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/mglib
vpathdir += $(BOXLIB_HOME)/mglib

ifdef FBOXLIB_HOME
  DEFINES += -DMG_USE_FBOXLIB

  include $(BOXLIB_HOME)/MGT_Solver/Make.package
  INCLUDE_LOCATIONS += $(BOXLIB_HOME)/MGT_Solver
  vpathdir += $(BOXLIB_HOME)/MGT_Solver

  include $(FBOXLIB_HOME)/mg/source/FParallelMG.mak
  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/mg/source
  Fdirs   := boxlib mg/source extern/SPARSKIT extern/LAPACK extern/BLAS
  include $(FBOXLIB_HOME)/extern/BLAS/Make.package
  Flocs   := $(foreach dir, $(Fdirs), $(FBOXLIB_HOME)/$(dir))
  vpathdir += $(Flocs)
endif

# We need these from amrlib 

CEXE_sources += BCRec.cpp 
CEXE_headers += BCRec.H BC_TYPES.H
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/amrlib
vpathdir += $(BOXLIB_HOME)/amrlib

# Needed from AmrVis
#CEXE_headers += AmrData.H
#INCLUDE_LOCATIONS += $(BOXLIB_HOME)/pAmrvis
#vpathdir += $(BOXLIB_HOME)/pAmrvis

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(BOXLIB_HOME)/mk/Make.rules
