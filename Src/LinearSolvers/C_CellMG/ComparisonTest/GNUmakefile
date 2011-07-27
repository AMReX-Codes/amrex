PBOXLIB_HOME := ../..
FBOXLIB_HOME := ../../../fParallel

PRECISION = DOUBLE
PROFILE   = FALSE

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

#COMP = Intel
#FCOMP= Intel

COMP =gcc
FCOMP=gfortran

USE_MPI=TRUE
USE_MPI=FALSE

#GPROF = TRUE

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

INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/mglib
vpathdir += $(PBOXLIB_HOME)/mglib

# We need these from amrlib 
CEXE_sources += BCRec.cpp 
CEXE_headers += BCRec.H BC_TYPES.H
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/amrlib
vpathdir += $(PBOXLIB_HOME)/amrlib

ifdef FBOXLIB_HOME
  DEFINES += -DMG_USE_FBOXLIB

  include $(PBOXLIB_HOME)/MGT_Solver/Make.package
  INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/MGT_Solver
  vpathdir          += $(PBOXLIB_HOME)/MGT_Solver

  include $(FBOXLIB_HOME)/mg/source/FParallelMG.mak
  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/mg/source
  vpathdir          += $(FBOXLIB_HOME)/mg/source

  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/extern/SPARSKIT
  vpathdir          += $(FBOXLIB_HOME)/extern/SPARSKIT

  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/extern/LAPACK
  vpathdir          += $(FBOXLIB_HOME)/extern/LAPACK

  include $(FBOXLIB_HOME)/extern/BLAS/Make.package
  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/extern/BLAS
  vpathdir          += $(FBOXLIB_HOME)/extern/BLAS

  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/boxlib
  vpathdir          += $(FBOXLIB_HOME)/boxlib
endif

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(PBOXLIB_HOME)/mk/Make.rules
