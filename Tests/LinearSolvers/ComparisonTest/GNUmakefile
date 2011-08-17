BOXLIB_HOME := /path/to/BoxLib

PRECISION = DOUBLE

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP =gcc
FCOMP=gfortran

USE_MPI=TRUE
USE_MPI=FALSE

EBASE = main

include $(BOXLIB_HOME)/mk/Make.defs

# multigrid test sources
CEXE_sources += $(EBASE).cpp

include Make.package
include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
include $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BoundaryLib
vpathdir += $(BOXLIB_HOME)/Src/C_BoundaryLib

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
vpathdir += $(BOXLIB_HOME)/Src/C_BaseLib

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG
vpathdir += $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG

# We need these from amrlib 
CEXE_sources += BCRec.cpp 
CEXE_headers += BCRec.H BC_TYPES.H
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/amrlib
vpathdir += $(BOXLIB_HOME)/amrlib

ifdef FBOXLIB_HOME
  DEFINES += -DMG_USE_FBOXLIB

  include $(BOXLIB_HOME)/Src/LinearSolvers/C_to_F_MG/Make.package
  INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/C_to_F_MG
  vpathdir          += $(BOXLIB_HOME)/Src/LinearSolvers/C_to_F_MG

  include $(BOXLIB_HOME)/Src/LinearSolvers/F_MG/FParallelMG.mak
  INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/F_MG
  vpathdir          += $(BOXLIB_HOME)/Src/LinearSolvers/F_MG

  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/extern/SPARSKIT
  vpathdir          += $(FBOXLIB_HOME)/extern/SPARSKIT

  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/extern/LAPACK
  vpathdir          += $(FBOXLIB_HOME)/extern/LAPACK

  include $(FBOXLIB_HOME)/extern/BLAS/Make.package
  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/extern/BLAS
  vpathdir          += $(FBOXLIB_HOME)/extern/BLAS

  INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib
  vpathdir          += $(BOXLIB_HOME)/Src/F_BaseLib
endif

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(BOXLIB_HOME)/mk/Make.rules
