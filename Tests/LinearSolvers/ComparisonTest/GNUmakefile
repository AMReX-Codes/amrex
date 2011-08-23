# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
BOXLIB_HOME ?= ../../..

USE_F90_MG := TRUE

PRECISION = DOUBLE

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP =gcc
FCOMP=gfortran

USE_MPI=TRUE
USE_MPI=FALSE

EBASE = main

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

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

# We need these from C_AMRLib 

CEXE_sources += BCRec.cpp 
CEXE_headers += BCRec.H BC_TYPES.H
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib
vpathdir += $(BOXLIB_HOME)/Src/C_AMRLib

ifdef USE_F90_MG
  DEFINES += -DUSE_F90_MG

  include $(BOXLIB_HOME)/Src/LinearSolvers/C_to_F_MG/Make.package
  INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/C_to_F_MG
  vpathdir          += $(BOXLIB_HOME)/Src/LinearSolvers/C_to_F_MG

  include $(BOXLIB_HOME)/Src/LinearSolvers/F_MG/FParallelMG.mak
  INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/F_MG
  vpathdir          += $(BOXLIB_HOME)/Src/LinearSolvers/F_MG

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

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
