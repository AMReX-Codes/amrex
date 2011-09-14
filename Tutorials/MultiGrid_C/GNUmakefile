# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
BOXLIB_HOME ?= ../../..

# where the hypre library is
HYPRE_DIR ?= $(HOME)/lib/hypre

PRECISION = DOUBLE

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 3

COMP =gcc
FCOMP=gfortran

# USE_MPI must be TRUE because of hypre
USE_MPI=TRUE

EBASE = main

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

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

CEXE_sources += BCRec.cpp 
CEXE_headers += BCRec.H BC_TYPES.H
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib
vpathdir += $(BOXLIB_HOME)/Src/C_AMRLib

include $(BOXLIB_HOME)/Src/LinearSolvers/C_to_F_MG/Make.package
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/C_to_F_MG
vpathdir          += $(BOXLIB_HOME)/Src/LinearSolvers/C_to_F_MG

include $(BOXLIB_HOME)/Src/LinearSolvers/F_MG/FParallelMG.mak
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/F_MG
vpathdir          += $(BOXLIB_HOME)/Src/LinearSolvers/F_MG

include $(BOXLIB_HOME)/Src/F_BaseLib/FParallelMG.mak
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib
vpathdir          += $(BOXLIB_HOME)/Src/F_BaseLib

# hypre driver
include HypreABecLap/Make.package
INCLUDE_LOCATIONS += HypreABecLap
vpathdir          += HypreABecLap

# hypre library
INCLUDE_LOCATIONS += $(HYPRE_DIR)/include
vpathdir          += $(HYPRE_DIR)/include
LIBRARY_LOCATIONS += $(HYPRE_DIR)/lib
LIBRARIES += -lHYPRE

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
