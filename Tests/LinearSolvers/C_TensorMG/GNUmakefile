BOXLIB_HOME := ../../..

PRECISION = DOUBLE

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 2
DIM	= 3

COMP    = Intel
FCOMP   = Intel

USE_MPI = FALSE

EBASE = main

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

CXXFLAGS +=
CXXOPTF  +=
CXXDEBF  +=

include ./Make.package

include $(BOXLIB_HOME)/Src/LinearSolvers/C_TensorMG/Make.package
include $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package
include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/C_TensorMG
vpathdir          += $(BOXLIB_HOME)/Src/LinearSolvers/C_TensorMG

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BoundaryLib
vpathdir          += $(BOXLIB_HOME)/Src/C_BoundaryLib

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
vpathdir          += $(BOXLIB_HOME)/Src/C_BaseLib

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib
vpathdir          += $(BOXLIB_HOME)/Src/C_AMRLib

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG
vpathdir          += $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG

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


