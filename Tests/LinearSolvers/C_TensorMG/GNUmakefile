AMREX_HOME := ../../..

PRECISION = DOUBLE

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 2
DIM	= 3

COMP    = gnu

USE_MPI = FALSE

EBASE = main

include $(AMREX_HOME)/Tools/C_mk/Make.defs

CXXFLAGS +=
CXXOPTF  +=
CXXDEBF  +=

include ./Make.package

include $(AMREX_HOME)/Src/LinearSolvers/C_TensorMG/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(AMREX_HOME)/Src/C_BaseLib/Make.package
include $(AMREX_HOME)/Src/C_BoundaryLib/Make.package

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/LinearSolvers/C_TensorMG
vpathdir          += $(AMREX_HOME)/Src/LinearSolvers/C_TensorMG

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/C_BoundaryLib
vpathdir          += $(AMREX_HOME)/Src/C_BoundaryLib

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/C_BaseLib
vpathdir          += $(AMREX_HOME)/Src/C_BaseLib

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/C_AMRLib
vpathdir          += $(AMREX_HOME)/Src/C_AMRLib

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/LinearSolvers/C_CellMG
vpathdir          += $(AMREX_HOME)/Src/LinearSolvers/C_CellMG

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

include $(AMREX_HOME)/Tools/C_mk/Make.rules


