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

all: $(executable)

include $(AMREX_HOME)/Tools/C_mk/Make.rules


