AMREX_HOME ?= ../..

DEBUG	= TRUE

DIM	= 2

COMP =gcc

USE_MPI=TRUE

EBASE = main

include $(AMREX_HOME)/Tools/C_mk/Make.defs

CEXE_sources += $(EBASE).cpp

include Make.package
include $(AMREX_HOME)/Src/C_BoundaryLib/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(AMREX_HOME)/Src/C_BaseLib/Make.package

all: $(executable)
	@echo SUCCESS

include $(AMREX_HOME)/Tools/C_mk/Make.rules
