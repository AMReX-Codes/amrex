AMREX_HOME ?= ../..

DEBUG	= TRUE

DIM	= 2

COMP =gcc

USE_MPI=TRUE

EBASE = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(AMREX_HOME)/Src/Base/Make.package

all: $(executable)
	@echo SUCCESS

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
