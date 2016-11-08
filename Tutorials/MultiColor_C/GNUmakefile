BOXLIB_HOME ?= ../..

DEBUG	= TRUE

DIM	= 2

COMP =gcc

USE_MPI=TRUE

EBASE = main

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

CEXE_sources += $(EBASE).cpp

include Make.package
include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
include $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

all: $(executable)
	@echo SUCCESS

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
