BOXLIB_HOME ?= ../..

DEBUG	= TRUE

DIM	= 2

COMP =gcc
FCOMP=gfortran

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

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG4
vpathdir += $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG4

vpath %.c   : . $(vpathdir) $(VPATH_LOCATIONS)
vpath %.h   : . $(vpathdir) $(VPATH_LOCATIONS)
vpath %.cpp : . $(vpathdir) $(VPATH_LOCATIONS)
vpath %.H   : . $(vpathdir) $(VPATH_LOCATIONS)
vpath %.F   : . $(vpathdir) $(VPATH_LOCATIONS)
vpath %.f   : . $(vpathdir) $(VPATH_LOCATIONS)
vpath %.f90 : . $(vpathdir) $(VPATH_LOCATIONS)

all: $(executable)
	@echo SUCCESS

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
