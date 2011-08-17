BOXLIB_HOME := ../..

PROFILE = FALSE
PRECISION = DOUBLE
DEBUG	= FALSE
DEBUG	= TRUE
DIM	= 2
DIM	= 3
COMP    = Intel
FCOMP   = Intel
USE_MPI = FALSE

USE_ARRAYVIEW = FALSE

EBASE = main

include $(BOXLIB_HOME)/mk/Make.defs

CXXFLAGS +=
CXXOPTF  +=
CXXDEBF  +=

include Make.package

include $(BOXLIB_HOME)/Src/LinearSolvers/C_TensorMG/Make.package
include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
include $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/tensorMG
vpathdir += $(BOXLIB_HOME)/tensorMG
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/bndrylib
vpathdir += $(BOXLIB_HOME)/bndrylib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/BoxLib
vpathdir += $(BOXLIB_HOME)/BoxLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/iamrlib
vpathdir += $(BOXLIB_HOME)/iamrlib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/mglib
vpathdir += $(BOXLIB_HOME)/mglib

# We need these from amrlib 

CEXE_sources += BCRec.cpp 
CEXE_headers += BCRec.H BC_TYPES.H
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/amrlib
vpathdir += $(BOXLIB_HOME)/amrlib

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(BOXLIB_HOME)/mk/Make.rules


