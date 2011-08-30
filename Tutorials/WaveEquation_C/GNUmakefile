# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
BOXLIB_HOME ?= /path/to/BoxLib

PRECISION = DOUBLE

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 2

COMP =gcc
FCOMP=gfortran

USE_MPI=FALSE
USE_MPI=TRUE

EBASE = main

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

include ./Make.package

include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BoundaryLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib

vpathdir += $(BOXLIB_HOME)/Src/C_BoundaryLib
vpathdir += $(BOXLIB_HOME)/Src/C_BaseLib
vpathdir += $(BOXLIB_HOME)/Src/C_AMRLib

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
