BOXLIB_HOME ?= ../../

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 3

COMP    = gcc
FCOMP   = gfortran


PRECISION = DOUBLE
USE_MPI   = FALSE
USE_OMP   = FALSE
EBASE     = main


include ./Make.package
include $(BOXLIB_HOME)/Tools/C_mk/Make.defs
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
