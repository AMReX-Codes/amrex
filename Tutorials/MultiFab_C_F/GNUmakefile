# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME ?= ../..

DEBUG     = TRUE
USE_MPI   = TRUE
USE_OMP   = FALSE
PROFILE   = FALSE
COMP      = g++
FCOMP     = gfortran
DIM       = 3
PRECISION = DOUBLE
EBASE     = main

include ./Make.package
include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
vpathdir += $(BOXLIB_HOME)/Src/C_BaseLib

include $(BOXLIB_HOME)/Src/F_BaseLib/FParallelMG.mak
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib
vpathdir += $(BOXLIB_HOME)/Src/F_BaseLib

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules

#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)
