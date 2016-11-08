# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME ?= ../..

DEBUG     = FALSE
USE_MPI   = FALSE
USE_OMP   = FALSE
PROFILE   = FALSE
COMP      = g++
DIM       = 2
PRECISION = DOUBLE
EBASE     = main

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

include ./Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
