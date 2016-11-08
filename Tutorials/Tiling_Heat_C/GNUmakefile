# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME ?= ../..

DEBUG     = FALSE
USE_MPI   = FALSE
USE_OMP   = TRUE
PROFILE   = FALSE
COMP      = Intel
DIM       = 3
PRECISION = DOUBLE
EBASE     = Heat

include ./Make.package
include $(BOXLIB_HOME)/Tools/C_mk/Make.defs
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
