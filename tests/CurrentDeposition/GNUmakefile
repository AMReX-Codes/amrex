BOXLIB_HOME ?= ../../../BoxLib
PICSAR_HOME ?= ../../../picsar

USE_PARTICLES = TRUE

DEBUG     = FALSE
USE_MPI   = TRUE
USE_OMP   = FALSE
PROFILE   = FALSE
COMP      = gnu
DIM       = 3
PRECISION = DOUBLE

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

include ./Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package
include $(PICSAR_HOME)/src/Make.package

default: $(executable)
	@echo SUCCESS

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
