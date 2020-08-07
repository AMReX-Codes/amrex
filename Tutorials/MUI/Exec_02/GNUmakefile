# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME = ../../../
MUI_HOME ?= ../../../../MUI/

DEBUG     = FALSE
USE_MPI   = TRUE
USE_OMP   = FALSE
COMP      = gnu
DIM       = 2

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ../Source_02/Make.package
VPATH_LOCATIONS  += ../Source_02
INCLUDE_LOCATIONS += ../Source_02
INCLUDE_LOCATIONS += $(MUI_HOME)

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

