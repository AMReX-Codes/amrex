# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME ?= ../..

DEBUG     = TRUE

USE_MPI   = FALSE
USE_OMP   = FALSE
USE_MPI3  = FALSE

MEM_PROFILE  = FALSE
TINY_PROFILE = FALSE

COMP      = g++
FCOMP     = gfortran

DIM       = 3

EBASE     = main

include ./Make.package
include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package
include $(BOXLIB_HOME)/Src/F_Interfaces/BaseLib/Make.package

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/F_Interfaces/BaseLib

VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_Interfaces/BaseLib

vpath %.c   . $(VPATH_LOCATIONS)
vpath %.h   . $(VPATH_LOCATIONS)
vpath %.cpp . $(VPATH_LOCATIONS)
vpath %.H   . $(VPATH_LOCATIONS)
vpath %.F   . $(VPATH_LOCATIONS)
vpath %.f   . $(VPATH_LOCATIONS)
vpath %.f90 . $(VPATH_LOCATIONS)
vpath %.F90 . $(VPATH_LOCATIONS)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules

