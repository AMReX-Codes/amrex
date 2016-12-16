# AMREX_HOME defines the directory in which we will find all the BoxLib code
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../..

DEBUG     = TRUE

USE_MPI   = TRUE
USE_OMP   = FALSE
USE_MPI3  = FALSE

MEM_PROFILE  = FALSE
TINY_PROFILE = FALSE

COMP      = g++

DIM       = 3

EBASE     = main

include ./Make.package
include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/F_Interfaces/BaseLib/Make.package

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/F_Interfaces/BaseLib

VPATH_LOCATIONS += $(AMREX_HOME)/Src/Base
VPATH_LOCATIONS += $(AMREX_HOME)/Src/F_Interfaces/BaseLib

vpath %.c   . $(VPATH_LOCATIONS)
vpath %.h   . $(VPATH_LOCATIONS)
vpath %.cpp . $(VPATH_LOCATIONS)
vpath %.H   . $(VPATH_LOCATIONS)
vpath %.F   . $(VPATH_LOCATIONS)
vpath %.f   . $(VPATH_LOCATIONS)
vpath %.f90 . $(VPATH_LOCATIONS)
vpath %.F90 . $(VPATH_LOCATIONS)

all: $(executable)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

