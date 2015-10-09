BOXLIB_HOME ?= ../../

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 3

COMP    = gcc
FCOMP   = gfortran

USE_PARTICLES = TRUE

PRECISION = DOUBLE
USE_MPI   = FALSE
USE_OMP   = FALSE
EBASE     = main

include ./Make.package
include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

# These are the directories in BoxLib
Pdirs   := C_BaseLib C_AMRLib C_BoundaryLib 
Ppack   += $(foreach dir, $(Pdirs), $(BOXLIB_HOME)/Src/$(dir)/Make.package)
Plocs   += $(foreach dir, $(Pdirs), $(BOXLIB_HOME)/Src/$(dir))

include $(Ppack)
INCLUDE_LOCATIONS += $(Plocs)
VPATH_LOCATIONS   += $(Plocs)

vpath %.c   : . $(VPATH_LOCATIONS)
vpath %.h   : . $(VPATH_LOCATIONS)
vpath %.cpp : . $(VPATH_LOCATIONS)
vpath %.H   : . $(VPATH_LOCATIONS)
vpath %.F   : . $(VPATH_LOCATIONS)
vpath %.f   : . $(VPATH_LOCATIONS)
vpath %.f90 : . $(VPATH_LOCATIONS)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
