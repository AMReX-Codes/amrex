AMREX_HOME ?= ../../../

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP    = pgi

TINY_PROFILE = TRUE
USE_PARTICLES = TRUE

PRECISION = DOUBLE

USE_MPI     = FALSE
USE_OMP     = FALSE
USE_CUDA    = TRUE
USE_CUDA_UM = TRUE

###################################################

EBASE     = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
