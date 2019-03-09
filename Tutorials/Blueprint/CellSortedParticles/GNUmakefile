AMREX_HOME ?= ../../../

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP    = gcc

TINY_PROFILE  = TRUE
USE_PARTICLES = TRUE

PRECISION     = DOUBLE

USE_MPI       = FALSE
USE_OMP       = FALSE

###################################################

EBASE     = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

ifeq ($(USE_CONDUIT),TRUE)
    include $(AMREX_HOME)/Src/Extern/Conduit/Make.package
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
