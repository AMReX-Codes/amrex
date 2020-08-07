AMREX_HOME = ../../

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP    = gcc

USE_PARTICLES = TRUE

PRECISION = DOUBLE

USE_MPI   = TRUE
USE_OMP   = TRUE

###################################################

EBASE     = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
