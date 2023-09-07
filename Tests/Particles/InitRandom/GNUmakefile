AMREX_HOME = ../../../

# DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP    = gnu

PRECISION = DOUBLE

USE_MPI   = TRUE
MPI_THREAD_MULTIPLE = FALSE

USE_OMP   = FALSE

TINY_PROFILE = TRUE

USE_PARTICLES = TRUE

###################################################

EBASE     = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
