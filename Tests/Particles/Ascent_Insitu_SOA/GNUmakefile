AMREX_HOME = ../../../

DEBUG	= FALSE

DIM	= 3

COMP    = gcc

USE_MPI   = TRUE
USE_OMP   = FALSE
USE_CUDA  = FALSE

TINY_PROFILE = TRUE
USE_PARTICLES = TRUE
USE_ASCENT = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package
include $(AMREX_HOME)/Src/Extern/Conduit/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
