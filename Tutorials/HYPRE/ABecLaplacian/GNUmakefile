AMREX_HOME ?= ../../..

# Set HYPRE_DIR here or as environment variable
# HYPRE_DIR ?=

USE_HYPRE = TRUE

DIM = 3

USE_MPI = TRUE

DEBUG = FALSE

TINY_PROFILE = TRUE

COMP = gcc

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/Extern/HYPRE/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
