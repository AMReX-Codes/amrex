AMREX_HOME ?= ../../../

USE_MPI      = FALSE
USE_OMP      = FALSE
TINY_PROFILE = FALSE
DIM          = 3

USE_CUDA = TRUE

# Also build the MLMG lifetime test.
DEFINES += -DAMREX_TEST_LINEAR_SOLVERS

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

Pdirs := Base Boundary LinearSolvers
Ppack += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)
include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
