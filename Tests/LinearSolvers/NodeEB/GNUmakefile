
DEBUG = FALSE

TEST = TRUE
USE_ASSERTION = TRUE

USE_EB = TRUE

USE_MPI  = TRUE
USE_OMP  = FALSE

USE_HYPRE  = FALSE
USE_PETSC  = FALSE

COMP = gnu

DIM = 3

ifeq ($(DIM),3)
  USE_ALGOIM = TRUE
else
  USE_ALGOIM = FALSE
endif

AMREX_HOME ?= ../../..

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include ./Make.package

Pdirs := Base Boundary AmrCore
Pdirs += EB
Pdirs += LinearSolvers/C_CellMG LinearSolvers/MLMG

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

