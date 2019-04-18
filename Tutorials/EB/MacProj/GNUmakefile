USE_EB = TRUE
DEBUG = FALSE
USE_MPI  = TRUE
USE_OMP  = FALSE

USE_HYPRE = FALSE

COMP = gnu

DIM = 2

AMREX_HOME ?= ../../..

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

Pdirs := Base Boundary AmrCore
Pdirs += EB
Pdirs += LinearSolvers/C_CellMG LinearSolvers/MLMG

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

