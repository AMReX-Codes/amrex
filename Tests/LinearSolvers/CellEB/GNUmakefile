DEBUG = FALSE
TEST = TRUE
USE_ASSERTION = TRUE

USE_EB = TRUE

USE_MPI  = TRUE
USE_OMP  = FALSE

USE_HYPRE = FALSE
USE_PETSC = FALSE


COMP = gnu

DIM = 2

AMREX_HOME ?= ../../..
ifeq ($(USE_PETSC),TRUE)
#Your PETSC Directory
    PETSC_DIR = /usr/local/petsc
#The Hypre directory for utilizing hypre functions
    HYPRE_DIR = /home/sireeves/hypre/src/hypre
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include ./Make.package

Pdirs := Base Boundary AmrCore
Pdirs += EB EB2
Pdirs += LinearSolvers/C_CellMG LinearSolvers/MLMG

ifeq ($(USE_HYPRE),TRUE)
    Pdirs += Extern/HYPRE
endif
ifeq ($(USE_PETSC), TRUE)
    Pdirs += Extern/PETSc
endif

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

