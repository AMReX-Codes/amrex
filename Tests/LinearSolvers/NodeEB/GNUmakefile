
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

AMREX_HOME ?= ../../..

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include ./Make.package

Pdirs := Base Boundary AmrCore
Pdirs += EB
Pdirs += LinearSolvers/C_CellMG LinearSolvers/MLMG

ifeq ($(DIM),3)
USE_ALGOIM = TRUE
else
USE_ALGOIM = FALSE
endif

ifeq ($(USE_ALGOIM),TRUE)
    ALGOIM_HOME ?= /path/to/algoim
    BLITZ_HOME  ?= /path/to/blitz
    INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/Algoim
    INCLUDE_LOCATIONS += $(ALGOIM_HOME)/src
    INCLUDE_LOCATIONS += $(BLITZ_HOME)/include
    VPATH_LOCATIONS   += $(AMREX_HOME)/Src/Extern/Algoim
    VPATH_LOCATIONS   += $(ALGOIM_HOME)/src
    DEFINES += -DUSE_ALGOIM
endif

ifeq ($(USE_ALGOIM),TRUE)
    Pdirs += Extern/Algoim
endif

ifeq ($(USE_HYPRE),TRUE)
    Pdirs += Extern/HYPRE
endif
ifeq ($(USE_PETSC), TRUE)
    Pdirs += Extern/PETSc
endif

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

