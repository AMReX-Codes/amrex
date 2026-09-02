DEBUG = FALSE

USE_MPI   = TRUE
USE_OMP   = FALSE
USE_HYPRE = FALSE
USE_EB    = TRUE

COMP = gnu

DIM = 3

AMREX_HOME = ../../..

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

Pdirs := Base Boundary LinearSolvers

ifeq ($(USE_EB),TRUE)
  USERSuffix += .EB
  Pdirs += AmrCore EB
endif

Ppack += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
