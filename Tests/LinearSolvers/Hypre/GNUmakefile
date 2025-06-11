DEBUG = FALSE
TEST = FALSE
USE_ASSERTION = FALSE

USE_EB = TRUE

USE_MPI  = TRUE
USE_OMP  = FALSE

USE_HYPRE = TRUE

COMP = gnu

DIM = 2

BL_NO_FORT = TRUE

AMREX_HOME = ../../..

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

Pdirs := Base Boundary AmrCore EB
Pdirs += LinearSolvers/MLMG
Pdirs += Extern/HYPRE

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

