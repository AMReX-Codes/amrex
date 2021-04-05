USE_MPI  = TRUE
USE_OMP  = FALSE

COMP = gnu

DIM = 3

DEBUG = FALSE

AMREX_HOME = ../../..

USE_EB = TRUE

USE_HYPRE = TRUE
USE_HYPRE = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

Pdirs := AmrCore 
Pdirs += Base 
Pdirs += Boundary 
Pdirs += EB
Pdirs += LinearSolvers/MLMG LinearSolvers/Projections

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
