DEBUG = FALSE
USE_MPI  = TRUE
USE_OMP  = FALSE
COMP = gnu
DIM = 3
BL_NO_FORT = TRUE

USE_CUDA  = FALSE
USE_SYCL  = FALSE
USE_HIP   = FALSE

TINY_PROFILE = FALSE

AMREX_HOME = ../../..

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

Pdirs 	:= Base Boundary LinearSolvers

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

