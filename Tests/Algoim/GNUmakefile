DEBUG = FALSE
TEST = TRUE
USE_ASSERTION = TRUE

USE_EB = TRUE

USE_MPI  = FALSE
USE_OMP  = FALSE

USE_CUDA = FALSE

USE_ALGOIM = TRUE

COMP = gnu

DIM = 3

AMREX_HOME ?= ../../..

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include ./Make.package

Pdirs := Base Boundary AmrCore EB

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

