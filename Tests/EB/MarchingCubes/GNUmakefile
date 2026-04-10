DEBUG = FALSE
TEST = TRUE
USE_ASSERTION = TRUE

USE_EB = TRUE

USE_MPI  = FALSE
USE_OMP  = FALSE

COMP = gnu

DIM = 3

BL_NO_FORT = TRUE

AMREX_HOME = ../../..

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

Pdirs := Base Boundary AmrCore EB

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

