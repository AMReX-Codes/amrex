AMREX_HOME ?= ../../..

DEBUG = FALSE

DIM = 2

COMP = gnu

USE_MPI = TRUE
USE_OMP = FALSE

BL_NO_FORT = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

Pdirs := Base Boundary AmrCore LinearSolvers
Ppack += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)
include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
