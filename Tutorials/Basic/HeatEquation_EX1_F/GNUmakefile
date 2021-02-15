# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME = ../../..

DEBUG     = TRUE

USE_MPI   = TRUE
USE_OMP   = FALSE
USE_MPI3  = FALSE

USE_F_INTERFACES = TRUE

MEM_PROFILE  = FALSE
TINY_PROFILE = FALSE

COMP      = gnu

DIM       = 3

EBASE     = main

include ./Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/F_Interfaces/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
