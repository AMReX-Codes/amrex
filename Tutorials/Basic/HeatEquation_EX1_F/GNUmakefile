# AMREX_HOME defines the directory in which we will find all the AMReX code
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../..

DEBUG     = TRUE

USE_MPI   = TRUE
USE_OMP   = FALSE
USE_MPI3  = FALSE

MEM_PROFILE  = FALSE
TINY_PROFILE = FALSE

COMP      = g++

DIM       = 3

EBASE     = main

include ./Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/F_Interfaces/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
