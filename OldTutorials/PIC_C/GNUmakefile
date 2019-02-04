AMREX_HOME ?= ../../

DEBUG	= FALSE

DIM	= 3

COMP    = gcc

USE_PARTICLES = TRUE

PRECISION = DOUBLE

USE_MPI   = TRUE
USE_OMP   = FALSE

###################################################

EBASE     = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

Pdirs 	:= Base Boundary AmrCore LinearSolvers/C_CellMG LinearSolvers/MLMG

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
