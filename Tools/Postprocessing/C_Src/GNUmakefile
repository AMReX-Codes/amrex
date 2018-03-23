AMREX_HOME ?= ../../..

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP    = gcc

PRECISION = DOUBLE

USE_MPI   = TRUE
USE_OMP   = FALSE

###################################################

EBASE     = particle_compare

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

executable = particle_compare.exe

include ./Make.package
include $(AMREX_HOME)/Tools/GNUMake/Make.rules
