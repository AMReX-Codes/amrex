BOXLIB_HOME ?= ../BoxLib
PICSAR_HOME = ../picsar

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 3

COMP    = gcc
FCOMP   = gfortran

USE_PARTICLES = TRUE

PRECISION = DOUBLE
USE_MPI   = TRUE

USE_OMP   = FALSE
EBASE     = main

include ./Make.package
include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package
include $(BOXLIB_HOME)/Src/C_ParticleLib/Make.package
include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package

#include $(BOXLIB_HOME)/Src/F_BaseLib/FParallelMG.mak

include $(PICSAR_HOME)/src/Make.package

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
