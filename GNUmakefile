BOXLIB_HOME = ../BoxLib
PICSAR_HOME = ../picsar

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP    = gcc
FCOMP   = gfortran

USE_PARTICLES = TRUE

PRECISION = DOUBLE
USE_MPI   = FALSE

USE_OMP   = FALSE
EBASE     = main

include ./Make.package
include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package
include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package

include $(BOXLIB_HOME)/Src/LinearSolvers/C_to_F_MG/Make.package
include $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(BOXLIB_HOME)/Src/LinearSolvers/F_MG/FParallelMG.mak

include $(BOXLIB_HOME)/Src/F_BaseLib/FParallelMG.mak

VPATH_LOCATIONS   += $(PICSAR_HOME)/src
INCLUDE_LOCATIONS += $(PICSAR_HOME)/src

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
