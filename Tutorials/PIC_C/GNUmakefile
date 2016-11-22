BOXLIB_HOME ?= ../../

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP    = gcc

USE_PARTICLES = TRUE

PRECISION = DOUBLE

USE_MPI   = TRUE
USE_OMP   = FALSE

###################################################

EBASE     = main

# location of finite-volume HPGMG on edison
HPGMG_DIR ?= ../../../hpgmg/finite-volume

USE_HPGMG = FALSE
HPGMG_FCYCLES = TRUE
HPGMG_POST_F_CYCLE_TYPE = F
HPGMG_HELMHOLTZ = FALSE
HPGMG_STENCIL_VARIABLE_COEFFICIENT = FALSE
HPGMG_USE_SUBCOMM = TRUE
HPGMG_BOTTOM_SOLVER= BICGSTAB
HPGMG_SMOOTHER = GSRB

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

include ./Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package
include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
include $(BOXLIB_HOME)/Src/C_ParticleLib/Make.package
include $(BOXLIB_HOME)/Src/LinearSolvers/C_to_F_MG/Make.package
include $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(BOXLIB_HOME)/Src/LinearSolvers/F_MG/FParallelMG.mak
include $(BOXLIB_HOME)/Src/F_BaseLib/FParallelMG.mak

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
