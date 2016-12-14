AMREX_HOME ?= ../../

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 3

COMP    = gcc

USE_PARTICLES = TRUE

PRECISION = DOUBLE
USE_MPI   = TRUE

USE_OMP   = FALSE

include $(AMREX_HOME)/Tools/C_mk/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/C_BaseLib/Make.package
include $(AMREX_HOME)/Src/C_ParticleLib/Make.package
include $(AMREX_HOME)/Src/C_BoundaryLib/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_to_F_MG/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/F_MG/FParallelMG.mak
include $(AMREX_HOME)/Src/F_BaseLib/FParallelMG.mak

include $(AMREX_HOME)/Tools/C_mk/Make.rules
