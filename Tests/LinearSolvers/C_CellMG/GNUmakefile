AMREX_HOME := ../../..

PRECISION = DOUBLE

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 3
DIM	= 2

COMP =g++
FCOMP=gfortran

USE_MPI=FALSE

EBASE = main

include $(AMREX_HOME)/Tools/C_mk/Make.defs

CEXE_sources += $(EBASE).cpp

include ./Make.package

include $(AMREX_HOME)/Src/C_BoundaryLib/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(AMREX_HOME)/Src/C_BaseLib/Make.package
include $(AMREX_HOME)/Src/C_AmrCoreLib/Make.package

include $(AMREX_HOME)/Src/LinearSolvers/C_to_F_MG/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/F_MG/FParallelMG.mak
include $(AMREX_HOME)/Src/F_BaseLib/FParallelMG.mak

include $(AMREX_HOME)/Tools/C_mk/Make.rules
