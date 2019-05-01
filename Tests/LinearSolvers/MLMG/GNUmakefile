AMREX_HOME ?= ../../../

DEBUG	?= FALSE
DIM	?= 3
COMP    ?= gnu

USE_MPI   ?= TRUE
USE_OMP   ?= FALSE

USE_HYPRE = FALSE

TINY_PROFILE ?= TRUE
PROFILE ?= FALSE
COMM_PROFILE ?= FALSE
TRACE_PROFILE ?= FALSE

DEFINES += -DAMREX_SOFT_PERF_COUNTERS

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/MLMG/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
