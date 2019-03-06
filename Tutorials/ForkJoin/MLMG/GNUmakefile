AMREX_HOME ?= ../../..

DEBUG ?= TRUE
DIM ?= 2
COMP ?= gnu
USE_MPI ?= TRUE
USE_OMP ?= FALSE

EBASE = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include Make.package
Pdirs := Base Boundary LinearSolvers/C_CellMG LinearSolvers/MLMG
Ppack += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)
include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
