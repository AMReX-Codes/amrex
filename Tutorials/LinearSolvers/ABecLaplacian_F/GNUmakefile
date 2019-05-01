DEBUG = FALSE

USE_MPI  = TRUE
USE_OMP  = FALSE
USE_F_INTERFACES = TRUE

COMP = gnu

DIM = 3

AMREX_HOME ?= ../../..

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

Pdirs 	:= Base Boundary AmrCore LinearSolvers/C_CellMG LinearSolvers/MLMG
Pdirs   += F_Interfaces/Base F_Interfaces/AmrCore F_Interfaces/LinearSolvers

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

