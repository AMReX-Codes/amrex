AMREX_HOME ?= ..

DEBUG = FALSE
DIM = 3
COMP = gcc
PRECISION = DOUBLE
USE_MPI = FALSE
USE_OMP = FALSE
USE_CUDA = FALSE

DEBUG=TRUE

programs ?= fboxinfo fvarnames

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

multiple_executables = $(addsuffix .$(machineSuffix).ex, $(programs))
default: $(multiple_executables)

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Extern/amrdata/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

