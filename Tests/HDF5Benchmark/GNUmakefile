AMREX_HOME = ../../

# DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP    = gnu

PRECISION = DOUBLE

USE_MPI   = TRUE
MPI_THREAD_MULTIPLE = TRUE

USE_OMP   = FALSE

TINY_PROFILE = TRUE

###################################################

EBASE     = main

USE_HDF5  = TRUE
HDF5_HOME = /gpfs/alpine/csc300/world-shared/gnu_build/hpc-io/hdf5/async_optional_build/hdf5

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

# To use HDF5 asynchronous I/O VOL connector, follow the instructions at https://github.com/hpc-io/vol-async, use dset_memcpy branch
USE_HDF5_ASYNC = FALSE
ABT_HOME       = /gpfs/alpine/csc300/world-shared/gnu_build/hpc-io/vol-async/argobots/build
ASYNC_HOME     = /gpfs/alpine/csc300/world-shared/gnu_build/hpc-io/vol-async-memcpy/src

ifeq ($(USE_HDF5_ASYNC), TRUE)
DEFINES   	  += -DAMREX_USE_HDF5_ASYNC 
INCLUDE_LOCATIONS += $(ABT_HOME)/include $(ASYNC_HOME)
LIBRARIES 	  += -L$(ABT_HOME)/lib -L$(ASYNC_HOME) -lh5async -lasynchdf5 -labt -Wl,-rpath=$(ABT_HOME)/lib  -Wl,-rpath=$(ASYNC_HOME)
endif
