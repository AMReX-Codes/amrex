AMREX_HOME = ../../

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP    = gnu

PRECISION = DOUBLE

USE_MPI   = TRUE
USE_OMP   = FALSE

TINY_PROFILE = TRUE

###################################################

EBASE     = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

USE_HDF5  = TRUE
HDF5_HOME = $(HOME)/hdf5/develop_build/hdf5

# To use HDF5 asynchronous I/O VOL connector, follow the instructions at https://bitbucket.hdfgroup.org/projects/HDF5VOL/repos/async/browse
USE_HDF5_ASYNC = FALSE
ABT_HOME       = $(HOME)/cori/argobots/install
ASYNC_HOME     = $(HOME)/hdf5vol/async/src
ifeq ($(USE_HDF5_ASYNC), TRUE)
DEFINES   	  += -DAMREX_USE_HDF5_ASYNC -DAMREX_MPI_THREAD_MULTIPLE
INCLUDE_LOCATIONS += $(ABT_HOME)/include $(ASYNC_HOME)
LIBRARIES 	  += -L$(ABT_HOME)/lib -L$(ASYNC_HOME) -lh5async -labt
endif
