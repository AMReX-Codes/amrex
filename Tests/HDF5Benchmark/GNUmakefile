AMREX_HOME = ../../

# DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP    = gnu

PRECISION = DOUBLE

USE_MPI   = TRUE
MPI_THREAD_MULTIPLE = FALSE

USE_OMP   = FALSE

TINY_PROFILE = TRUE

USE_PARTICLES = TRUE

###################################################

EBASE     = main

USE_HDF5  = TRUE
# HDF5_HOME = /gpfs/alpine/csc300/world-shared/gnu_build/hdf5/develop_build/hdf5
HDF5_HOME = $(OLCF_HDF5_ROOT)

# To use HDF5 with ZFP compression, follow the instructions at https://h5z-zfp.readthedocs.io to compile ZFP and H5Z_ZFP plugin, or use the following gcc build on Summit
USE_HDF5_ZFP = FALSE
H5Z_HOME  = /gpfs/alpine/csc300/world-shared/gnu_build/H5Z-ZFP/install
ZFP_HOME  = /gpfs/alpine/csc300/world-shared/gnu_build/zfp

# To use HDF5 with SZ compression, follow the instructions at https://github.com/szcompressor/SZ/tree/master/hdf5-filter/H5Z-SZ to compile SZ and H5Z_SZ plugin, or use the following gcc build on Summit
USE_HDF5_SZ = FALSE
SZ_HOME  = /gpfs/alpine/csc300/world-shared/gnu_build/SZ/install

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

# To use HDF5 asynchronous I/O VOL connector, follow the instructions at https://github.com/hpc-io/vol-async, add -DENABLE_WRITE_MEMCPY=1 to CFLAGS in src/Makefile, or use the following gcc build on Summit
USE_HDF5_ASYNC = FALSE
ABT_HOME       = /gpfs/alpine/csc300/world-shared/gnu_build/argobots/install
ASYNC_HOME     = /gpfs/alpine/csc300/world-shared/gnu_build/vol-async-memcpy/src

ifeq ($(USE_HDF5_ASYNC), TRUE)
DEFINES   	  += -DAMREX_USE_HDF5_ASYNC 
INCLUDE_LOCATIONS += $(ABT_HOME)/include $(ASYNC_HOME)
LIBRARIES 	  += -L$(ABT_HOME)/lib -L$(ASYNC_HOME) -lh5async -lasynchdf5 -labt -Wl,-rpath=$(ABT_HOME)/lib  -Wl,-rpath=$(ASYNC_HOME)
endif
