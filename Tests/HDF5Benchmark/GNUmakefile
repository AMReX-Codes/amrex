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
HDF5_HOME = /lustre/orion/geo130/world-shared/gnu/hdf5.1.14.1.sync.cmake
# HDF5_HOME = $(HDF5_ROOT)
# HDF5_HOME = /lustre/orion/geo130/world-shared/gnu/hdf5.1.14.1.async.cmake

# To use HDF5 with ZFP compression, follow the instructions at https://h5z-zfp.readthedocs.io to compile ZFP and H5Z_ZFP plugin, or use the following gcc build on Frontier
USE_HDF5_ZFP = FALSE
H5Z_HOME  = /lustre/orion/geo130/world-shared/gnu/H5Z-ZFP
ZFP_HOME  = /lustre/orion/geo130/world-shared/gnu/zfp

# To use HDF5 with SZ compression, follow the instructions at https://github.com/szcompressor/SZ/tree/master/hdf5-filter/H5Z-SZ to compile SZ and H5Z_SZ plugin, or use the following gcc build on Frontier
USE_HDF5_SZ = FALSE
SZ_HOME  = /lustre/orion/geo130/world-shared/gnu/SZ

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

# To use HDF5 asynchronous I/O VOL connector, follow the instructions at https://github.com/hpc-io/vol-async, add -DENABLE_WRITE_MEMCPY=1 to CFLAGS in src/Makefile, or use the following gcc build on Frontier
USE_HDF5_ASYNC = FALSE
ABT_HOME       = /lustre/orion/geo130/world-shared/gnu/argobots
ASYNC_HOME     = /lustre/orion/geo130/world-shared/gnu/vol-async-memcpy

ifeq ($(USE_HDF5_ASYNC), TRUE)
DEFINES   	  += -DAMREX_USE_HDF5_ASYNC 
LIBRARIES 	  += -L$(ABT_HOME)/lib -L$(ASYNC_HOME)/lib -lh5async -labt -Wl,-rpath=$(ABT_HOME)/lib  -Wl,-rpath=$(ASYNC_HOME)/lib
endif
