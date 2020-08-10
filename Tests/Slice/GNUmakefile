AMREX_HOME = ../../

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 3

COMP    = gnu

PRECISION = DOUBLE

USE_MPI   = TRUE
USE_OMP   = TRUE
USE_CUDA   = FALSE

TINY_PROFILE = TRUE

###################################################

EBASE     = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

HDF5_HOME   = /global/homes/s/sbyna/software/hdf5develop/install
USE_HDF5 = FALSE
