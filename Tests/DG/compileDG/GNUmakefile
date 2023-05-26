AMREX_HOME ?= ../../..

DEBUG = TRUE

USE_MPI = TRUE
USE_OMP = FALSE

DIM = 2

EBASE := main

# Set pre-processor directives
# This is also where, e.g., USE_HDF5 and USE_PETSC are used
include $(AMREX_HOME)/Tools/GNUMake/Make.defs

# Include basic AMReX framework
include $(AMREX_HOME)/Src/Base/Make.package

ifneq ("$(wildcard $(AMREX_HOME)/Src/DG/Make.package)","")
  # Requires amrex fork from github.com:dunhamsj/amrex
  # Branch: MeshRefinement_DG
  include $(AMREX_HOME)/Src/DG/Make.package
endif

include Make.package

# Include rules for compilation
include $(AMREX_HOME)/Tools/GNUMake/Make.rules
