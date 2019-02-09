# AMREX_HOME defines the directory in which we will find all the BoxLib code -- if this is already defined
#    then the line below will be ignored
AMREX_HOME ?= ../../..

# HYPRE_DIR defines the directory in which you can find hypre -- if this is already defined 
#    then the line below will be ignored
HYPRE_DIR ?= $(HOME)/lib/hypre

# If you try to build with USE_MPI = TRUE and the make system can not find your mpi files, it may be because
#    MPI_HOME is not already defined on your machine.  In that case, you may need to set
#    the following variables according to the path of MPI on your particular machine.  
#    These lines are included here as example -- please uncomment and modify as appropriate.
#    We also recommend the CPPFLAG below.
#    MPI_HOME=/usr/local/mpich2
#    LIBRARY_LOCATIONS += $(MPI_HOME)/lib
#    INCLUDE_LOCATIONS += $(MPI_HOME)/include
#    CPPFLAGS += -DMPICH_SKIP_MPICXX

PRECISION = DOUBLE

DEBUG	= TRUE
#DEBUG	= FALSE

DIM	= 3

COMP =gcc

USE_HYPRE = FALSE

# Make BoxLib_C bottom CG solver use the old convergence criteria for comparison with BoxLib_F
CPPFLAGS += -DCG_USE_OLD_CONVERGENCE_CRITERIA

USE_MPI=TRUE

EBASE = LST

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += main.cpp

include Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package

ifeq ($(USE_HYPRE), TRUE)
  # USE_MPI must be TRUE if USE_HYPRE
  USE_MPI=TRUE

  # hypre driver
  include HypreABecLap/Make.package
  INCLUDE_LOCATIONS += HypreABecLap
  VPATH_LOCATIONS   += HypreABecLap

  # hypre library
  INCLUDE_LOCATIONS += $(HYPRE_DIR)/include
  VPATH_LOCATIONS   += $(HYPRE_DIR)/include
  LIBRARY_LOCATIONS += $(HYPRE_DIR)/lib
  LIBRARIES += -lHYPRE

  DEFINES += -DUSEHYPRE
endif

all: $(executable)
	@echo SUCCESS

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
