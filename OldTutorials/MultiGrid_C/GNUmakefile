# AMREX_HOME defines the directory in which we will find all the BoxLib code -- if this is already defined
#    then the line below will be ignored
AMREX_HOME ?= ../..

# HYPRE_DIR defines the directory in which you can find hypre -- if this is already defined 
#    then the line below will be ignored
HYPRE_DIR ?= $(HOME)/lib/hypre

# location of finite-volume HPGMG, if you decide to use it
HPGMG_DIR ?= $(HOME)/hpgmg/finite-volume

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
DEBUG	= FALSE

DIM	= 3
DIM	= 2

COMP =gcc

USE_HYPRE = FALSE

USE_HPGMG = TRUE
USE_HPGMG = FALSE
HPGMG_FCYCLES = TRUE
HPGMG_POST_F_CYCLE_TYPE = V
HPGMG_HELMHOLTZ = TRUE
HPGMG_STENCIL_VARIABLE_COEFFICIENT = FALSE
HPGMG_USE_SUBCOMM = TRUE
HPGMG_BOTTOM_SOLVER= BICGSTAB
HPGMG_SMOOTHER = GSRB

# Make BoxLib_C bottom CG solver use the old convergence criteria for comparison with BoxLib_F
CPPFLAGS += -DCG_USE_OLD_CONVERGENCE_CRITERIA

# USE_MPI must be TRUE if USE_HYPRE
USE_MPI=TRUE
USE_OMP=TRUE
USE_OMP=FALSE

EBASE = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += $(EBASE).cpp

include Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG4/Make.package
include $(AMREX_HOME)/Src/Base/Make.package

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Boundary
vpathdir += $(AMREX_HOME)/Src/Boundary

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
vpathdir += $(AMREX_HOME)/Src/Base

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/LinearSolvers/C_CellMG
vpathdir += $(AMREX_HOME)/Src/LinearSolvers/C_CellMG

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/LinearSolvers/C_CellMG4
vpathdir += $(AMREX_HOME)/Src/LinearSolvers/C_CellMG4

ifeq ($(USE_HYPRE), TRUE)
  # hypre driver
  include HypreABecLap/Make.package
  INCLUDE_LOCATIONS += HypreABecLap
  vpathdir          += HypreABecLap

  # hypre library
  INCLUDE_LOCATIONS += $(HYPRE_DIR)/include
  vpathdir          += $(HYPRE_DIR)/include
  LIBRARY_LOCATIONS += $(HYPRE_DIR)/lib
  LIBRARIES += -lHYPRE

  DEFINES += -DUSEHYPRE
endif

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)
vpath %.F90 : . $(vpathdir)

# The HPGMG hooks in Make.defs populate VPATH_LOCATIONS, not vpathdir, so we
# need to grab these too. In fact everything in Make.defs uses VPATH_LOCATIONS
# instead of vpathdir. Perhaps they should be consolidated ...
vpath %.c   : . $(VPATH_LOCATIONS)
vpath %.h   : . $(VPATH_LOCATIONS)
vpath %.cpp : . $(VPATH_LOCATIONS)
vpath %.H   : . $(VPATH_LOCATIONS)
vpath %.F   : . $(VPATH_LOCATIONS)
vpath %.f   : . $(VPATH_LOCATIONS)
vpath %.f90 : . $(VPATH_LOCATIONS)
vpath %.F90 : . $(VPATH_LOCATIONS)

all: $(executable)
	@echo SUCCESS

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
