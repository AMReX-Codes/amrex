# BOXLIB_HOME defines the directory in which we will find all the BoxLib code -- if this is already defined
#    then the line below will be ignored
BOXLIB_HOME = ../..

# If you try to build with USE_MPI = TRUE and the make system can not find your mpi files, it may be because
#    MPI_HOME is not already defined on your machine.  In that case, you may need to set
#    the following variables according to the path of MPI on your particular machine.  
#    These lines are included here as example -- please uncomment and modify as appropriate.
#    We also recommend the CPPFLAG below.
#    MPI_HOME=/usr/local/mpich2
#    LIBRARY_LOCATIONS += $(MPI_HOME)/lib
#    INCLUDE_LOCATIONS += $(MPI_HOME)/include
#    CPPFLAGS += -DMPICH_SKIP_MPICXX

DEBUG	     = FALSE
USE_MPI      = TRUE
USE_OMP      = FALSE
USE_IPM      = FALSE
PROFILE      = FALSE
COMM_PROFILE = FALSE
USE_DUMPI    = FALSE
COMP         = gcc
FCOMP        = gfortran
#COMP         = Intel
#FCOMP        = Intel
DIM	     = 3
PRECISION    = DOUBLE
EBASE        = main

# Make BoxLib_C bottom CG solver use the old convergence criteria for comparison with BoxLib_F
CPPFLAGS += -DCG_USE_OLD_CONVERGENCE_CRITERIA


include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

CEXE_sources += $(EBASE).cpp

include Make.package
include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
include $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

all: $(executable)
	@echo SUCCESS

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
