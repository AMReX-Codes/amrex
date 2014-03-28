# BOXLIB_HOME defines the directory in which we will find the BoxLib code
# If you have BOXLIB_HOME set as an environment variable then that definition
#    will override this one.
BOXLIB_HOME ?= /path/to/boxlib

# These are additional external (to BoxLib) libraries that you will need.
TRILINOS_DIR	 ?= /path/to/trilinos-11.0.3-Source
PARMETIS_LIB_DIR ?= /path/to/ParMetis-3.2.0
METIS_LIB_DIR	 ?= /path/to/ParmMetis-3.2.0
BLAS_LIB_DIR	 ?= /path/to/BLAS
LAPACK_LIB_DIR	 ?= /path/to/lapack-3.4.2/lib

PRECISION = DOUBLE

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 2
DIM	= 3

COMP =gcc
FCOMP=gfortran

# USE_MPI must be TRUE to run in parallel 
USE_MPI=TRUE
USE_OMP=FALSE

CPPFLAGS += -DMPICH_SKIP_MPICXX 

EBASE = driver

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

CEXE_sources += $(EBASE).cpp

include Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib

vpathdir += $(BOXLIB_HOME)/Src/C_BaseLib
vpathdir += $(BOXLIB_HOME)/Src/C_AMRLib

# trilinos library
INCLUDE_LOCATIONS += $(TRILINOS_DIR)/include
vpathdir          += $(TRILINOS_DIR)/include
LIBRARY_LOCATIONS += $(TRILINOS_DIR)/lib
LIBRARY_LOCATIONS  += $(PARMETIS_LIB_DIR)
LIBRARY_LOCATIONS  += $(METIS_LIB_DIR)
LIBRARY_LOCATIONS  += $(BLAS_LIB_DIR)
LIBRARY_LOCATIONS  += $(LAPACK_LIB_DIR)

LIBRARIES += -lbelosepetra -lbelos -lml -lgaleri -lamesos -lifpack -lepetraext -lepetra -lteuchos -ltriutils -laztecoo -lisorropia -lzoltan -lparmetis -lmetis -lblas -llapacke -llapack  -lgfortran

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)
	@echo SUCCESS

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
