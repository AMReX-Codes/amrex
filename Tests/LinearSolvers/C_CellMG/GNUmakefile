PRECISION = DOUBLE
PROFILE = TRUE
PROFILE = FALSE
DEBUG	= FALSE
DEBUG	= TRUE
DIM	= 3
DIM	= 2
COMP = g++
COMP = KCC
USE_MPI=TRUE
USE_MPI=FALSE
USE_THREADS=TRUE
USE_THREADS=FALSE
NAMESPACE=FALSE
STRICTLY=TRUE

USE_ARRAYVIEW = TRUE
USE_ARRAYVIEW = FALSE

USE_HYPRE=TRUE
USE_HYPRE=FALSE

USE_BL3 = TRUE
USE_BL3 = FALSE

PBOXLIB_HOME := ../..

EBASE = main
LBASE = 

include $(PBOXLIB_HOME)/mk/Make.defs

# multigrid test sources
CEXE_sources += $(EBASE).cpp

ifeq ($(USE_HYPRE),TRUE)
  HYPRE_DIR := $(PBOXLIB_HOME)/hypre
  LIBRARY_LOCATIONS += $(HYPRE_DIR)/lib
  INCLUDE_LOCATIONS += $(HYPRE_DIR)/include
  INCLUDE_LOCATIONS += $(HYPRE_DIR)/HypreABec
  vpathdir += $(HYPRE_DIR)/HypreABec
  LIBRARIES +=  -lHYPRE_ls			\
	 -lHYPRE_mv				\
	 -lHYPRE_parcsr_ls			\
	 -lHYPRE_DistributedMatrixPilutSolver 	\
	 -lHYPRE_ParaSails 			\
	 -lHYPRE_MatrixMatrix 			\
	 -lHYPRE_DistributedMatrix 		\
	 -lHYPRE_IJ_mv				\
	 -lHYPRE_parcsr_mv			\
	 -lHYPRE_seq_mv				\
	 -lHYPRE_utilities
  CPPFLAGS += -DMG_USE_HYPRE
  include $(HYPRE_DIR)/HypreABec/Make.package
endif


include Make.package
include $(PBOXLIB_HOME)/bndrylib/Make.package
include $(PBOXLIB_HOME)/mglib/Make.package
include $(PBOXLIB_HOME)/BoxLib/Make.package
include $(PBOXLIB_HOME)/util/Make.package

INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/bndrylib
vpathdir += $(PBOXLIB_HOME)/bndrylib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/BoxLib
vpathdir += $(PBOXLIB_HOME)/BoxLib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/iamrlib
vpathdir += $(PBOXLIB_HOME)/iamrlib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/mglib
vpathdir += $(PBOXLIB_HOME)/mglib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/util
vpathdir += $(PBOXLIB_HOME)/util


# We need these from amrlib 

CEXE_sources += BCRec.cpp 
CEXE_headers += BCRec.H BC_Types.H
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/amrlib
vpathdir += $(PBOXLIB_HOME)/amrlib

# Needed from AmrVis
CEXE_headers += AmrData.H
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/pAmrvis
vpathdir += $(PBOXLIB_HOME)/pAmrvis

# BL3 Support for threads and profiling
ifeq ($(USE_BL3),TRUE)
  CPPFLAGS += -DBL3_PROFILING
  ifeq ($(USE_MPI),TRUE)
    CPPFLAGS += -DBL3_USE_MPI
    CEXE_sources += BL3Mpi.cpp 
  endif
  ifeq ($(USE_THREADS),TRUE)
    DEFINES += -DBL3_PTHREADS
  endif
  include $(PBOXLIB_HOME)/BoxLib3/Make.package.mini
  INCLUDE_LOCATIONS += $(PBOXLIB_HOME)
  vpathdir += $(PBOXLIB_HOME)/BoxLib3
  ifeq ($(COMP),g++)
    INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/BoxLib3/std
  endif
endif

ifeq ($(DEBUG),TRUE)
  ifeq ($(USE_ARRAYVIEW),TRUE)
    ARRAYVIEWDIR  = $(PBOXLIB_HOME)/util
    INCLUDE_LOCATIONS += $(ARRAYVIEWDIR)
    CEXE_headers += DatasetClient.H ArrayView.H
    CEXE_sources += DatasetClient.cpp
    DEFINES += -DBL_USE_ARRAYVIEW
    #DEFINES += -DBL_ARRAYVIEW_TAGBOX
  endif
  CEXE_headers += TV_TempWrite.H
endif

ifeq ($(MACHINE),OSF1)
  3RD=1
  3RD=
  ifdef 3RD
    # FOR RUNNING 3RD ONLY
    LDFLAGS += --link_command_prefix 3rd
    #CXXDEBF = +K0 --link_command_prefix 3rd -non_shared
    LDFLAGS += -non_shared -v
    #LIBRARIES += -ldnet_stub
    #FDEBF += -automatic
    # FOR RUNNING 3RD ONLY
  endif
  DEFINES += -DBL3_NANOSLEEP
  LIBRARIES+=-lrt
  FFLAGS += -automatic
endif

CXXFLAGS +=
CXXOPTF +=
CXXDEBF +=

CFLAGS +=
COPTF +=
CDEBF +=

vpath %.cpp : . $(vpathdir)
vpath %.F   : . $(vpathdir)

all: $(executable)

include $(PBOXLIB_HOME)/mk/Make.rules
