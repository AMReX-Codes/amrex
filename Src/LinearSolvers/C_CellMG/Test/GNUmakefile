PRECISION = DOUBLE
PROFILE = FALSE
DEBUG	= FALSE
DEBUG	= TRUE
DIM	= 2
DIM	= 3
COMP = g++
COMP = KCC
USE_MPI=FALSE
USE_MPI=TRUE
NAMESPACE=FALSE

USE_ARRAYVIEW = TRUE
USE_ARRAYVIEW = FALSE

USE_HYPRE=FALSE
USE_HYPRE=TRUE

PBOXLIB_HOME := ../..

EBASE = mg
LBASE = 

include $(PBOXLIB_HOME)/mk/Make.defs

# multigrid test sources
CEXE_sources += main.cpp

ifeq ($(USE_HYPRE),TRUE)
  HYPRE_DIR := $(PBOXLIB_HOME)/hypre
  LIBRARY_LOCATIONS += $(HYPRE_DIR)/lib
  INCLUDE_LOCATIONS += $(HYPRE_DIR)/include
  INCLUDE_LOCATIONS += $(HYPRE_DIR)/mglib_interface
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
  include $(HYPRE_DIR)/mglib_interface/Make.package
endif

include $(PBOXLIB_HOME)/bndrylib/Make.package
include $(PBOXLIB_HOME)/mglib/Make.package
include $(PBOXLIB_HOME)/pBoxLib_2/Make.package
include $(PBOXLIB_HOME)/util/Make.package

INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/bndrylib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/pBoxLib_2
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/mglib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/util

# We need these from amrlib 

CEXE_sources += BCRec.cpp 
CEXE_headers += BCRec.H BC_Types.H
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/amrlib

# Needed from AmrVis
CEXE_headers += AmrData.H
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/pAmrvis


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
endif

CXXFLAGS +=
CXXOPTF +=
CXXDEBF +=

CFLAGS +=
COPTF +=
CDEBF +=

ifeq ($(USE_HYPRE), TRUE)
  vpath %.cpp : . $(HYPRE_DIR)/mglib_interface $(PBOXLIB_HOME)/pBoxLib_2 $(PBOXLIB_HOME)/amrlib $(PBOXLIB_HOME)/bndrylib $(PBOXLIB_HOME)/mglib $(PBOXLIB_HOME)/util
  vpath %.F   : . $(HYPRE_DIR)/mglib_interface $(PBOXLIB_HOME)/pBoxLib_2 $(PBOXLIB_HOME)/amrlib $(PBOXLIB_HOME)/bndrylib $(PBOXLIB_HOME)/mglib
else
  vpath %.cpp : . $(PBOXLIB_HOME)/pBoxLib_2 $(PBOXLIB_HOME)/amrlib $(PBOXLIB_HOME)/bndrylib $(PBOXLIB_HOME)/mglib $(PBOXLIB_HOME)/util
  vpath %.F   : . $(PBOXLIB_HOME)/pBoxLib_2 $(PBOXLIB_HOME)/amrlib $(PBOXLIB_HOME)/bndrylib $(PBOXLIB_HOME)/mglib
endif

all: $(executable)

include $(PBOXLIB_HOME)/mk/Make.rules
