#
# $Id: GNUmakefile,v 1.1 2000-03-23 22:09:44 sstanley Exp $
#
PBOXLIB_HOME = ../../

TOP = $(PBOXLIB_HOME)
#
# Variables for the user to set ...
#
PRECISION     = DOUBLE
DEBUG	      = FALSE
PROFILE       = FALSE
DIM	      = 3
COMP          = KCC
#
# Control Make Procedure
#
BUILD_INPLACE:=FALSE
BUILD_INPLACE:=TRUE


#
# Base name of the executable.
#
EBASE = crsGrids
CEXE_sources += $(EBASE).cpp

DEFINES += -DBL_PARALLEL_IO
include $(TOP)/mk/Make.defs ./Make.package

#
### Library and Include Search Paths ###
#######################################
ifeq ($(BUILD_INPLACE), TRUE)
  include $(TOP)/pBoxLib_2/Make.package

  INCLUDE_LOCATIONS += . 
  INCLUDE_LOCATIONS += $(TOP)/pBoxLib_2 
else
  INCLUDE_LOCATIONS += . 
  INCLUDE_LOCATIONS += $(TOP)/include
  LIBRARY_LOCATIONS += $(TOP)/lib/$(machineSuffix)
  LIBRARIES         += -lbox$(DIM)d
endif

#
#
#
ifeq ($(MACHINE),OSF1)
  #
  # Some additional stuff for our preferred development/debugging environment.
  #
  ifeq ($(PRECISION),DOUBLE)
    FFLAGS += -real_size 64
  endif
  FDEBF += -C
  FDEBF += -warn argument_checking
  FDEBF += -warn declarations
  ifneq ($(FC), f90)
    FDEBF += -warn truncated_source
    FDEBF += -warn unused
  endif
endif
#
# For Running 3rd Only
#
3RD = 1
3RD =
ifdef 3RD
LDFLAGS += --link_command_prefix 3rd
LDFLAGS += -non_shared -v
endif

#
# Set Search Paths
#
ifeq ($(BUILD_INPLACE), TRUE)
  vpath %.H   . $(TOP)/pBoxLib_2
  vpath %.cpp . $(TOP)/pBoxLib_2
  vpath %.F   . $(TOP)/pBoxLib_2
else
  vpath %.H   . $(TOP)/include
  vpath %.cpp . 
  vpath %.F   . 
  vpath %.a   $(LIBRARY_LOCATIONS)
endif


#
### Make Rules ###
##################

all: $(executable)

ifeq ($(MACHINE),T3E)
  #
  # We don't explicity set the path of the mpi directory on the T3E.
  #
  $(executable): $(filter-out -lmpi, $(LIBRARIES))
else
  $(executable): $(LIBRARIES)
endif

ifneq ($(BUILD_INPLACE), TRUE)
#
# Build and install all libraries needed in an appropriate order.
#
libs:
	cd $(TOP)/pBoxLib_2; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install

#
# Cleanup libraries.
#
cleanlibs:
	cd $(TOP)/pBoxLib_2; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) clean
endif

include $(TOP)/mk/Make.rules
