#
# $Id: GNUmakefile,v 1.5 2000-03-29 23:05:54 sstanley Exp $
#
PBOXLIB_HOME = ../../

TOP = $(PBOXLIB_HOME)
#
# Variables for the user to set ...
#
PRECISION     = DOUBLE
DEBUG	      = TRUE
PROFILE       = FALSE
DIM	      = 3
COMP          = KCC
USE_MPI       = TRUE

#
# Control Make Procedure
#
BUILD_INPLACE:=FALSE
BUILD_INPLACE:=TRUE

#
# Base name of the executable.
#
EBASE = viewMF

include $(TOP)/mk/Make.defs ./Make.package

CEXE_sources += $(EBASE).cpp

DEFINES += -DBL_PARALLEL_IO

#
# Local INCLUDES
#
ifeq ($(BUILD_INPLACE), TRUE)
  include $(TOP)/pBoxLib_2/Make.package
  INCLUDE_LOCATIONS += . $(TOP)/pBoxLib_2 $(TOP)/pAmrvis ../
else
  INCLUDE_LOCATIONS += . $(TOP)/include $(TOP)/pAmrvis ../
  LIBRARY_LOCATIONS += $(TOP)/lib/$(machineSuffix)
  LIBRARIES         += -lbox$(DIM)d
endif

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

ifeq ($(MACHINE), T3E)
ifeq ($(WHICHT3E), NERSC)
INCLUDE_LOCATIONS += /opt/ctl/cvt/3.1.0.0/include/X11
INCLUDE_LOCATIONS += /opt/ctl/cvt/3.1.0.0/include/Xm
INCLUDE_LOCATIONS += /opt/ctl/cvt/3.1.0.0/include/X11/Xaw
LIBRARY_LOCATIONS += /opt/ctl/cvt/ctv/lib
LIBRARIES += -lSM -lICE -lXt -lX11
DEFINES += -D_CRAY
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

ifeq ($(BUILD_INPLACE), TRUE)
  vpath %.F   .. $(TOP)/pAmrvis $(TOP)/pBoxLib_2
  vpath %.H   .. $(TOP)/pAmrvis $(TOP)/pBoxLib_2
  vpath %.cpp .. $(TOP)/pAmrvis $(TOP)/pBoxLib_2
else
  vpath %.F   .. $(TOP)/pAmrvis
  vpath %.H   .. $(TOP)/pAmrvis
  vpath %.cpp .. $(TOP)/pAmrvis
  vpath %.a   $(LIBRARY_LOCATIONS)
endif

all: $(executable)

ifeq ($(MACHINE),T3E)
#
# We don't explicity set the path of the mpi directory on the T3E.
#
$(executable): $(filter-out -lmpi, $(LIBRARIES))
else
$(executable): $(LIBRARIES)
endif
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

include $(TOP)/mk/Make.rules
