#
# $Id: GNUmakefile,v 1.3 1999-05-24 18:11:46 car Exp $
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
USE_MPI       = FALSE

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
INCLUDE_LOCATIONS += . $(TOP)/include $(TOP)/pAmrvis ../
LIBRARY_LOCATIONS += $(TOP)/lib/$(machineSuffix)
LIBRARIES         += -lbox$(DIM)d

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

vpath %.F   .. $(TOP)/pAmrvis
vpath %.H   .. $(TOP)/pAmrvis
vpath %.cpp .. $(TOP)/pAmrvis
vpath %.a   $(LIBRARY_LOCATIONS)

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
