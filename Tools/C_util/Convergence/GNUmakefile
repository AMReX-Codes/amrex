#
# $Id: GNUmakefile,v 1.4 1999-06-04 01:03:25 sstanley Exp $
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
USE_MPI       = FALSE

BUILD_INPLACE:=FALSE
BUILD_INPLACE:=TRUE
#
# Use hgproj-serial -- only for testing.
#
# Touch Projection.cpp if you want to change the state of USE_HGPROJ_SERIAL.
#
USE_HGPROJ_SERIAL = TRUE
#
# The only stencils in which we're interested when using USE_HGPROJ_SERIAL.
#
ifeq ($(USE_HGPROJ_SERIAL),TRUE)
DEFINES += -DBL_USE_HGPROJ_SERIAL
ifeq ($(DIM),2)
PRVERSION = v5
else
PRVERSION = v7
endif
endif
#
# Base name of the executable.
#
EBASE = DiffSameGrid

include $(TOP)/mk/Make.defs ./Make.package

CEXE_sources += $(EBASE).cpp

DEFINES += -DBL_PARALLEL_IO


#
# Local INCLUDES
#
INCLUDE_LOCATIONS += ../../pAmrvis
vpath %.H   .. ../../pAmrvis
vpath %.cpp .. ../../pAmrvis
vpath %.F   .. ../../pAmrvis



ifeq ($(BUILD_INPLACE), TRUE)
  include $(TOP)/pBoxLib_2/Make.package
  include $(TOP)/bndrylib/Make.package
  INCLUDE_LOCATIONS += . ..
  INCLUDE_LOCATIONS += $(TOP)/iamrlib
  INCLUDE_LOCATIONS += $(TOP)/mglib
  INCLUDE_LOCATIONS += $(TOP)/tensorMG
  INCLUDE_LOCATIONS += $(TOP)/amrlib
  INCLUDE_LOCATIONS += $(TOP)/bndrylib
  INCLUDE_LOCATIONS += $(TOP)/hgproj
  INCLUDE_LOCATIONS += $(TOP)/pBoxLib_2
else
  ifeq ($(USE_HGPROJ_SERIAL),FALSE)
    INCLUDE_LOCATIONS += . .. $(TOP)/include
    LIBRARY_LOCATIONS += $(TOP)/lib/$(machineSuffix)
    LIBRARIES         += -liamr$(DIM)d -lmg$(DIM)d -lamr$(DIM)d -lbndry$(DIM)d -lproj$(DIM)d -lbox$(DIM)d
  else
    INCLUDE_LOCATIONS += . .. $(TOP)/hgproj-serial
    INCLUDE_LOCATIONS += $(TOP)/hgproj-serial/include/$(DIM)d.$(PRVERSION)
    INCLUDE_LOCATIONS += $(TOP)/include
    LIBRARY_LOCATIONS += $(TOP)/hgproj-serial/lib/$(machineSuffix) $(TOP)/lib/$(machineSuffix)
    LIBRARIES         += -liamr$(DIM)d -lmg$(DIM)d -lamr$(DIM)d -lbndry$(DIM)d -lproj$(DIM)d.$(PRVERSION) -lbox$(DIM)d
  endif
endif

ifeq ($(USE_ARRAYVIEW),TRUE)
DEFINES += -DBL_USE_ARRAYVIEW
DEFINES += -DBL_ARRAYVIEW_TAGBOX
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
  vpath %.cpp . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib   \
                $(TOP)/bndrylib $(TOP)/hgproj $(TOP)/pBoxLib_2
  vpath %.H   . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib   \
                $(TOP)/bndrylib $(TOP)/hgproj $(TOP)/pBoxLib_2
  vpath %.h   . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib   \
                $(TOP)/bndrylib $(TOP)/hgproj $(TOP)/pBoxLib_2
  vpath %.F   . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib   \
                $(TOP)/bndrylib $(TOP)/hgproj $(TOP)/pBoxLib_2
else
  vpath %.H   ..
  vpath %.cpp ..
  vpath %.a   $(LIBRARY_LOCATIONS)
endif

all: $(executable)

ifeq ($(MACHINE),T3E)
ifeq ($(WHICHT3E), NERSC)
DEFINES += -D_CRAY
endif

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
	cd $(TOP)/bndrylib;  $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install
#
# Cleanup libraries.
#
cleanlibs:
	cd $(TOP)/pBoxLib_2; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) clean
	cd $(TOP)/bndrylib;  $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) clean
endif

include $(TOP)/mk/Make.rules
