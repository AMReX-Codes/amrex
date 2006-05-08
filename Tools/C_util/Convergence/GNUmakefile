#
# $Id: GNUmakefile,v 1.6 2006-05-08 21:41:33 lijewski Exp $
#
PBOXLIB_HOME = ../../

TOP = $(PBOXLIB_HOME)
#
# Variables for the user to set ...
#
PRECISION     = DOUBLE
DEBUG	      = TRUE
PROFILE       = FALSE
DIM	      = 2
COMP          = KCC
USE_MPI       = FALSE
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

DEFINES           += -DBL_NOLINEVALUES

CEXE_sources += $(EBASE).cpp
#
# Local INCLUDES
#
INCLUDE_LOCATIONS += ../../pAmrvis

vpath %.H   .. ../../pAmrvis
vpath %.cpp .. ../../pAmrvis
vpath %.F   .. ../../pAmrvis
vpath %.f   .. ../../pAmrvis

include $(TOP)/BoxLib/Make.package
include $(TOP)/bndrylib/Make.package

INCLUDE_LOCATIONS += . ..
INCLUDE_LOCATIONS += $(TOP)/iamrlib
INCLUDE_LOCATIONS += $(TOP)/mglib
INCLUDE_LOCATIONS += $(TOP)/tensorMG
INCLUDE_LOCATIONS += $(TOP)/amrlib
INCLUDE_LOCATIONS += $(TOP)/bndrylib
INCLUDE_LOCATIONS += $(TOP)/hgproj
INCLUDE_LOCATIONS += $(TOP)/BoxLib

ifeq ($(USE_ARRAYVIEW),TRUE)
DEFINES += -DBL_USE_ARRAYVIEW
DEFINES += -DBL_ARRAYVIEW_TAGBOX
endif

vpath %.cpp . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib   \
              $(TOP)/bndrylib $(TOP)/hgproj $(TOP)/BoxLib
vpath %.H   . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib   \
              $(TOP)/bndrylib $(TOP)/hgproj $(TOP)/BoxLib
vpath %.h   . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib   \
              $(TOP)/bndrylib $(TOP)/hgproj $(TOP)/BoxLib
vpath %.F   . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib   \
              $(TOP)/bndrylib $(TOP)/hgproj $(TOP)/BoxLib
vpath %.f   . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib   \
              $(TOP)/bndrylib $(TOP)/hgproj $(TOP)/BoxLib

all: $(executable)

$(executable): $(LIBRARIES)

include $(TOP)/mk/Make.rules
