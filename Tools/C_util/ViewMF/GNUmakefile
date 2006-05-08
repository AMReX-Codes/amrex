#
# $Id: GNUmakefile,v 1.7 2006-05-08 21:59:46 lijewski Exp $
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
COMP          = g++
USE_MPI       = FALSE
#
# Base name of the executable.
#
EBASE = viewMFdiff

include $(TOP)/mk/Make.defs ./Make.package

CEXE_sources += $(EBASE).cpp

DEFINES += -DBL_NOLINEVALUES

include $(TOP)/BoxLib/Make.package
INCLUDE_LOCATIONS += . $(TOP)/BoxLib $(TOP)/pAmrvis ../

vpath %.F   .. $(TOP)/pAmrvis $(TOP)/BoxLib
vpath %.f   .. $(TOP)/pAmrvis $(TOP)/BoxLib
vpath %.H   .. $(TOP)/pAmrvis $(TOP)/BoxLib
vpath %.cpp .. $(TOP)/pAmrvis $(TOP)/BoxLib

all: $(executable)

include $(TOP)/mk/Make.rules
