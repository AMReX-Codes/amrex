#
# $Id: GNUmakefile,v 1.2 2006-05-08 22:07:45 lijewski Exp $
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
# Base name of the executable.
#
EBASE = crsGrids
CEXE_sources += $(EBASE).cpp

DEFINES += -DBL_PARALLEL_IO
include $(TOP)/mk/Make.defs ./Make.package

include $(TOP)/BoxLib/Make.package

INCLUDE_LOCATIONS += . 
INCLUDE_LOCATIONS += $(TOP)/BoxLib 

vpath %.H   . $(TOP)/BoxLib
vpath %.cpp . $(TOP)/BoxLib
vpath %.F   . $(TOP)/BoxLib
vpath %.f   . $(TOP)/BoxLib

all: $(executable)

$(executable): $(LIBRARIES)

include $(TOP)/mk/Make.rules
