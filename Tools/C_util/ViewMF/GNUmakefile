BOXLIB_HOME = ../..

TOP = $(BOXLIB_HOME)
#
# Variables for the user to set ...
#
PRECISION     = DOUBLE
DEBUG	      = TRUE
DIM	      = 3
COMP          = g++
USE_MPI       = FALSE
#
# Base name of the executable.
#
EBASE = viewMF
EBASE = mfMinMax

include $(TOP)/mk/Make.defs ./Make.package

CEXE_sources += $(EBASE).cpp

DEFINES += -DBL_NOLINEVALUES

include $(TOP)/BoxLib/Make.package
INCLUDE_LOCATIONS += . $(TOP)/BoxLib $(TOP)/pAmrvis ..

vpath %.F   .. $(TOP)/pAmrvis $(TOP)/BoxLib
vpath %.f   .. $(TOP)/pAmrvis $(TOP)/BoxLib
vpath %.H   .. $(TOP)/pAmrvis $(TOP)/BoxLib
vpath %.cpp .. $(TOP)/pAmrvis $(TOP)/BoxLib

all: $(executable)

include $(TOP)/mk/Make.rules
