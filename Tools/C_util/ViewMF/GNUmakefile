AMREX_HOME = ../..

TOP = $(AMREX_HOME)
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

include $(TOP)/amrex/Make.package
INCLUDE_LOCATIONS += . $(TOP)/amrex $(TOP)/pAmrvis ..

vpath %.F   .. $(TOP)/pAmrvis $(TOP)/amrex
vpath %.f   .. $(TOP)/pAmrvis $(TOP)/amrex
vpath %.H   .. $(TOP)/pAmrvis $(TOP)/amrex
vpath %.cpp .. $(TOP)/pAmrvis $(TOP)/amrex

all: $(executable)

include $(TOP)/mk/Make.rules
