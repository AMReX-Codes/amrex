
TOP = $(AMREX_HOME)
USE_EB = FALSE
#
# Variables for the user to set ...
#
PRECISION     = DOUBLE
DEBUG	      = FALSE
DIM	          = 3
COMP          = g++
USE_MPI       = TRUE 

#
# Base name of the executable.
#
EBASE = AugmentPlotfile
CEXE_sources += $(EBASE).cpp
FEXE_headers += $(EBASE)_F.H
f90EXE_sources += $(EBASE)_F.f90 


include $(AMREX_HOME)/Tools/GNUMake/Make.defs


include ../Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/Extern/amrdata/Make.package


INCLUDE_LOCATIONS += .
INCLUDE_LOCATIONS += ..
vpathdir += ..

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
vpathdir += $(AMREX_HOME)/Src/Base

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Boundary
vpathdir += $(AMREX_HOME)/Src/Boundary

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Amr
vpathdir += $(AMREX_HOME)/Src/Amr

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata
vpathdir += $(AMREX_HOME)/Src/Extern/amrdata

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/ProfParser
vpathdir += $(AMREX_HOME)/Src/Extern/ProfParser


vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

#include $(TOP)/Tools/GNUMake/Make.defs ./Make.package

#ifeq ($(USE_ARRAYVIEW),TRUE)
#  DEFINES += -DBL_USE_ARRAYVIEW
#  DEFINES += -DBL_ARRAYVIEW_TAGBOX
#endif





