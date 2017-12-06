
TOP = $(AMREX_HOME)
USE_EB =FALSE
#
# Variables for the user to set ...
#
PRECISION     = DOUBLE
DEBUG	      = TRUE
DEBUG	      = FALSE
#DIM	      = 3
DIM	      = 2
COMP          = g++
#USE_MPI       = TRUE 
USE_MPI       = FALSE
#
# Base name of the executable.
#
EBASE  = sphere

CEXE_sources += $(EBASE).cpp

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

DEFINES += -DBL_NOLINEVALUES -DBL_PARALLEL_IO


include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/EB/Make.package
include $(AMREX_HOME)/Src/GeometryShop/Make.package
include $(AMREX_HOME)/Src/Amr/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package

include $(AMREX_HOME)/Src/Extern/amrdata/Make.package
include $(AMREX_HOME)/Tools/C_util/Make.package
include ../src/Make.package

INCLUDE_LOCATIONS += .
vpathdir += .

INCLUDE_LOCATIONS += ../src
vpathdir += ../src

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
vpathdir += $(AMREX_HOME)/Src/Base


INCLUDE_LOCATIONS += $(AMREX_HOME)/Tools/C_util
vpathdir += $(AMREX_HOME)/Tools/C_util

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Boundary
vpathdir += $(AMREX_HOME)/Src/Boundary

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/EB
vpathdir += $(AMREX_HOME)/Src/EB

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/GeometryShop
vpathdir += $(AMREX_HOME)/Src/GeometryShop

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Amr
vpathdir += $(AMREX_HOME)/Src/Amr


INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/AmrCore
vpathdir += $(AMREX_HOME)/Src/AmrCore

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata
vpathdir += $(AMREX_HOME)/Src/Extern/amrdata


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





