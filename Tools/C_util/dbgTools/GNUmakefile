AMREX_HOME = ../../

TOP = $(AMREX_HOME)
#
# Variables for the user to set ...
#
PRECISION     = DOUBLE
DEBUG	      = FALSE
DIM	      = 3
COMP          = KCC
#
# Base name of the executable.
#
EBASE = crsGrids
CEXE_sources += $(EBASE).cpp

DEFINES += -DBL_PARALLEL_IO
include $(TOP)/Tools/GNUMake/Make.defs ./Make.package

include $(TOP)/Src/Base/Make.package

INCLUDE_LOCATIONS += . 
INCLUDE_LOCATIONS += $(TOP)/Src/Base 

vpath %.H   . $(TOP)/Source/Base
vpath %.cpp . $(TOP)/Source/Base
vpath %.F   . $(TOP)/Source/Base
vpath %.f   . $(TOP)/Source/Base

all: $(executable)

$(executable): $(LIBRARIES)

include $(TOP)/mk/Make.rules
