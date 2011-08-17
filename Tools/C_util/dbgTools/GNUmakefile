BOXLIB_HOME = ../../

TOP = $(BOXLIB_HOME)
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
include $(TOP)/Tools/C_mk/Make.defs ./Make.package

include $(TOP)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += . 
INCLUDE_LOCATIONS += $(TOP)/Src/C_BaseLib 

vpath %.H   . $(TOP)/Source/C_BaseLib
vpath %.cpp . $(TOP)/Source/C_BaseLib
vpath %.F   . $(TOP)/Source/C_BaseLib
vpath %.f   . $(TOP)/Source/C_BaseLib

all: $(executable)

$(executable): $(LIBRARIES)

include $(TOP)/mk/Make.rules
