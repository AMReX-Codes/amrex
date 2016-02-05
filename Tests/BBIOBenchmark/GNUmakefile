BOXLIB_HOME ?= ../../../BoxLib

PROFILE   = FALSE
PRECISION = DOUBLE
DEBUG     = TRUE
DEBUG     = FALSE
DIM       = 3
COMP      = Intel
FCOMP     = Intel
COMP      = g++
FCOMP     = gfortran

USE_MPI=FALSE
USE_MPI=TRUE

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

EBASE = bbiotest

HERE = .

include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib

CEXE_sources += BBIOTestDriver.cpp BBIOTest.cpp

vpath %.H   $(HERE) $(BOXLIB_HOME)/Src/C_BaseLib
vpath %.cpp $(HERE) $(BOXLIB_HOME)/Src/C_BaseLib
vpath %.F   $(HERE) $(BOXLIB_HOME)/Src/C_BaseLib
vpath %.f   $(HERE) $(BOXLIB_HOME)/Src/C_BaseLib
vpath %.f90 $(HERE) $(BOXLIB_HOME)/Src/C_BaseLib

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules

