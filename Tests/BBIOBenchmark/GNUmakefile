AMREX_HOME = ../../../amrex

PROFILE   = FALSE
PRECISION = DOUBLE
DEBUG     = TRUE
DEBUG     = FALSE
DIM       = 3
COMP      = Intel
COMP      = g++

USE_MPI=FALSE
USE_MPI=TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

EBASE = bbiotest

HERE = .

include $(AMREX_HOME)/Src/Base/Make.package

INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base

CEXE_sources += BBIOTestDriver.cpp BBIOTest.cpp

vpath %.H   $(HERE) $(AMREX_HOME)/Src/Base
vpath %.cpp $(HERE) $(AMREX_HOME)/Src/Base
vpath %.F   $(HERE) $(AMREX_HOME)/Src/Base
vpath %.f   $(HERE) $(AMREX_HOME)/Src/Base
vpath %.f90 $(HERE) $(AMREX_HOME)/Src/Base

all: $(executable)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

