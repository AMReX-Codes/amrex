AMREX_HOME = ../../../amrex

PRECISION = DOUBLE
DEBUG     = TRUE
DEBUG     = FALSE
DIM       = 3
COMP      = gnu

USE_MPI=FALSE
USE_MPI=TRUE

USE_OMP=FALSE

USE_CXX11 = TRUE

PROFILE       = TRUE
TRACE_PROFILE = TRUE
COMM_PROFILE  = TRUE
MEM_PROFILE   = TRUE
#TINY_PROFILE  = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

EBASE = iotest

HERE = .

include $(AMREX_HOME)/Src/Base/Make.package

INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base

CEXE_sources += IOTestDriver.cpp IOTest.cpp

vpath %.H   $(HERE) $(AMREX_HOME)/Src/Base
vpath %.cpp $(HERE) $(AMREX_HOME)/Src/Base
vpath %.F   $(HERE) $(AMREX_HOME)/Src/Base
vpath %.f   $(HERE) $(AMREX_HOME)/Src/Base
vpath %.f90 $(HERE) $(AMREX_HOME)/Src/Base

all: $(executable)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

