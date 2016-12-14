AMREX_HOME = ../../../amrex

PRECISION = DOUBLE
DEBUG     = FALSE
DEBUG     = TRUE
DIM       = 3
COMP      = g++
FCOMP     = gfortran

USE_MPI=FALSE
USE_MPI=TRUE

USE_CXX11 = TRUE

#PROFILE       = TRUE
#TRACE_PROFILE = TRUE
#COMM_PROFILE  = FALSE

include $(AMREX_HOME)/Tools/C_mk/Make.defs

EBASE = iotest

HERE = .

include $(AMREX_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/C_BaseLib

CEXE_sources += IOTestDriver.cpp IOTest.cpp

vpath %.H   $(HERE) $(AMREX_HOME)/Src/C_BaseLib
vpath %.cpp $(HERE) $(AMREX_HOME)/Src/C_BaseLib
vpath %.F   $(HERE) $(AMREX_HOME)/Src/C_BaseLib
vpath %.f   $(HERE) $(AMREX_HOME)/Src/C_BaseLib
vpath %.f90 $(HERE) $(AMREX_HOME)/Src/C_BaseLib

all: $(executable)

include $(AMREX_HOME)/Tools/C_mk/Make.rules

