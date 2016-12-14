AMREX_HOME ?= ../../../amrex

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

include $(AMREX_HOME)/Tools/C_mk/Make.defs

EBASE = bbiotest

HERE = .

include $(AMREX_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/C_BaseLib

CEXE_sources += BBIOTestDriver.cpp BBIOTest.cpp

vpath %.H   $(HERE) $(AMREX_HOME)/Src/C_BaseLib
vpath %.cpp $(HERE) $(AMREX_HOME)/Src/C_BaseLib
vpath %.F   $(HERE) $(AMREX_HOME)/Src/C_BaseLib
vpath %.f   $(HERE) $(AMREX_HOME)/Src/C_BaseLib
vpath %.f90 $(HERE) $(AMREX_HOME)/Src/C_BaseLib

all: $(executable)

include $(AMREX_HOME)/Tools/C_mk/Make.rules

