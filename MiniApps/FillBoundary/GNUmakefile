BOXLIB_HOME = ../../../BoxLib.commprof

DEBUG        = FALSE
USE_MPI      = TRUE
USE_OMP      = FALSE
USE_IPM      = FALSE
PROFILE      = TRUE
COMM_PROFILE = TRUE
USE_DUMPI    = FALSE
COMP         = g++
FCOMP        = gfortran
#COMP         = Intel
#FCOMP        = Intel
DIM          = 3
PRECISION    = DOUBLE
EBASE        = fbtest

CEXE_sources += MultiFabFillBoundaryTest.cpp

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs
include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BoundaryLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib

ifeq ($(COMP), g++)
  CXXFLAGS += -std=c++11
  DEFINES += -DBL_CXX11
endif

ifeq ($(USE_DUMPI), TRUE)
  LIBRARY_LOCATIONS += /project/projectdirs/m888/DUMPI_LIB/lib
  LIBRARIES += -ldumpi 
endif

vpathdir += $(BOXLIB_HOME)/Src/C_BoundaryLib
vpathdir += $(BOXLIB_HOME)/Src/C_BaseLib
vpathdir += $(BOXLIB_HOME)/Src/C_AMRLib

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
