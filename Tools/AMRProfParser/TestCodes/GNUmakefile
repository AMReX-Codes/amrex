AMREX_HOME ?= ../../../../amrex

PROFILE        = TRUE
COMM_PROFILE   = TRUE
TRACE_PROFILE  = TRUE
PRECISION = DOUBLE
DEBUG     = FALSE
DEBUG     = TRUE
DIM       = 3
DIM       = 2
COMP      = gnu 

USE_MPI=TRUE

EBASE = pieces
EBASE = blprofcalltimes
EBASE = stest0
EBASE = AMRProfTest0
EBASE = TokenizeTest
EBASE = ProfWaitTest

ifeq ($(EBASE), blprofcalltimes)
  CEXE_sources += BLProfCallTimes.cpp
#  DEFINES += -DBL_PROF_TRACKREDUCTIONS
endif

ifeq ($(EBASE), AMRProfTest0)
  CEXE_sources += AMRProfTest0.cpp
endif

ifeq ($(EBASE), TokenizeTest)
  CEXE_sources += TokenizeTest.cpp
  USE_MPI=FALSE
  PROFILE        = FALSE
  COMM_PROFILE   = FALSE
  TRACE_PROFILE  = FALSE
endif

ifeq ($(EBASE), stest0)
  CEXE_sources += SendTest0.cpp
endif

ifeq ($(EBASE), pieces)
  DIM = 2
  CEXE_sources += Pieces.cpp
  CEXE_sources += RegionsProfStats.cpp
  CEXE_sources += CommProfStats.cpp
  CEXE_sources += BLProfStats.cpp
  CEXE_sources += BLProfUtilities.cpp
  FEXE_sources += AVGDOWN_$(DIM)D.F
endif

ifeq ($(EBASE), ProfWaitTest)
  CEXE_sources += ProfWaitTest.cpp
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

HERE = .

include $(AMREX_HOME)/Src/Base/Make.package

#CXXFLAGS += -std=c++0x

INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += ..
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base

VPATH_LOCATIONS   += $(HERE)
VPATH_LOCATIONS   += ..
VPATH_LOCATIONS   += $(AMREX_HOME)/Src/Base
#VPATH_LOCATIONS   += $(AMREX_HOME)/Src/Amr

vpath %.H   $(VPATH_LOCATIONS)
vpath %.cpp $(VPATH_LOCATIONS)
vpath %.F   $(VPATH_LOCATIONS)
vpath %.f   $(VPATH_LOCATIONS)
vpath %.f90 $(VPATH_LOCATIONS)

all: $(executable)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
