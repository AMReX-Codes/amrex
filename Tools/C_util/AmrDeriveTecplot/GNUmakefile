AMREX_HOME ?= ../../../../amrex

BL_NOFAST   = TRUE
EBASE       = AmrDeriveTecplot
DEBUG	    = TRUE
PROFILE     = FALSE
USE_THREADS = FALSE
DIM         = 2
USE_MPI     = FALSE
AIXLARGEMEM = FALSE
COMP        = g++
BINARY_OUTPUT = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

DEFINES     += -DBL_NOLINEVALUES

CEXE_sources += $(EBASE).cpp
CEXE_headers += DataServices.H AmrData.H XYPlotDataList.H AmrvisConstants.H
CEXE_sources += DataServices.cpp AmrData.cpp

ifeq ($(DIM),3)
FEXE_sources += FABUTIL_3D.F
else
ifeq ($(DIM),2)
FEXE_sources += FABUTIL_2D.F
endif
endif



DEFINES += -DBL_PARALLEL_IO

include $(AMREX_HOME)/Src/Base/Make.package

INCLUDE_LOCATIONS += .
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata

PATHDIRS  = .
PATHDIRS += $(AMREX_HOME)/Src/Base
PATHDIRS += $(AMREX_HOME)/Src/Extern/amrdata


ifeq ($(BINARY_OUTPUT),TRUE)
  DEFINES += -DUSE_TEC_BIN_IO
  TECHOME=/usr/local/tecplot
  INCLUDE_LOCATIONS += ${TECHOME}/include
  LIBRARIES += $(TECHOME)/lib/tecio64.a
endif

vpath %.f   $(PATHDIRS)
vpath %.F   $(PATHDIRS)
vpath %.h   $(PATHDIRS)
vpath %.H   $(PATHDIRS)
vpath %.cpp $(PATHDIRS)
vpath %.c   $(PATHDIRS)
vpath %.f90 $(PATHDIRS)

all: $(executable)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules



