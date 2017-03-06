
AMREX_HOME ?= ../../../../../amrex

BL_NOFAST=TRUE

EBASE = HyptermOnly

BL_HAS_FORT = TRUE


DEBUG	  = TRUE
DEBUG	  = FALSE
PROFILE   = FALSE
USE_THREADS =TRUE
USE_THREADS =FALSE
DIM       = 2
DIM       = 3
USE_MPI     =FALSE
USE_MPI     =TRUE
AIXLARGEMEM =TRUE
AIXLARGEMEM =FALSE
COMP      = g++
USE_ARRAYVIEW = FALSE
USE_OMP = TRUE

HERE = .

include $(AMREX_HOME)/Tools/GNUMake/Make.defs


override CXXOPTF = -O3 -msse3
override COPTF   = -O3 -msse3
override FOPTF   = -O3 -msse3
override fOPTF   = -O3 -msse3

ifeq ($(USE_OMP),TRUE)
  override CXXOPTF += -fopenmp
  override COPTF   += -fopenmp
  override FOPTF   += -fopenmp
  override fOPTF   += -fopenmp
endif


CEXE_sources += $(EBASE).cpp

include ./Make.package

DEFINES += -DBL_PARALLEL_IO

include $(AMREX_HOME)/Src/Base/Make.package

INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Amr
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Boundary


include $(AMREX_HOME)/Src/Boundary/Make.package
INCLUDES += $(AMREX_HOME)/Src/Boundary

PATHDIRS  = $(HERE)
PATHDIRS += $(AMREX_HOME)/Src/Base
PATHDIRS += $(AMREX_HOME)/Src/Amr
PATHDIRS += $(AMREX_HOME)/Src/Extern/amrdata
PATHDIRS += $(AMREX_HOME)/Src/Boundary



ifeq (${BL_HAS_FORT}, TRUE)
FEXE_sources += ${EBASE}_F.F
endif

vpath %.f   $(PATHDIRS)
vpath %.F   $(PATHDIRS)
vpath %.h   $(PATHDIRS)
vpath %.H   $(PATHDIRS)
vpath %.cpp $(PATHDIRS)

all: $(executable)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules


