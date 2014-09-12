
BOXLIB_HOME ?= ../../../../../BoxLib

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
FCOMP     = gfortran
USE_ARRAYVIEW = FALSE
USE_OMP = TRUE

HERE = .

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs


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

include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Extern/amrdata
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BoundaryLib


include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
INCLUDES += $(BOXLIB_HOME)/Src/C_BoundaryLib

PATHDIRS  = $(HERE)
PATHDIRS += $(BOXLIB_HOME)/Src/C_BaseLib
PATHDIRS += $(BOXLIB_HOME)/Src/C_AMRLib
PATHDIRS += $(BOXLIB_HOME)/Src/Extern/amrdata
PATHDIRS += $(BOXLIB_HOME)/Src/C_BoundaryLib



ifeq (${BL_HAS_FORT}, TRUE)
FEXE_sources += ${EBASE}_F.F
endif

vpath %.f   $(PATHDIRS)
vpath %.F   $(PATHDIRS)
vpath %.h   $(PATHDIRS)
vpath %.H   $(PATHDIRS)
vpath %.cpp $(PATHDIRS)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules


