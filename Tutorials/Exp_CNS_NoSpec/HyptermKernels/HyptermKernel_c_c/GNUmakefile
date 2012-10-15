
BOXLIB_HOME ?= $(HOME)/Development/BoxLib

EBASE = HyptermOnly

PRECISION = DOUBLE
PROFILE   = TRUE
PROFILE   = FALSE

COMP      = gcc
FCOMP     = nofort

DEBUG     = TRUE
DEBUG     = FALSE
DIM       = 3

USE_MPI=TRUE
USE_OMP=TRUE


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


HERE = .

INCLUDE_LOCATIONS += $(HERE)

cEXE_sources += bench.c timer.x86.c
cEXE_sources += FakeWriteMultifab.c


vpath %.cpp $(HERE)
vpath %.a $(LIBRARY_LOCATIONS)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules


