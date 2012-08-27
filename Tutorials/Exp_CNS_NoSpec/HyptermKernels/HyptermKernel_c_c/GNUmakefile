
BOXLIB_HOME ?= $(HOME)/Development/BoxLib

EBASE = htccc

PRECISION = DOUBLE
PROFILE   = TRUE
PROFILE   = FALSE

COMP      = gcc
FCOMP     = nofort

DEBUG     = TRUE
DEBUG     = FALSE
DIM       = 3

USE_MPI=FALSE


include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

HERE = .

INCLUDE_LOCATIONS += $(HERE)

cEXE_sources += bench.c timer.x86.c
cEXE_sources += FakeWriteMultifab.c

COPTF += -O3

vpath %.cpp $(HERE)
vpath %.a $(LIBRARY_LOCATIONS)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules


