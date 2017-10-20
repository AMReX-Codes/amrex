
AMREX_HOME ?= ../../../../../

# compilation options
COMP    = g++

USE_MPI   = TRUE
USE_UPCXX = FALSE

USE_OMP   = FALSE

PROFILE   = FALSE

PRECISION = DOUBLE

DEBUG     = FALSE

DIM      = 3
include ../../../arch.common
INCLUDE_LOCATIONS += $(INCLUDE)
INCLUDE_LOCATIONS +=  $(ROOT_PATH)/graph
INCLUDE_LOCATIONS +=  $(ROOT_PATH)/AMFIter
override XTRALIBS+= $(ROOT_PATH)/graph/graph.a $(RTS_DIR)/rts.a $(ROOT_PATH)/AMFIter/AMReX_AsyncMFIter.a -lpthread

Bpack := ./Make.package

EBASE = SMC

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

all: $(executable)
	@echo SUCCESS

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

#include $(AMREX_HOME)/Tools/GNUMake/Make.rules
#include $(AMREX_HOME)/Tools/GNUMake/Make.defs
#include ./Make.package
#include $(AMREX_HOME)/Src/Base/Make.package

all: $(executable)
	@echo SUCCESS




#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)

