# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../../../

DEBUG     = FALSE
USE_MPI   = TRUE
USE_OMP   = TRUE
COMP      = gnu
DIM       = 3

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

include ../../../arch.common
INCLUDE_LOCATIONS += $(INCLUDE)
INCLUDE_LOCATIONS +=  $(ROOT_PATH)/graph
INCLUDE_LOCATIONS +=  $(ROOT_PATH)/AMFIter
override XTRALIBS+= $(ROOT_PATH)/graph/graph.a $(RTS_DIR)/rts.a $(ROOT_PATH)/AMFIter/AMReX_AsyncMFIter.a
