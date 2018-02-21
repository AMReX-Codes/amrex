AMREX_HOME ?= ../../../../../../..
USE_EB = FALSE
PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

DIM        = 2
#DIM       = 3

COMP	   = gnu

USE_PARTICLES = TRUE

USE_MPI    = TRUE
USE_OMP    = TRUE

Bpack   := ./Make.package 
Blocs   := . 


include ../../../../../arch.common
INCLUDE_LOCATIONS += $(INCLUDE)
INCLUDE_LOCATIONS +=  $(ROOT_PATH)/graph
INCLUDE_LOCATIONS +=  $(ROOT_PATH)/AMFIter
override XTRALIBS+= $(ROOT_PATH)/graph/graph.a $(RTS_DIR)/rts.a $(ROOT_PATH)/AMFIter/AMReX_AsyncMFIter.a


include ../Make.Adv
