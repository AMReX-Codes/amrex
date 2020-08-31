AMREX_HOME = ../../../

DEBUG     = FALSE
USE_MPI   = TRUE
USE_OMP   = FALSE
USE_CUDA  = TRUE
COMP      = pgi
DIM       = 3

Bpack   := ./Make.package
Blocs   := .

EBASE := main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

TOP := $(AMREX_HOME)/Tests/Stream
include $(TOP)/Make.package
INCLUDE_LOCATIONS += $(TOP)
VPATH_LOCATIONS   += $(TOP)

include $(AMREX_HOME)/Src/Base/Make.package

all: $(executable)
	@echo SUCCESS

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
