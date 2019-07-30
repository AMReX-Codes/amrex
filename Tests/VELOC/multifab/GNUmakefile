AMREX_HOME ?= ../../..

DEBUG = FALSE
DIM = 3
COMP = gnu

USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
TINY_PROFILE = TRUE
USE_VELOC = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

ifeq ($(USE_VELOC),TRUE)
  VELOC_DIR ?= /path/to/veloc/installation/dir
  DEFINES += -DAMREX_USE_VELOC
  VPATH_LOCATIONS += $(VELOC_DIR)/include
  INCLUDE_LOCATIONS += $(VELOC_DIR)/include
  LIBRARY_LOCATIONS += $(VELOC_DIR)/lib
  libraries += -Xlinker -rpath -Xlinker $(VELOC_DIR)/lib -lveloc-client
endif

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
