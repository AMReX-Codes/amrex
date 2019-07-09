AMREX_HOME ?= ../../..

DEBUG = FALSE
DIM = 3
COMP = gnu

USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
TINY_PROFILE = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

VELOC_DIR ?= /path/to/veloc/installation/dir
VPATH_LOCATIONS += $(VELOC_DIR)/include
INCLUDE_LOCATIONS += $(VELOC_DIR)/include
LIBRARY_LOCATIONS += $(VELOC_DIR)/lib
LIBRARIES += -Xlinker -rpath -Xlinker $(VELOC_DIR)/lib -lveloc-client

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
