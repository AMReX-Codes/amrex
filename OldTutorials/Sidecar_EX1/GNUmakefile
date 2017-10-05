# AMREX_HOME defines the directory in which we will find all the BoxLib code
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../..

DEBUG     = FALSE
DEBUG     = TRUE
USE_MPI   = TRUE
USE_OMP   = FALSE
PROFILE   = FALSE
COMP      = g++
#COMP      = intel
DIM       = 3
PRECISION = DOUBLE
EBASE     = sidecarresizetest
EBASE     = gridmovetest
EBASE     = destmftest
EBASE     = testranksets
EBASE     = nsidecarstest

DO_IN_TRANSIT_ANALYSIS = TRUE
ANALYSIS_DIR = $(AMREX_HOME)/OldTutorials/Sidecar_EX1

include ./Make.package
include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include $(AMREX_HOME)/Src/Base/Make.package

ifeq ($(DO_IN_TRANSIT_ANALYSIS), TRUE)
  DEFINES += -DIN_TRANSIT
  INCLUDE_LOCATIONS += $(ANALYSIS_DIR)
  vpathdir += $(ANALYSIS_DIR)
endif

all: $(executable)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
