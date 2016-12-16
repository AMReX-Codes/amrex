# AMREX_HOME defines the directory in which we will find all the BoxLib code
# If you set AMREX_HOME as an environment variable, this line will be ignored
#AMREX_HOME ?= ../../../amrex
AMREX_HOME = ../..

DEBUG     = FALSE
DEBUG     = TRUE
USE_MPI   = TRUE
USE_OMP   = FALSE
PROFILE   = FALSE
COMP      = g++
DIM       = 3
PRECISION = DOUBLE
EBASE     = sidecarresizetest
EBASE     = gridmovetest
EBASE     = destmftest
EBASE     = testranksets
EBASE     = nsidecarstest

DO_IN_TRANSIT_ANALYSIS = TRUE
ANALYSIS_DIR = $(AMREX_HOME)/Tutorials/Sidecar_EX1

include ./Make.package
include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/Base/Make.package

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Boundary
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Amr

ifeq ($(DO_IN_TRANSIT_ANALYSIS), TRUE)
  DEFINES += -DIN_TRANSIT
  INCLUDE_LOCATIONS += $(ANALYSIS_DIR)
  vpathdir += $(ANALYSIS_DIR)
endif

vpathdir += $(AMREX_HOME)/Src/Boundary
vpathdir += $(AMREX_HOME)/Src/Base
vpathdir += $(AMREX_HOME)/Src/Amr

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
