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
include $(AMREX_HOME)/Tools/C_mk/Make.defs
include $(AMREX_HOME)/Src/C_BoundaryLib/Make.package
include $(AMREX_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/C_BoundaryLib
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/C_BaseLib
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/C_AMRLib

ifeq ($(DO_IN_TRANSIT_ANALYSIS), TRUE)
  DEFINES += -DIN_TRANSIT
  INCLUDE_LOCATIONS += $(ANALYSIS_DIR)
  vpathdir += $(ANALYSIS_DIR)
endif

vpathdir += $(AMREX_HOME)/Src/C_BoundaryLib
vpathdir += $(AMREX_HOME)/Src/C_BaseLib
vpathdir += $(AMREX_HOME)/Src/C_AMRLib

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(AMREX_HOME)/Tools/C_mk/Make.rules
