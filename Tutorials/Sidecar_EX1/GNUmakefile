# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
#BOXLIB_HOME ?= ../../../BoxLib
BOXLIB_HOME = ../..

DEBUG     = FALSE
DEBUG     = TRUE
USE_MPI   = TRUE
USE_OMP   = FALSE
PROFILE   = FALSE
COMP      = g++
FCOMP     = gfortran
DIM       = 3
PRECISION = DOUBLE
EBASE     = sidecarresizetest
EBASE     = gridmovetest
EBASE     = destmftest
EBASE     = testranksets
EBASE     = nsidecarstest

DO_IN_TRANSIT_ANALYSIS = TRUE
ANALYSIS_DIR = $(BOXLIB_HOME)/Tutorials/Sidecar_EX1

include ./Make.package
include $(BOXLIB_HOME)/Tools/C_mk/Make.defs
include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BoundaryLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib

ifeq ($(DO_IN_TRANSIT_ANALYSIS), TRUE)
  DEFINES += -DIN_TRANSIT
  INCLUDE_LOCATIONS += $(ANALYSIS_DIR)
  vpathdir += $(ANALYSIS_DIR)
endif

vpathdir += $(BOXLIB_HOME)/Src/C_BoundaryLib
vpathdir += $(BOXLIB_HOME)/Src/C_BaseLib
vpathdir += $(BOXLIB_HOME)/Src/C_AMRLib

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
