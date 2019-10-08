AMREX_HOME ?= ../../..

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 2
DIM	= 3

COMP    = gcc

PRECISION = DOUBLE

#########################################################
#NOTE: this version of diffmultifab only works in serial;
#      we must build it with USE_MPI = FALSE
#########################################################

USE_MPI   = FALSE
USE_OMP   = FALSE

TEST=TRUE
USE_ASSERTION=TRUE

###################################################

EBASE     = diffmultifab

# If NEEDS_f90_SRC=TRUE, look for ${EBASE}_nd.f90
NEEDS_f90_SRC = FALSE
#NEEDS_f90_SRC = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += ${EBASE}.cpp

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
include $(AMREX_HOME)/Src/Base/Make.package
vpathdir += $(AMREX_HOME)/Src/Base

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata
include $(AMREX_HOME)/Src/Extern/amrdata/Make.package
vpathdir += $(AMREX_HOME)/Src/Extern/amrdata

INCLUDE_LOCATIONS += $(AMREX_HOME)/Tools/C_util
include $(AMREX_HOME)/Tools/C_util/Make.package
vpathdir += $(AMREX_HOME)/Tools/C_util

ifeq ($(NEEDS_f90_SRC),TRUE)
  f90EXE_sources += ${EBASE}_nd.f90
endif

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
