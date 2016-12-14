# AMREX_HOME defines the directory in which we will find all the BoxLib code
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../amrex

DEBUG     = FALSE
DEBUG     = TRUE
USE_MPI   = TRUE
USE_OMP   = FALSE
PROFILE   = TRUE
COMM_PROFILE  = TRUE
TRACE_PROFILE = TRUE
COMP      = g++
#COMP      = Intel
DIM       = 3
PRECISION = DOUBLE
EBASE     = mfreadwrite
EBASE     = gridmovetest

include ./Make.package
include $(AMREX_HOME)/Tools/C_mk/Make.defs
include $(AMREX_HOME)/Src/C_BoundaryLib/Make.package
include $(AMREX_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/C_BoundaryLib
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/C_BaseLib
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/C_AMRLib

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
