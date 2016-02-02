# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME ?= ../../../BoxLib

DEBUG     = FALSE
DEBUG     = TRUE
USE_MPI   = TRUE
USE_OMP   = FALSE
PROFILE   = TRUE
COMM_PROFILE  = TRUE
TRACE_PROFILE = TRUE
COMP      = g++
FCOMP     = gfortran
#COMP      = Intel
#FCOMP     = Intel
DIM       = 3
PRECISION = DOUBLE
EBASE     = mfreadwrite
EBASE     = gridmovetest
USE_CXX11     = TRUE

include ./Make.package
include $(BOXLIB_HOME)/Tools/C_mk/Make.defs
include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BoundaryLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib

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
