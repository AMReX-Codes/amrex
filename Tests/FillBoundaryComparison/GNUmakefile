
DIM          = 3

COMP         = Intel
FCOMP        = Intel

DEBUG        = FALSE

USE_MPI      = TRUE
USE_OMP      = FALSE
USE_UPCXX    = TRUE

BOXLIB_HOME = ../..

EBASE = main

include ./Make.package

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib

vpathdir += $(BOXLIB_HOME)/Src/C_BaseLib

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
