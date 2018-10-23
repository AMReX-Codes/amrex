
DIM          = 3

COMP         = gnu

DEBUG        = FALSE

USE_MPI      = TRUE
USE_OMP      = FALSE

AMREX_HOME = ../..

EBASE = main

include ./Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base

vpathdir += $(AMREX_HOME)/Src/Base

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
