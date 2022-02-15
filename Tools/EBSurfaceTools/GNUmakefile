
TOP = $(AMREX_HOME)
USE_EB = TRUE
#
# Variables for the user to set ...
#
PRECISION     = DOUBLE
DEBUG	      = TRUE
DEBUG	      = FALSE
DIM	      = 3
COMP          = g++
USE_MPI       = FALSE

EBASE  = ConvertEBSurface

CEXE_sources += $(EBASE).cpp

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

DEFINES += -DBL_NOLINEVALUES -DBL_PARALLEL_IO

Bdirs := Base EB AmrCore Boundary

Bpack	+= $(foreach dir, $(Bdirs), $(TOP)/Src/$(dir)/Make.package)
Blocs	+= $(foreach dir, $(Bdirs), $(TOP)/Src/$(dir))

include $(Bpack)

INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

INCLUDE_LOCATIONS += .
vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

