#
# $Id: GNUmakefile,v 1.12 2010-12-22 23:46:58 gpau Exp $
#
PBOXLIB_HOME = ../../

TOP = $(PBOXLIB_HOME)
#
# Variables for the user to set ...
#
PRECISION     = DOUBLE
DEBUG	      = TRUE
DEBUG	      = FALSE
PROFILE       = FALSE
DIM	      = 2
COMP          = Intel
FCOMP         = Intel
USE_MPI       = TRUE

#
# Base name of the executable.
#
EBASE = DiffUniform
#EBASE = PltFileNorm
#EBASE = PltFileScalConvRate
#EBASE = DiffFab
#EBASE = DiffSameGrid
#EBASE = PltFileNorm
#EBASE = PltFileNormB
#EBASE = PltFileScalConvRate
#EBASE = DiffSameGrid
#EBASE = DiffSameGrid2
#EBASE = DiffUniform
#EBASE = DiffSameGridRefined

#EBASE = DiffFab

CEXE_sources += $(EBASE).cpp

include $(TOP)/mk/Make.defs ./Make.package

ifeq ($(USE_ARRAYVIEW),TRUE)
  DEFINES += -DBL_USE_ARRAYVIEW
  DEFINES += -DBL_ARRAYVIEW_TAGBOX
endif

fincludes=$(includes)

DEFINES += -DBL_NOLINEVALUES -DBL_PARALLEL_IO

Bdirs := BoxLib bndrylib

Bpack	:= $(foreach dir, $(Bdirs), $(TOP)/$(dir)/Make.package)
Blocs	:= $(foreach dir, $(Bdirs), $(TOP)/$(dir))

VPATH_LOCATIONS   += . .. $(Blocs)
INCLUDE_LOCATIONS += . .. $(Blocs) $(TOP)/pAmrvis

include $(Bpack)

vpath %.c   $(VPATH_LOCATIONS)
vpath %.cpp $(VPATH_LOCATIONS) $(TOP)/pAmrvis
vpath %.h   $(VPATH_LOCATIONS)
vpath %.H   $(VPATH_LOCATIONS) 
vpath %.F   $(VPATH_LOCATIONS) $(TOP)/pAmrvis
vpath %.f   $(VPATH_LOCATIONS)
vpath %.f90 $(VPATH_LOCATIONS)

all: $(executable)

$(executable): $(LIBRARIES)

include $(TOP)/mk/Make.rules
