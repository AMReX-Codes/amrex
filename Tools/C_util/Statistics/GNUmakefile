#
# $Id: GNUmakefile,v 1.5 2009-10-07 21:24:49 gpau Exp $
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
#DIM           = 3
COMP          = Intel
FCOMP         = Intel
USE_MPI       = FALSE
USE_MPI       = TRUE
#
# The only stencils in which we're interested when using USE_HGPROJ_SERIAL.
#
ifeq ($(USE_HGPROJ_SERIAL),TRUE)
DEFINES += -DBL_USE_HGPROJ_SERIAL
ifeq ($(DIM),2)
PRVERSION = v5
else
PRVERSION = v7
endif
endif
#
# Base name of the executable.
#
EBASE = PltFileXAve
#EBASE = PltFileStat
#EBASE = PltFileList
CEXE_sources += $(EBASE).cpp

include $(TOP)/mk/Make.defs ./Make.package
XTRALIBS += -ldl

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
