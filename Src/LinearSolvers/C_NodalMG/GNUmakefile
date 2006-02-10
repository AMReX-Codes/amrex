#
# $Id: GNUmakefile,v 1.54 2006-02-10 17:20:22 car Exp $
#
PBOXLIB_HOME = ..

TOP = $(PBOXLIB_HOME)

PRECISION = DOUBLE
DEBUG     = FALSE
DEBUG     = TRUE
PROFILE   = TRUE
PROFILE   = FALSE
GPROF     =
GPROF     = TRUE
DIM       = 2
DIM       = 3
COMP      = cxx
COMP	  = g++
COMP	  = Intel
FCOMP	  = g77
FCOMP     = Intel
USE_MPI   = TRUE
USE_MPI   = FALSE

# To build the executable: EBASE=proj

EBASE = proj

include $(TOP)/mk/Make.defs
include ./Make.package

include $(TOP)/bndrylib/Make.package
include $(TOP)/BoxLib/Make.package

CEXE_sources      += proj.cpp

INCLUDE_LOCATIONS += $(TOP)/bndrylib
INCLUDE_LOCATIONS += $(TOP)/BoxLib

vpath %.cpp	. $(TOP)/BoxLib $(TOP)/bndrylib
vpath %.H       . $(TOP)/BoxLib $(TOP)/bndrylib
vpath %.F	. $(TOP)/BoxLib $(TOP)/bndrylib
vpath %.f	. $(TOP)/BoxLib $(TOP)/bndrylib

all: $(executable)

include $(TOP)/mk/Make.rules
