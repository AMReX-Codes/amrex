#
# $Id: GNUmakefile,v 1.56 2010-02-13 22:52:49 lijewski Exp $
#
PBOXLIB_HOME = ..

TOP = $(PBOXLIB_HOME)

PRECISION = DOUBLE
DEBUG     = TRUE
DEBUG     = FALSE
PROFILE   = TRUE
PROFILE   = FALSE
GPROF     = TRUE
GPROF     =
DIM       = 2
DIM       = 3
COMP      = PathScale
FCOMP     = PathScale
USE_MPI   = FALSE
USE_MPI   = TRUE
USE_OMP   = TRUE

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
