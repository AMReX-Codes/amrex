PBOXLIB_HOME := ../..

PROFILE = FALSE
PRECISION = DOUBLE
DEBUG	= FALSE
DEBUG	= TRUE
DIM	= 2
DIM	= 3
COMP    = Intel
FCOMP   = Intel
USE_MPI = FALSE

USE_ARRAYVIEW = FALSE

EBASE = main

include $(PBOXLIB_HOME)/mk/Make.defs

CXXFLAGS +=
CXXOPTF  +=
CXXDEBF  +=

include Make.package

include $(PBOXLIB_HOME)/tensorMG/Make.package
include $(PBOXLIB_HOME)/bndrylib/Make.package
include $(PBOXLIB_HOME)/mglib/Make.package
include $(PBOXLIB_HOME)/BoxLib/Make.package

INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/tensorMG
vpathdir += $(PBOXLIB_HOME)/tensorMG
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/bndrylib
vpathdir += $(PBOXLIB_HOME)/bndrylib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/BoxLib
vpathdir += $(PBOXLIB_HOME)/BoxLib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/iamrlib
vpathdir += $(PBOXLIB_HOME)/iamrlib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/mglib
vpathdir += $(PBOXLIB_HOME)/mglib

# We need these from amrlib 

CEXE_sources += BCRec.cpp 
CEXE_headers += BCRec.H BC_TYPES.H
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/amrlib
vpathdir += $(PBOXLIB_HOME)/amrlib

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(PBOXLIB_HOME)/mk/Make.rules


