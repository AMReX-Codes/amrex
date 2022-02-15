# --------------------------------------------------------------
# GNUmakefile  for the amrprofparser
# --------------------------------------------------------------
AMREX_HOME ?= ../../../amrex
POSTPROC_HOME = 

# NOTE: USE_PROFPARSER should always be on.
# This IS the profparser!
USE_PROFPARSER  = TRUE
PRECISION       = DOUBLE
COMM_PROFILE    = FALSE 
GPROF           = FALSE
COMP            = intel
COMP           = gnu
DEBUG           = FALSE
DEBUG           = TRUE

# ---- compile in 2d even for data from a 3d simulation
DIM           = 2

USE_MPI       = TRUE
USE_MPI       = FALSE
USE_OMP       = FALSE
EBASE         = amrprofparser
BL_NO_FORT    = FALSE
USE_CXX11     = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include $(AMREX_HOME)/Src/Base/Make.package

CEXE_sources += BLProfParser.cpp

INCLUDE_LOCATIONS += .
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/ProfParser
INCLUDE_LOCATIONS += $(POSTPROC_HOME)

VPATH_LOCATIONS   += .
VPATH_LOCATIONS   += $(AMREX_HOME)/Src/Base
VPATH_LOCATIONS   += $(AMREX_HOME)/Src/Extern/amrdata
VPATH_LOCATIONS   += $(AMREX_HOME)/Src/Extern/ProfParser
VPATH_LOCATIONS   += $(POSTPROC_HOME)


ifeq ($(__gcc_major_version),4)
  ifeq ($(__gcc_minor_version),9)
    LIBRARIES += -lquadmath
  endif
endif

ifeq ($(__gcc_major_version),5)
  #ifeq ($(__gcc_minor_version),9)
    LIBRARIES += -lquadmath
  #endif
endif

#SED0 = | sed 's/\#define vout/\/\//'
#SED1 = | sed 's/vout/\/\//'
SED0 = 
SED1 = 

vpath %.c   $(VPATH_LOCATIONS)
vpath %.cpp $(VPATH_LOCATIONS)
vpath %.h   $(VPATH_LOCATIONS)
vpath %.H   $(VPATH_LOCATIONS)
vpath %.F   $(VPATH_LOCATIONS)
vpath %.f   $(VPATH_LOCATIONS)
vpath %.f90 $(VPATH_LOCATIONS)
vpath %.l   $(VPATH_LOCATIONS)
vpath %.y   $(VPATH_LOCATIONS)

all:	$(executable)

$(executable):

cclean:
	rm -f amrprofparser BLProfParser.lex.yy.cpp BLProfParser.tab.cpp \
	      BLProfParser.tab.H core


include $(AMREX_HOME)/Tools/GNUMake/Make.rules


