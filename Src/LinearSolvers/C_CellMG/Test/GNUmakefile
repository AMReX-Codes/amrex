PRECISION = DOUBLE
DEBUG	= FALSE
DEBUG	= TRUE
PROFILE = FALSE

DIM	= 2
DIM	= 3

COMP = KCC

PBOXLIB_HOME = ../..

USE_NETCDF=FALSE
USE_ARRAYVIEW = TRUE
USE_ARRAYVIEW = FALSE

EBASE = mg
LBASE = 

HERE = .

include $(PBOXLIB_HOME)/mk/Make.defs

CPPFLAGS += -DBL_USE_NEW_HFILES

include $(HERE)/Make.package

ifeq ($(DEBUG),TRUE)
ifeq ($(USE_ARRAYVIEW),TRUE)
ARRAYVIEWDIR  = $(PBOXLIB_HOME)/util
INCLUDE_LOCATIONS += $(ARRAYVIEWDIR)
CEXE_headers += DatasetClient.H ArrayView.H
CEXE_sources += DatasetClient.cpp
DEFINES += -DBL_USE_ARRAYVIEW
#DEFINES += -DBL_ARRAYVIEW_TAGBOX
endif
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/util
CEXE_headers += TV_TempWrite.H
endif

ifeq ($(USE_NETCDF),TRUE)
LIBRARIES += /usr/people/stevens/bin/libnetcdf.a
INCLUDE_LOCATIONS += /usr/people/stevens/bin
endif


###### exception library (for newest bsplib)
ifeq ($(MACHINE), OSF1)
LIBRARY_LOCATIONS += /usr/ccs/lib/cmplrs/cc
LIBRARIES += -lexc 
endif

3RD = 1
3RD=
ifdef 3RD
# FOR RUNNING 3RD ONLY
LDFLAGS += --link_command_prefix 3rd
#CXXDEBF = +K0 --link_command_prefix 3rd -non_shared
LDFLAGS += -non_shared -v
#LIBRARIES += -ldnet_stub
#FDEBF += -automatic
# FOR RUNNING 3RD ONLY
endif

CXXFLAGS +=
CXXOPTF +=
CXXDEBF +=

CFLAGS +=
COPTF +=
CDEBF +=

FOPTF = -fast

vpath %.cpp : . $(PBOXLIB_HOME)/pBoxLib_2 $(PBOXLIB_HOME)/amrlib $(PBOXLIB_HOME)/bndrylib $(PBOXLIB_HOME)/mglib $(UTILITY_DIR) 
vpath %.F   : . $(PBOXLIB_HOME)/pBoxLib_2 $(PBOXLIB_HOME)/amrlib $(PBOXLIB_HOME)/bndrylib $(PBOXLIB_HOME)/mglib

all: $(executable)

godzillaLink:
	rsh godzilla 'cd $(PWD); make $(executable)

include $(PBOXLIB_HOME)/mk/Make.rules
