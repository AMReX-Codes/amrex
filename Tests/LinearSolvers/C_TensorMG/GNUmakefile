PRECISION = DOUBLE
DEBUG	= FALSE
DEBUG	= TRUE
PROFILE = FALSE

DIM	= 3
DIM	= 2

COMP = KCC

PBOXLIB_HOME = ../..

USE_MPI = FALSE
USE_NETCDF=FALSE
USE_ARRAYVIEW = FALSE
USE_ARRAYVIEW = TRUE

EBASE = mcmg
LBASE = 

HERE = .

include $(PBOXLIB_HOME)/mk/Make.defs

include $(HERE)/Make.package

# Get BoxLib stuff
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/include
LIBRARY_LOCATIONS += $(PBOXLIB_HOME)/lib/$(machineSuffix)
LIBRARIES += -lmcmg$(DIM)d -lmg$(DIM)d -lamr$(DIM)d -lbndry$(DIM)d -lbox$(DIM)d

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


all: $(executable)

$(executable): $(LIBRARIES)

libs:
	cd $(PBOXLIB_HOME)/pBoxLib_2; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install
	cd $(PBOXLIB_HOME)/bndrylib; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install
	cd $(PBOXLIB_HOME)/amrlib; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install
	cd $(PBOXLIB_HOME)/mglib; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install
	cd $(PBOXLIB_HOME)/tensorMG; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install

vpath %.cpp : . .. 
vpath %.F   : . .. 
vpath %.a   $(LIBRARY_LOCATIONS)

include $(PBOXLIB_HOME)/mk/Make.rules
