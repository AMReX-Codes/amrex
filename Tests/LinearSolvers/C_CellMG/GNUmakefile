PRECISION = DOUBLE
DEBUG	= TRUE
DEBUG	= FALSE
PROFILE = FALSE

DIM	= 2
DIM	= 3

COMP = KCC

PBOXLIB_HOME = ../..

USE_BSP=FALSE
USE_BSP=TRUE
USE_NETCDF=FALSE
USE_ARRAYVIEW = FALSE
USE_ARRAYVIEW = TRUE

#
# What is the type of BSP device?
#
ifeq ($(USE_BSP),TRUE)
BSP_DEVICE = SHMEM_SYSV
endif

EBASE = mg
LBASE = 

HERE = .

include $(PBOXLIB_HOME)/mk/Make.defs

CPPFLAGS += -DBL_USE_NEW_HFILES

ifeq ($(USE_BSP), TRUE)
DEFINES += -DBL_USE_BSP
ifeq ($(BSP_MACHINE), OSF1)
BSP_HOME = /usr/people/vince/Parallel/BSP/BSP
endif
ifeq ($(MACHINE), T3E)
ifeq ($(WHICHT3E), NERSC)
DEFINES += -DBL_T3E_NERSC
BSP_HOME = /u1/vince/BSP
endif
ifeq ($(WHICHT3E), NAVO)
DEFINES += -DBL_T3E_NAVO
BSP_HOME = /home/Cvince/BSP
endif
endif
endif

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

ifeq ($(USE_BSP), TRUE)
INCLUDE_LOCATIONS += $(BSP_HOME)/include
LIBRARY_LOCATIONS += $(BSP_HOME)/lib/$(BSP_MACHINE)
LIBRARY_LOCATIONS += $(BSP_HOME)/lib/$(BSP_MACHINE)/$(BSP_DEVICE)
endif

ifeq ($(USE_NETCDF),TRUE)
LIBRARIES += /usr/people/stevens/bin/libnetcdf.a
INCLUDE_LOCATIONS += /usr/people/stevens/bin
endif

ifeq ($(USE_BSP), TRUE)
ifeq ($(DEBUG), TRUE)
LIBRARIES += -lbspcore_O0 -lbsplevel1_O0
else
ifeq ($(BSP_MACHINE), OSF1)
LIBRARIES += -lbspcore_O2 -lbsplevel1_O0
else
LIBRARIES += -lbspcore_O2 -lbsplevel1_O2
endif
endif
endif

###### exception library (for newest bsplib)
ifeq ($(BSP_MACHINE), OSF1)
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

ifeq ($(USE_BSP),TRUE)
vpath %.cpp : . $(PBOXLIB_HOME)/pBoxLib_2 $(PBOXLIB_HOME)/amrlib $(PBOXLIB_HOME)/bndrylib $(PBOXLIB_HOME)/mglib $(UTILITY_DIR) $(PRELOAD_DIR)
else
vpath %.cpp : . $(PBOXLIB_HOME)/pBoxLib_2 $(PBOXLIB_HOME)/amrlib $(PBOXLIB_HOME)/bndrylib $(PBOXLIB_HOME)/mglib $(UTILITY_DIR) 
endif
vpath %.F   : . $(PBOXLIB_HOME)/pBoxLib_2 $(PBOXLIB_HOME)/amrlib $(PBOXLIB_HOME)/bndrylib $(PBOXLIB_HOME)/mglib

all: $(executable)

godzillaLink:
	rsh godzilla 'cd $(PWD); make $(executable)

include $(PBOXLIB_HOME)/mk/Make.rules
