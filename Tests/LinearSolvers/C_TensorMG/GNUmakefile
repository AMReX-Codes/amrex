PRECISION = DOUBLE
DEBUG	= FALSE
DEBUG	= TRUE
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

EBASE = mcmg
LBASE = 

HERE = .

include $(PBOXLIB_HOME)/mk/Make.defs

CPPFLAGS += -DBL_USE_NEW_HFILES

INCLUDE_LOCATIONS += . $(PBOXLIB_HOME)/pBoxLib_2 $(PBOXLIB_HOME)/amrlib
INCLUDE_LOCATIONS +=  $(PBOXLIB_HOME)/tensorMG
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/bndrylib $(PBOXLIB_HOME)/mglib

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

#
# Uncomment this is you want to use Parallel I/O.
#
DEFINES += -DBL_PARALLEL_IO

ifeq ($(USE_ARRAYVIEW),TRUE)
#ARRAYVIEWDIR  = /usr/local/ccse/ArrayView
#ARRAYVIEWDIR  = /usr/people/vince/Visualization/ArrayView
ARRAYVIEWDIR  = .
INCLUDE_LOCATIONS += $(ARRAYVIEWDIR)
#LIBRARY_LOCATIONS += $(ARRAYVIEWDIR) 
#LIBRARIES += -larrayview$(DIM)d.$(machineSuffix) 
DEFINES += -DBL_USE_ARRAYVIEW
DEFINES += -DBL_ARRAYVIEW_TAGBOX
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

#CXXFLAGS = -g --diag_suppress 177
#CXXFLAGS = --strict_warnings
ifeq ($(BSP_MACHINE), OSF1)
#FFLAGS += -real_size 64
#FDEBF += -C
#FDEBF += - -warn argument_checking
#FDEBF += - -warn declarations
#FDEBF += - -warn truncated_source
#FDEBF += - -warn unused
#FOPTF  = -fast -O5 -tune ev5
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

include $(HERE)/Make.package

FOPTF = -fast

vpath %.cpp : . $(PBOXLIB_HOME)/pBoxLib_2 $(PBOXLIB_HOME)/amrlib $(PBOXLIB_HOME)/tensorMG $(PBOXLIB_HOME)/bndrylib $(PBOXLIB_HOME)/mglib
vpath %.F   : . $(PBOXLIB_HOME)/pBoxLib_2 $(PBOXLIB_HOME)/amrlib $(PBOXLIB_HOME)/tensorMG $(PBOXLIB_HOME)/bndrylib $(PBOXLIB_HOME)/mglib

all: $(executable)

godzillaLink:
	rsh godzilla 'cd $(PWD); make $(executable)

include $(PBOXLIB_HOME)/mk/Make.rules
