#
# Generic setup for using Cray
#
CXX = CC
CC  = cc
FC  = ftn
F90 = ftn

CXXFLAGS =
CFLAGS   =
FFLAGS   =
F90FLAGS =

########################################################################

cray_version = $(shell $(CXX) -V 2>&1 | grep 'Version')

COMP_VERSION = $(cray_version)

########################################################################

ifeq ($(DEBUG),TRUE)

  GENERIC_COMP_FLAGS += -K trap=fp

  CXXFLAGS += -g -O0
  CFLAGS   += -g -O0
  FFLAGS   += -g -O0 -e i
  F90FLAGS += -g -O0 -e i

else

  GENERIC_COMP_FLAGS += -h list=a

  CXXFLAGS += -O2
  CFLAGS   += -O2
  FFLAGS   += -O2
  F90FLAGS += -O2

endif

########################################################################

CXXFLAGS += -h std=c++11
CFLAGS   += -h c99

F90FLAGS += -N 255 -em
FFLAGS   += -N 255 -em

FMODULES = -I $(fmoddir) -J $(fmoddir)

########################################################################

ifneq ($(USE_OMP),TRUE)
  GENERIC_COMP_FLAGS += -h noomp
endif

ifneq ($(USE_ACC),TRUE)
  GENERIC_COMP_FLAGS += -h noacc
endif

CXXFLAGS += $(GENERIC_COMP_FLAGS)
CFLAGS   += $(GENERIC_COMP_FLAGS)
FFLAGS   += $(GENERIC_COMP_FLAGS)
F90FLAGS += $(GENERIC_COMP_FLAGS)
