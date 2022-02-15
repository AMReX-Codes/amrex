#
# Generic setup for using ncc
#
CXX = nc++
CC  = ncc
FC  = nfort
F90 = nfort

CXXFLAGS =
CFLAGS   =
FFLAGS   =
F90FLAGS =

########################################################################

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0 -fno-inline -ftrace -Wall -Wunused
  CFLAGS   += -g -O0 -fno-inline -ftrace -Wall -Wunused

  FFLAGS   += -g -O0 -fcheck=bounds -ftrace -Wuninitialized
  F90FLAGS += -g -O0 -fcheck=bounds -ftrace -Wuninitialized

else

  CXXFLAGS += -g -O3
  CFLAGS   += -g -O3
  FFLAGS   += -g -O3
  F90FLAGS += -g -O3

endif

########################################################################

FMODULES += -module $(fmoddir) -I$(fmoddir)

########################################################################

GENERIC_COMP_FLAGS =

ifeq ($(USE_OMP),TRUE)
  GENERIC_COMP_FLAGS += -fopenmp
endif

CXXFLAGS += $(GENERIC_COMP_FLAGS)
CFLAGS   += $(GENERIC_COMP_FLAGS)
FFLAGS   += $(GENERIC_COMP_FLAGS)
F90FLAGS += $(GENERIC_COMP_FLAGS)
LDFLAGS  += -cxxlib

########################################################################

