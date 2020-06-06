#
# Generic setup for using NAG
#
CXX = g++
CC  = gcc
FC  = nagfor
F90 = nagfor

CXXFLAGS =
CFLAGS   =
FFLAGS   =
F90FLAGS =

########################################################################

gcc_version       = $(shell $(CXX) -dumpversion | head -1 | sed -e 's;.*  *;;')
gcc_major_version = $(shell $(CXX) -dumpversion | head -1 | sed -e 's;.*  *;;' | sed -e 's;\..*;;')
gcc_minor_version = $(shell $(CXX) -dumpversion | head -1 | sed -e 's;.*  *;;' | sed -e 's;[^.]*\.;;' | sed -e 's;\..*;;')

COMP_VERSION = $(gcc_version)

########################################################################

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0 -fno-inline -ggdb -Wall -Wno-sign-compare -ftrapv
  CFLAGS   += -g -O0 -fno-inline -ggdb -Wall -Wno-sign-compare -ftrapv

  FFLAGS   += -g -O0 -nan -C
  F90FLAGS += -g -O0 -nan -C

else

  CXXFLAGS += -g -O3
  CFLAGS   += -g -O3
  FFLAGS   += -g -O3
  F90FLAGS += -g -O3

endif


ifeq ($(USE_GPROF),TRUE)

  CXXFLAGS += -pg
  CFLAGS += -pg
  FFLAGS += -pg
  F90FLAGS += -pg

endif

########################################################################

ifdef CXXSTD
  CXXSTD := $(strip $(CXXSTD))
  ifeq ($(shell expr $(gcc_major_version) \< 5),1)
    ifeq ($(CXXSTD),c++14)
      $(error C++14 support requires GCC 5 or newer.)
    endif
  endif
  CXXFLAGS += -std=$(CXXSTD)
else
  ifeq ($(gcc_major_version),4)
    CXXFLAGS += -std=c++11
  else ifeq ($(gcc_major_version),5)
    CXXFLAGS += -std=c++14
  endif
endif

CFLAGS   += -std=gnu99

FFLAGS   += -mismatch
F90FLAGS += -mismatch -u

FMODULES = -mdir $(fmoddir) -I $(fmoddir)

########################################################################

GENERIC_COMP_FLAGS =
GENERIC_FORT_FLAGS =

ifeq ($(THREAD_SANITIZER),TRUE)
  GENERIC_COMP_FLAGS += -fsanitize=thread
endif
ifeq ($(FSANITIZER),TRUE)
  GENERIC_COMP_FLAGS += -fsanitize=address -fsanitize=undefined
endif

ifeq ($(USE_OMP),TRUE)
  GENERIC_COMP_FLAGS += -fopenmp
  GENERIC_FORT_FLAGS += -openmp
endif

CXXFLAGS += $(GENERIC_COMP_FLAGS)
CFLAGS   += $(GENERIC_COMP_FLAGS)
FFLAGS   += $(GENERIC_FORT_FLAGS)
F90FLAGS += $(GENERIC_FORT_FLAGS)

########################################################################
