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

# CRAY_CC_VERSION is defined by the 'cce' module. Starting with CCE 9, Cray
# changed the C and C++ compilers to clang/LLVM based, so all of the
# options/flags for those compilers changed. But the Cray Fortran compiler is
# still based on CCE 8 and so it has the same options as before.
COMP_VERSION = $(shell echo $(CRAY_CC_VERSION) | cut -f 1 -d .)

########################################################################

ifeq ($(DEBUG),TRUE)

  ifeq ($(COMP_VERSION),9)
    CXXFLAGS += -g -O0
    CFLAGS   += -g -O0
    FFLAGS   += -g -O0 -e i -K trap=fp
    F90FLAGS += -g -O0 -e i -K trap=fp
  else
    GENERIC_COMP_FLAGS += -K trap=fp

    CXXFLAGS += -g -O0
    CFLAGS   += -g -O0
    FFLAGS   += -g -O0 -e i
    F90FLAGS += -g -O0 -e i
  endif

else
  ifeq ($(COMP_VERSION),9)
    # The LLVM optimizer is not as aggressive as the native Cray optimizer from
    # CCE <= 8. So we adjust some flags to achieve similar optimization. See
    # this page:
    # http://pubs.cray.com/content/S-5212/9.0/cray-compiling-environment-cce-release-overview/cce-900-software-enhancements
    CXXFLAGS += -O2 -ffast-math -fsave-loopmark -fsave-decompile
    CFLAGS   += -O2 -ffast-math -fsave-loopmark -fsave-decompile
    FFLAGS   += -O2 -h list=a
    F90FLAGS += -O2 -h list=a
  else
    GENERIC_COMP_FLAGS += -h list=a

    CXXFLAGS += -O2
    CFLAGS   += -O2
    FFLAGS   += -O2
    F90FLAGS += -O2
  endif

endif

########################################################################

ifeq ($(COMP_VERSION),9)
  CXXFLAGS += -std=c++11
  CFLAGS   += -std=c99
else
  CXXFLAGS += -h std=c++11
  CFLAGS   += -h c99
endif

F90FLAGS += -N 255 -em
FFLAGS   += -N 255 -em

FMODULES = -I $(fmoddir) -J $(fmoddir)

########################################################################

ifeq ($(USE_OMP),TRUE)
  # Starting in CCE 9, OpenMP is disabled by default in each of C/C++/Fortran
  # compilers.
  ifeq ($(COMP_VERSION),9)
    CXXFLAGS += -fopenmp
    CFLAGS   += -fopenmp
    FFLAGS   += -h omp
    F90FLAGS += -h omp
  else
    GENERIC_COMP_FLAGS += -h omp
  endif
else
  ifneq ($(COMP_VERSION),9)
    GENERIC_COMP_FLAGS += -h noomp
  endif
endif

ifeq ($(USE_ACC),TRUE)
  # OpenACC is removed from CCE altogether in CCE 9.
  ifeq ($(COMP_VERSION),9)
    $(error OpenACC has been removed from CCE 9.)
  endif
else
  ifneq ($(COMP_VERSION),9)
    GENERIC_COMP_FLAGS += -h noacc
  endif
endif

CXXFLAGS += $(GENERIC_COMP_FLAGS)
CFLAGS   += $(GENERIC_COMP_FLAGS)
FFLAGS   += $(GENERIC_COMP_FLAGS)
F90FLAGS += $(GENERIC_COMP_FLAGS)
