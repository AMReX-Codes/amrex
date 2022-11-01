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

AMREX_CCOMP = cray
AMREX_FCOMP = cray

########################################################################

ifneq ($(shell CC --version | grep -E "LLVM|clang"),)
  CRAY_IS_CLANG_BASED = TRUE
else
  CRAY_IS_CLANG_BASED = FALSE
endif

ifeq ($(CRAY_IS_CLANG_BASED),FALSE)
  # -MMD -MP not supported
  USE_LEGACY_DEPFLAGS = TRUE
  DEPFLAGS =
  LEGACY_DEPFLAGS = -M
endif

########################################################################

ifeq ($(DEBUG),TRUE)

  ifeq ($(CRAY_IS_CLANG_BASED),TRUE)
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
  ifeq ($(CRAY_IS_CLANG_BASED),TRUE)
    # The LLVM optimizer is not as aggressive as the native Cray optimizer from
    # CCE <= 8. So we adjust some flags to achieve similar optimization. See
    # this page:
    # http://pubs.cray.com/content/S-5212/9.0/cray-compiling-environment-cce-release-overview/cce-900-software-enhancements
    CXXFLAGS += -O3 -ffast-math #-fsave-loopmark -fsave-decompile
    CFLAGS   += -O3 -ffast-math #-fsave-loopmark -fsave-decompile
    FFLAGS   += -O3 -h list=a
    F90FLAGS += -O3 -h list=a
  else
    GENERIC_COMP_FLAGS += -h list=a

    CXXFLAGS += -O2
    CFLAGS   += -O2
    FFLAGS   += -O2
    F90FLAGS += -O2
  endif

endif

########################################################################

ifdef CXXSTD
  CXXSTD := $(strip $(CXXSTD))
else
  CXXSTD := c++17
endif

ifeq ($(CRAY_IS_CLANG_BASED),TRUE)
  CXXFLAGS += -std=$(CXXSTD)
  CFLAGS   += -std=c11
else
  CXXFLAGS += -h std=$(CXXSTD)
  CFLAGS   += -h c11
endif

F90FLAGS += -N 255 -em
FFLAGS   += -N 255 -em

FMODULES = -I $(fmoddir) -J $(fmoddir)

########################################################################

ifeq ($(USE_OMP),TRUE)
  # Starting in CCE 9, OpenMP is disabled by default in each of C/C++/Fortran
  # compilers.
  ifeq ($(CRAY_IS_CLANG_BASED),TRUE)
    CXXFLAGS += -fopenmp
    CFLAGS   += -fopenmp
    FFLAGS   += -h omp
    F90FLAGS += -h omp
  else
    GENERIC_COMP_FLAGS += -h omp
  endif
else
  ifeq ($(CRAY_IS_CLANG_BASED),FALSE)
    GENERIC_COMP_FLAGS += -h noomp
  endif
endif

ifeq ($(USE_ACC),TRUE)
  # OpenACC is removed from CCE altogether in CCE 9.
  ifeq ($(CRAY_IS_CLANG_BASED),TRUE)
    $(error OpenACC has been removed from CCE >= 9.)
  endif
else
  ifeq ($(CRAY_IS_CLANG_BASED),FALSE)
    GENERIC_COMP_FLAGS += -h noacc
  endif
endif

CXXFLAGS += $(GENERIC_COMP_FLAGS)
CFLAGS   += $(GENERIC_COMP_FLAGS)
FFLAGS   += $(GENERIC_COMP_FLAGS)
F90FLAGS += $(GENERIC_COMP_FLAGS)
