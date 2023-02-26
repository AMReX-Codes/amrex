#
# Generic setup for using Intel LLVM compiler
#
CXX = icpx
CC  = icx
FC  = ifx
F90 = ifx

CXXFLAGS =
CFLAGS   =
FFLAGS   =
F90FLAGS =

AMREX_CCOMP = intel-llvm
AMREX_FCOMP = intel-llvm
lowercase_comp = intel-llvm

########################################################################

intel_version = $(shell $(CXX) -dumpversion)

COMP_VERSION = $(intel_version)

########################################################################

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0 -ftrapv
  CFLAGS   += -g -O0 -ftrapv
  FFLAGS   += -g -O0 -ftrapuv -check bounds,pointers,uninit -traceback
  F90FLAGS += -g -O0 -ftrapuv -check bounds,pointers,uninit -traceback

else

  CXXFLAGS += -g1 -O3
  CFLAGS   += -g1 -O3
  FFLAGS   += -g -O3
  F90FLAGS += -g -O3

endif

########################################################################

ifeq ($(WARN_ALL),TRUE)
  warning_flags = -Wall -Wextra -Wno-sign-compare -Wunreachable-code -Wnull-dereference
  warning_flags += -Wfloat-conversion -Wextra-semi

  ifneq ($(USE_CUDA),TRUE)
    warning_flags += -Wpedantic
  endif

  ifneq ($(WARN_SHADOW),FALSE)
    warning_flags += -Wshadow
  endif

  CXXFLAGS += $(warning_flags) -Woverloaded-virtual -Wnon-virtual-dtor
  CFLAGS += $(warning_flags)
endif

ifeq ($(WARN_ERROR),TRUE)
  CXXFLAGS += -Werror
  CFLAGS += -Werror
endif

CXXFLAGS += -Wno-tautological-constant-compare
CFLAGS += -Wno-tautological-constant-compare

########################################################################

ifdef CXXSTD
  CXXSTD := $(strip $(CXXSTD))
  CXXFLAGS += -std=$(CXXSTD)
else
  CXXFLAGS += -std=c++17
endif

CFLAGS   += -std=c11

F90FLAGS += -implicitnone

FMODULES = -module $(fmoddir) -I$(fmoddir)

########################################################################

# ../../../../Src/Base/AMReX_GpuUtility.H:148:16: warning: comparison with NaN
# always evaluates to false in fast floating point modes
# [-Wtautological-constant-compare]
#        return std::isnan(m);
#               ^~~~~~~~~~~~~
GENERIC_COMP_FLAGS =

ifeq ($(USE_OMP),TRUE)
  GENERIC_COMP_FLAGS += -qopenmp
endif

CXXFLAGS += $(GENERIC_COMP_FLAGS) -pthread
CFLAGS   += $(GENERIC_COMP_FLAGS)
FFLAGS   += $(GENERIC_COMP_FLAGS)
F90FLAGS += $(GENERIC_COMP_FLAGS)

########################################################################

ifneq ($(BL_NO_FORT),TRUE)
  override XTRALIBS += -lifcore

  ifeq ($(USE_OMP),TRUE)
    override XTRALIBS += -lifcoremt
  endif

  LINK_WITH_FORTRAN_COMPILER ?= $(USE_F_INTERFACES)

  ifeq ($(LINK_WITH_FORTRAN_COMPILER),TRUE)
    override XTRALIBS += -lstdc++
  endif
endif
