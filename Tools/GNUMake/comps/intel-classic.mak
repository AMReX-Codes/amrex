#
# Generic setup for using Intel classic compiler
#
CXX = icpc
CC  = icc
FC  = ifort
F90 = ifort

CXXFLAGS =
CFLAGS   =
FFLAGS   =
F90FLAGS =

AMREX_CCOMP = intel-classic
AMREX_FCOMP = intel-classic
lowercase_comp = intel-classic

########################################################################

intel_version = $(shell $(CXX) -dumpversion)

COMP_VERSION = $(intel_version)

########################################################################

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0 -traceback -Wcheck
  CFLAGS   += -g -O0 -traceback -Wcheck
  FFLAGS   += -g -O0 -traceback -check bounds,uninit,pointers
  F90FLAGS += -g -O0 -traceback -check bounds,uninit,pointers

else

  CXXFLAGS += -g -O2 -ip -qopt-report=5 -qopt-report-phase=vec
  CFLAGS   += -g -O2 -ip -qopt-report=5 -qopt-report-phase=vec
  FFLAGS   += -g -O2 -ip -qopt-report=5 -qopt-report-phase=vec
  F90FLAGS += -g -O2 -ip -qopt-report=5 -qopt-report-phase=vec

endif

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

GENERIC_COMP_FLAGS =

ifeq ($(USE_OMP),TRUE)
  GENERIC_COMP_FLAGS += -qopenmp
endif

CXXFLAGS += $(GENERIC_COMP_FLAGS) -pthread
CFLAGS   += $(GENERIC_COMP_FLAGS)
FFLAGS   += $(GENERIC_COMP_FLAGS)
F90FLAGS += $(GENERIC_COMP_FLAGS)

########################################################################

override XTRALIBS += -lifcore

ifeq ($(USE_OMP),TRUE)
  override XTRALIBS += -lifcoremt
endif

LINK_WITH_FORTRAN_COMPILER ?= $(USE_F_INTERFACES)

ifeq ($(LINK_WITH_FORTRAN_COMPILER),TRUE)
  override XTRALIBS += -lstdc++
endif
