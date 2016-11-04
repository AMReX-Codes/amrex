#
# Generic setup for using Intel
#
CXX = icpc
CC  = icc
FC  = ifort
F90 = ifort

CXXFLAGS =
CFLAGS   =
FFLAGS   =
F90FLAGS =

########################################################################

intel_version := $(shell $(CXX) -V 2>&1 1>/dev/null | grep Version)

########################################################################

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0 -Wcheck
  CFLAGS   += -g -O0 -Wcheck
  FFLAGS   += -g -O0 -CB     # -u -check all -warn all
  F90FLAGS += -g -O0 -CB     # -u -check all -warn all

else

  CXXFLAGS += -g -O2 -ip -qopt-report=5 -qopt-report-phase=vec
  CFLAGS   += -g -O2 -ip -qopt-report=5 -qopt-report-phase=vec
  FFLAGS   += -g -O2 -ip -qopt-report=5 -qopt-report-phase=vec
  F90FLAGS += -g -O2 -ip -qopt-report=5 -qopt-report-phase=vec

endif

########################################################################

CXXFLAGS += -std=c++11
CFLAGS   += -std=c99

F90FLAGS += -module $(fmoddir) -I$(fmoddir)
FFLAGS   += -module $(fmoddir) -I$(fmoddir)

########################################################################

GENERIC_COMP_FLAGS =

ifeq ($(USE_OMP),TRUE)
  GENERIC_COMP_FLAGS += -openmp
endif

CXXFLAGS += $(GENERIC_COMP_FLAGS)
CFLAGS   += $(GENERIC_COMP_FLAGS)
FFLAGS   += $(GENERIC_COMP_FLAGS)
F90FLAGS += $(GENERIC_COMP_FLAGS)

########################################################################

override XTRALIBS += -lifcore

