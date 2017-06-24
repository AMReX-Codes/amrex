#
# Generic setup for using PGI
#
CXX = pgc++
CC  = pgcc
FC  = pgfortran
F90 = pgfortran

CXXFLAGS =
CFLAGS   =
FFLAGS   =
F90FLAGS =

########################################################################

pgi_version := $(shell $(CXX) -V 2>&1 | grep 'target')

COMP_VERSION := $(pgi_version)

COMP_MAJOR_VERSION := $(shell echo "$(COMP_VERSION)" | awk '{print $$2}' | awk '{print substr ($$0, 0, 2)}')
COMP_MINOR_VERSION := $(shell echo "$(COMP_VERSION)" | awk '{print $$2}' | awk '{print substr ($$0, 4, 10)}')

########################################################################

ifeq ($(DEBUG),TRUE)

  # 2016-12-02: pgi 16.10 doesn't appear to like -traceback together with c++11

  CXXFLAGS += -g -O0 -Mbounds
  CFLAGS   += -g -O0 -Mbounds
  FFLAGS   += -g -O0 -Mbounds -Ktrap=divz,inv -Mchkptr
  F90FLAGS += -g -O0 -Mbounds -Ktrap=divz,inv -Mchkptr

else

  # 2017-06-23: PGI 17.4 causes the CUDA code to break when using -O2/-fast

  ifeq ($(USE_CUDA),TRUE)
    PGI_OPT := -O
  else
    PGI_OPT := -fast
  endif

  CXXFLAGS += -gopt $(PGI_OPT)
  CFLAGS   += -gopt $(PGI_OPT)
  FFLAGS   += -gopt $(PGI_OPT)
  F90FLAGS += -gopt $(PGI_OPT)

endif

########################################################################

CXXFLAGS += --c++11
CFLAGS   += -c99

F90FLAGS += -module $(fmoddir) -I$(fmoddir) -Mdclchk
FFLAGS   += -module $(fmoddir) -I$(fmoddir) -Mextend

########################################################################

GENERIC_COMP_FLAGS =

ifeq ($(USE_OMP),TRUE)
  GENERIC_COMP_FLAGS += -mp=nonuma -Minfo=mp
endif

ifeq ($(USE_ACC),TRUE)
  GENERIC_COMP_FLAGS += -acc -Minfo=acc -ta=nvidia -lcudart -mcmodel=medium
else
  GENERIC_COMP_FLAGS += -noacc
endif

ifeq ($(USE_CUDA),TRUE)
  CXXFLAGS += -Mcuda=cuda8.0
  CFLAGS   += -Mcuda=cuda8.0
  FFLAGS   += -Mcuda=cuda8.0 -Mnomain
  F90FLAGS += -Mcuda=cuda8.0 -Mnomain

  override XTRALIBS += -lstdc++
endif

CXXFLAGS += $(GENERIC_COMP_FLAGS)
CFLAGS   += $(GENERIC_COMP_FLAGS)
FFLAGS   += $(GENERIC_COMP_FLAGS)
F90FLAGS += $(GENERIC_COMP_FLAGS)

########################################################################

# Because we do not have a Fortran main

ifeq ($(which_computer),$(filter $(which_computer),summitdev))
override XTRALIBS += -pgf90libs -L /sw/summitdev/gcc/5.4.0new/lib64/ -latomic
else
override XTRALIBS += -pgf90libs -latomic
endif
