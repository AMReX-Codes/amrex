
ifeq ($(USE_OMP),TRUE)
  CXX = xlC_r
  CC  = xlc_r
else ifeq ($(USE_CUDA),TRUE)
  CXX = xlC_r
  CC  = xlc_r
else
  CXX = xlC
  CC  = xlc
endif

CXXFLAGS =
CFLAGS   =

########################################################################

ibm_version  = $(shell $(CXX) --version | head -1)

COMP_VERSION = $(ibm_version)

########################################################################

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0
  CFLAGS   += -g -O0

else

  CXXFLAGS += -g -O2
  CFLAGS   += -g -O2

endif

########################################################################

CXXFLAGS += -std=c++1y
CFLAGS   += -std=gnu99

########################################################################

GENERIC_IBM_FLAGS = -qmoddir=$(fmoddir) -I $(fmoddir)

ifeq ($(USE_OMP),TRUE)
  GENERIC_COMP_FLAGS += -qsmp=omp
endif

CXXFLAGS += $(GENERIC_IBM_FLAGS)
CFLAGS   += $(GENERIC_IBM_FLAGS)

CPP_PREFIX = -WF,

ifeq ($(USE_CUDA),TRUE)
  include $(AMREX_HOME)/Tools/GNUMake/comps/gnu.mak

  # Force immediate expansion of the GCC defines,
  # since after this point GCC will no longer be
  # the actual compiler defined in CXX.

  DEFINES := $(DEFINES)

  CXXFLAGS := -Wno-deprecated-gpu-targets -x cu --std=c++11 -ccbin=$(CXX) -Xcompiler='$(CXXFLAGS)'
  CFLAGS := -Wno-deprecated-gpu-targets -x c -ccbin=$(CC) -Xcompiler='$(CFLAGS)'

  HOST_CXX := $(CXX)
  HOST_CC := $(CC)

  CXX := nvcc
  CC := nvcc
endif



ifeq ($(USE_OMP),TRUE)
  FC  = xlf_r
  F90 = xlf_r
else ifeq ($(USE_CUDA),TRUE)
  FC  = xlf_r
  F90 = xlf_r
else
  FC  = xlf
  F90 = xlf
endif

FFLAGS =
F90FLAGS =

ifeq ($(DEBUG),TRUE)

  FFLAGS   += -g -O0
  F90FLAGS += -g -O0

else

  FFLAGS   += -g -O2
  F90FLAGS += -g -O2

endif

F90FLAGS += -qlanglvl=extended -qxlf2003=polymorphic

FFLAGS   += -WF,-C!
F90FLAGS += -WF,-C!

FFLAGS   += -qfixed=72

FFLAGS   += $(GENERIC_IBM_FLAGS)
F90FLAGS += $(GENERIC_IBM_FLAGS)

override XTRALIBS = -L$(OLCF_XLF_ROOT)/lib -L$(OLCF_XLC_ROOT)/lib -lstdc++ -libmc++ -lxlf90_r -lm -lxlfmath

ifeq ($(USE_OMP),TRUE)
  override XTRALIBS += -lxlsmp
endif

ifeq ($(USE_MPI),TRUE)
  override XTRALIBS += $(shell mpifort -showme:link)
endif

FORTLINK = LOWERCASE

ifeq ($(USE_CUDA),TRUE)
  F90FLAGS += -qcuda
  FFLAGS += -qcuda
endif
