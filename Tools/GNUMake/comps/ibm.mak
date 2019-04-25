
ifndef AMREX_CCOMP
  AMREX_CCOMP = ibm
endif

ifndef AMREX_FCOMP
  AMREX_FCOMP = ibm
endif

########################################################################

ibm_version  = $(shell $(CXX) --version | head -1)

COMP_VERSION = $(ibm_version)

########################################################################

GENERIC_IBM_FLAGS =

ifeq ($(USE_OMP),TRUE)
  GENERIC_IBM_FLAGS += -qsmp=omp
endif

########################################################################
########################################################################
########################################################################

ifeq ($(AMREX_CCOMP),ibm)

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

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0
  CFLAGS   += -g -O0

else

  CXXFLAGS += -g -O2 -qsimd=auto -qmaxmem=-1
  CFLAGS   += -g -O2 -qsimd=auto -qmaxmem=-1

endif

########################################################################

CXXFLAGS += -std=c++1y
CFLAGS   += -std=gnu99

########################################################################

CXXFLAGS += -Wunknown-pragmas
CFLAGS   += -Wunknown-pragmas

########################################################################

CXXFLAGS += $(GENERIC_IBM_FLAGS)
CFLAGS   += $(GENERIC_IBM_FLAGS)

endif # AMREX_CCOMP == ibm

########################################################################
########################################################################
########################################################################

ifeq ($(AMREX_FCOMP),ibm)

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

FMODULES = -qmoddir=$(fmoddir) -I $(fmoddir)

FFLAGS   += $(GENERIC_IBM_FLAGS)
F90FLAGS += $(GENERIC_IBM_FLAGS)

CPP_PREFIX = -WF,

override XTRALIBS = -lstdc++ -libmc++ -lxlf90_r -lm -lxlfmath

ifeq ($(USE_OMP),TRUE)
  override XTRALIBS += -lxlsmp
endif

ifeq ($(USE_MPI),TRUE)
  override XTRALIBS += $(shell mpifort -showme:link)
endif

FORTLINK = LOWERCASE

ifeq ($(USE_CUDA),TRUE)
  F90FLAGS += -qcuda -qtgtarch=sm_$(CUDA_ARCH)
  FFLAGS += -qcuda -qtgtarch=sm_$(CUDA_ARCH)

  ifdef CUDA_MAXREGCOUNT
    F90FLAGS += -Xptxas -maxrregcount=$(CUDA_MAXREGCOUNT)
    FFLAGS   += -Xptxas -maxrregcount=$(CUDA_MAXREGCOUNT)
  endif

  DEFINES += -DAMREX_USE_CUDA_FORTRAN

  LINK_WITH_FORTRAN_COMPILER = TRUE
endif

LINK_WITH_FORTRAN_COMPILER ?= $(USE_F_INTERFACES)

endif
