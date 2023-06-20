
ifndef AMREX_CCOMP
  AMREX_CCOMP = nvhpc
endif

ifndef AMREX_FCOMP
  AMREX_FCOMP = nvhpc
endif

########################################################################

nvhpc_version = $(shell $(CXX) -V 2>&1 | grep 'target' | sed 's|.*$(CXX) \([0-9\.]*\).*|\1|')
nvhpc_major_version = $(shell echo $(nvhpc_version) | cut -f1 -d.)
nvhpc_minor_version = $(shell echo $(nvhpc_version) | cut -f2 -d.)

gcc_version       = $(shell g++ -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;')
gcc_major_version = $(shell g++ -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;' | sed -e 's;\..*;;')
gcc_minor_version = $(shell g++ -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;' | sed -e 's;[^.]*\.;;' | sed -e 's;\..*;;')

COMP_VERSION = $(nvhpc_version)

# -MP not supported by nvhpc and -MMD's output is put in the wrong directory
USE_LEGACY_DEPFLAGS = TRUE
DEPFLAGS =

########################################################################

GENERIC_NVHPC_FLAGS =

ifeq ($(USE_OMP),TRUE)
  GENERIC_NVHPC_FLAGS += -mp -Minfo=mp
endif

ifeq ($(USE_OMP_OFFLOAD),TRUE)
  # CUDA + NVIDIA OpenMP-offload requires CUDA 11.0 or later:
  # "nvfortran-Fatal-Use of -Mcuda and -mp=gpu requires CUDA 11.0 or later"
  # My Cori GPU interactive tests worked without needing to specify
  # a CUDA version: I used nvhpc/21.3 and cuda/11.1.1 modules.
  # (The USE_ACC code path below should be revisited)
  GENERIC_NVHPC_FLAGS += -mp=gpu -Minfo=mp
  ifneq ($(CUDA_ARCH),)
    GENERIC_NVHPC_FLAGS += -gpu=cc$(CUDA_ARCH)
  endif
endif

ifeq ($(USE_ACC),TRUE)
  GENERIC_NVHPC_FLAGS += -acc=gpu -Minfo=accel -mcmodel=medium
  ifneq ($(CUDA_ARCH),)
    GENERIC_NVHPC_FLAGS += -gpu=cc$(CUDA_ARCH)
  endif
endif

# Note that -O2 is the default optimization level for NVHPC

NVHPC_OPT := -O2 -fast

########################################################################
########################################################################
########################################################################

ifeq ($(AMREX_CCOMP),nvhpc)

CXX = nvc++
CC  = nvc

########################################################################

CXXFLAGS =
CFLAGS   =

# Allow -gopt to be disabled to work around a compiler bug on P9.

NVHPC_GOPT ?= TRUE

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0 -Mbounds
  CFLAGS   += -g -O0 -Mbounds

else

  CXXFLAGS += $(NVHPC_OPT)
  CFLAGS   += $(NVHPC_OPT)

  ifeq ($(NVHPC_GOPT),TRUE)

    CXXFLAGS += -gopt
    CFLAGS   += -gopt

  endif

endif

# The logic here should be consistent with what's in nvcc.mak
ifdef CXXSTD
  CXXSTD := $(strip $(CXXSTD))
  ifeq ($(shell expr $(gcc_major_version) \< 8),1)
    $(error GCC >= 8 required.)
  endif
  CXXFLAGS += -std=$(CXXSTD)
else
  CXXFLAGS += -std=c++17
endif

CFLAGS   += -c11

CXXFLAGS += $(GENERIC_NVHPC_FLAGS)
CFLAGS   += $(GENERIC_NVHPC_FLAGS)

else # AMREX_CCOMP == nvhpc

# If we're using OpenACC but also CUDA, then nvcc will be the C++ compiler. If
# we want to call the OpenACC API from C++ then we need to make sure we have
# the includes for it, because NVHPC may not be the host compiler for nvcc.

ifeq ($(USE_ACC),TRUE)
  NVHPC_BIN_LOCATION := $(shell nvc++ -show 2>&1 | grep CPPCOMPDIR | awk '{print $$7}' | cut -c2-)
  NVHPC_LOCATION := $(shell dirname $(NVHPC_BIN_LOCATION))
  INCLUDE_LOCATIONS += $(NVHPC_LOCATION)/etc/include_acc
endif

endif # AMREX_CCOMP == nvhpc

########################################################################
########################################################################
########################################################################

ifeq ($(AMREX_FCOMP),nvhpc)

#
# Now set the Fortran flags. Since this is done after the GNU include
# in the CUDA version, all of the GNU specific options are overridden.
#

FC  = nvfortran
F90 = nvfortran

FFLAGS   =
F90FLAGS =

ifeq ($(DEBUG),TRUE)

  FFLAGS   += -g -O0 -Mbounds -Ktrap=divz,inv -Mchkptr
  F90FLAGS += -g -O0 -Mbounds -Ktrap=divz,inv -Mchkptr

else

  FFLAGS   += $(NVHPC_OPT)
  F90FLAGS += $(NVHPC_OPT)

  ifeq ($(NVHPC_GOPT),TRUE)

    FFLAGS   += -gopt
    F90FLAGS += -gopt

  endif

endif

# Note that we do not have a Fortran main

ifneq ($(USE_F_INTERFACES),TRUE)
  F90FLAGS += -Mnomain
  FFLAGS   += -Mnomain
endif

ifeq ($(USE_CUDA),TRUE)

  F90FLAGS += -gpu=cc$(CUDA_ARCH),fastmath
  FFLAGS   += -gpu=cc$(CUDA_ARCH),fastmath

  ifneq ($(DEBUG),TRUE)
    F90FLAGS += -Mcuda=lineinfo
    FFLAGS   += -Mcuda=lineinfo
  endif

  ifeq ($(CUDA_VERBOSE),TRUE)
    F90FLAGS += -gpu=keepptx
    FFLAGS   += -gpu=keepptx
  endif

  F90FLAGS += CUDA_HOME=$(COMPILE_CUDA_PATH)
  FFLAGS   += CUDA_HOME=$(COMPILE_CUDA_PATH)

  ifdef CUDA_MAXREGCOUNT
    F90FLAGS += -gpu=maxregcount:$(CUDA_MAXREGCOUNT)
    FFLAGS   += -gpu=maxregcount:$(CUDA_MAXREGCOUNT)
  endif

  ifeq ($(USE_MPI),TRUE)
  ifneq ($(findstring Open MPI, $(shell mpicxx -showme:version 2>&1)),)
    OMPI_FCFLAGS_ORIG = $(shell mpif90 -showme:compile)
    export OMPI_FCFLAGS := $(subst -pthread,-lpthread,$(OMPI_FCFLAGS_ORIG))
  endif
  endif

endif


########################################################################

F90FLAGS += -Mdclchk
FFLAGS   += -Mextend

FMODULES = -module $(fmoddir) -I$(fmoddir)

########################################################################

FFLAGS   += $(GENERIC_NVHPC_FLAGS)
F90FLAGS += $(GENERIC_NVHPC_FLAGS)

########################################################################

# Add -lrt for the missing "aio_return" symbol
# /usr/common/software/sles15_cgpu/nvhpc/21.3/Linux_x86_64/21.3/compilers/lib/libnvf.a:async.o:                 U aio_return
override XTRALIBS += -lstdc++ -latomic -lnvf -lrt

LINK_WITH_FORTRAN_COMPILER ?= $(USE_F_INTERFACES)

endif # AMREX_FCOMP == nvhpc
