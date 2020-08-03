
ifndef AMREX_CCOMP
  AMREX_CCOMP = hpcsdk
endif

ifndef AMREX_FCOMP
  AMREX_FCOMP = hpcsdk
endif

########################################################################

hpcsdk_version = $(shell $(CXX) -V 2>&1 | grep 'target' | sed 's|.*$(CXX) \([0-9\.]*\).*|\1|')
hpcsdk_major_version = $(shell echo $(hpcsdk_version) | cut -f1 -d.)
hpcsdk_minor_version = $(shell echo $(hpcsdk_version) | cut -f2 -d.)

gcc_version       = $(shell g++ -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;')
gcc_major_version = $(shell g++ -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;' | sed -e 's;\..*;;')
gcc_minor_version = $(shell g++ -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;' | sed -e 's;[^.]*\.;;' | sed -e 's;\..*;;')

COMP_VERSION = $(hpcsdk_version)

########################################################################

GENERIC_HPCSDK_FLAGS =

ifeq ($(USE_OMP),TRUE)
  GENERIC_HPCSDK_FLAGS += -mp -Mconcur=nonuma -Minfo=mp
endif

ifeq ($(USE_ACC),TRUE)
  GENERIC_HPCSDK_FLAGS += -acc=gpu -Minfo=accel -mcmodel=medium
  ifneq ($(CUDA_ARCH),)
    # We use 10.1 because nvcc defaults to 10.1 if it can't detect a GPU
    # driver. And in Cori GPU interactive jobs, nvcc can't see the GPU driver
    # unless it is executed within an `srun`. Most people do not execute `srun
    # make`, so nvcc defaults to 10.1. But HPC SDK defaults to CUDA 11, so we
    # get link errors between nvcc and the HPC SDK if we don't do something
    # about this. The easiest fix is to simply force HPC SDK to use CUDA 10.1
    # to match the blind nvcc.
    GENERIC_HPCSDK_FLAGS += -acc=gpu -gpu=cc$(CUDA_ARCH),cuda10.1
  else
    GENERIC_HPCSDK_FLAGS += -acc=gpu
  endif
endif

# Note that -O2 is the default optimization level for HPCSDK

HPCSDK_OPT := -O2 -fast

########################################################################
########################################################################
########################################################################

ifeq ($(AMREX_CCOMP),hpcsdk)

CXX = nvc++
CC  = nvc

########################################################################

CXXFLAGS =
CFLAGS   =

# Allow -gopt to be disabled to work around a compiler bug on P9.

HPCSDK_GOPT ?= TRUE

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0 -Mbounds
  CFLAGS   += -g -O0 -Mbounds

else

  CXXFLAGS += $(HPCSDK_OPT)
  CFLAGS   += $(HPCSDK_OPT)

  ifeq ($(HPCSDK_GOPT),TRUE)

    CXXFLAGS += -gopt
    CFLAGS   += -gopt

  endif

endif

# The logic here should be consistent with what's in nvcc.mak
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

CFLAGS   += -c99

CXXFLAGS += $(GENERIC_HPCSDK_FLAGS)
CFLAGS   += $(GENERIC_HPCSDK_FLAGS)

else # AMREX_CCOMP == hpcsdk

# If we're using OpenACC but also CUDA, then nvcc will be the C++ compiler. If
# we want to call the OpenACC API from C++ then we need to make sure we have
# the includes for it, because HPCSDK may not be the host compiler for nvcc.

ifeq ($(USE_ACC),TRUE)
  HPCSDK_BIN_LOCATION := $(shell nvc++ -show 2>&1 | grep CPPCOMPDIR | awk '{print $$7}' | cut -c2-)
  HPCSDK_LOCATION := $(shell dirname $(HPCSDK_BIN_LOCATION))
  INCLUDE_LOCATIONS += $(HPCSDK_LOCATION)/etc/include_acc
endif

endif # AMREX_CCOMP == hpcsdk

########################################################################
########################################################################
########################################################################

ifeq ($(AMREX_FCOMP),hpcsdk)

#
# Now set the Fortran flags. Since this is done after the GNU include
# in the CUDA version, all of the GNU specific options are overriden.
#

FC  = nvfortran
F90 = nvfortran

FFLAGS   =
F90FLAGS =

ifeq ($(DEBUG),TRUE)

  FFLAGS   += -g -O0 -Mbounds -Ktrap=divz,inv -Mchkptr
  F90FLAGS += -g -O0 -Mbounds -Ktrap=divz,inv -Mchkptr

else

  FFLAGS   += $(HPCSDK_OPT)
  F90FLAGS += $(HPCSDK_OPT)

  ifeq ($(HPCSDK_GOPT),TRUE)

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

  F90FLAGS += -gpu=cc$(CUDA_ARCH),fastmath,cuda10.1
  FFLAGS   += -gpu=cc$(CUDA_ARCH),fastmath,cuda10.1

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

  DEFINES += -DAMREX_USE_CUDA_FORTRAN

  LINK_WITH_FORTRAN_COMPILER = TRUE

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

FFLAGS   += $(GENERIC_HPCSDK_FLAGS)
F90FLAGS += $(GENERIC_HPCSDK_FLAGS)

########################################################################

override XTRALIBS += -lstdc++ -latomic

LINK_WITH_FORTRAN_COMPILER ?= $(USE_F_INTERFACES)

endif # AMREX_FCOMP == hpcsdk
