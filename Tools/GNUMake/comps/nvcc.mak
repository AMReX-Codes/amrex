# Store the CUDA toolkit version.

nvcc_version       := $(shell nvcc --version | grep "release" | awk 'BEGIN {FS = ","} {print $$2}' | awk '{print $$2}')
nvcc_major_version := $(shell nvcc --version | grep "release" | awk 'BEGIN {FS = ","} {print $$2}' | awk '{print $$2}' | awk 'BEGIN {FS = "."} {print $$1}')
nvcc_minor_version := $(shell nvcc --version | grep "release" | awk 'BEGIN {FS = ","} {print $$2}' | awk '{print $$2}' | awk 'BEGIN {FS = "."} {print $$2}')

# Disallow CUDA toolkit versions < 8.0.

nvcc_major_lt_8 = $(shell expr $(nvcc_major_version) \< 8)
ifeq ($(nvcc_major_lt_8),1)
  $(error Your nvcc version is $(nvcc_version). This is unsupported. Please use CUDA toolkit version 8.0 or newer.)
endif

nvcc_forward_unknowns = 0
ifeq ($(shell expr $(nvcc_major_version) \= 10),1)
ifeq ($(shell expr $(nvcc_minor_version) \>= 2),1)
  nvcc_forward_unknowns = 1
endif
endif
ifeq ($(shell expr $(nvcc_major_version) \>= 11),1)
  nvcc_forward_unknowns = 1
endif

ifeq ($(shell expr $(nvcc_major_version) \< 10),1)
  DEPFLAGS = -M  # -MM not supported in < 10
endif

ifeq ($(shell expr $(nvcc_major_version) \= 10),1)
ifeq ($(shell expr $(nvcc_minor_version) \= 0),1)
  DEPFLAGS = -M  # -MM not supported in 10.0
endif
endif

#
# nvcc compiler driver does not always accept pgc++
# as a host compiler at present. However, if we're using
# OpenACC, we proabably need to use PGI.
#

NVCC_HOST_COMP ?= $(AMREX_CCOMP)

lowercase_nvcc_host_comp = $(shell echo $(NVCC_HOST_COMP) | tr A-Z a-z)

ifeq ($(lowercase_nvcc_host_comp),$(filter $(lowercase_nvcc_host_comp),gcc gnu g++))
  lowercase_nvcc_host_comp = gnu
  AMREX_CCOMP = gnu
  ifndef GNU_DOT_MAK_INCLUDED
    include $(AMREX_HOME)/Tools/GNUMake/comps/gnu.mak
  endif
endif

ifeq ($(lowercase_nvcc_host_comp),gnu)
  ifdef CXXSTD
    CXXSTD := $(strip $(CXXSTD))
    ifeq ($(shell expr $(gcc_major_version) \< 5),1)
      ifeq ($(CXXSTD),c++14)
        $(error C++14 support requires GCC 5 or newer.)
      endif
    endif
  else
    ifeq ($(gcc_major_version),4)
      CXXSTD := c++11
    else ifeq ($(gcc_major_version),5)
      CXXSTD := c++14
    else
      CXXSTD := c++14
    endif
  endif

  CXXFLAGS += -std=$(CXXSTD)

  NVCC_CCBIN ?= g++

  CXXFLAGS_FROM_HOST := -ccbin=$(NVCC_CCBIN) -Xcompiler='$(CXXFLAGS)' --std=$(CXXSTD)
  CFLAGS_FROM_HOST := $(CXXFLAGS_FROM_HOST)
  ifeq ($(USE_OMP),TRUE)
     LIBRARIES += -lgomp
  endif
else ifeq ($(lowercase_nvcc_host_comp),pgi)
  ifdef CXXSTD
    CXXSTD := $(strip $(CXXSTD))
    ifeq ($(shell expr $(gcc_major_version) \< 5),1)
      ifeq ($(CXXSTD),c++14)
        $(error C++14 support requires GCC 5 or newer.)
      endif
    endif
  else
    ifeq ($(gcc_major_version),4)
      CXXSTD := c++11
    else ifeq ($(gcc_major_version),5)
      CXXSTD := c++14
    else
      CXXSTD := c++14
    endif
  endif

  CXXFLAGS += -std=$(CXXSTD)

  NVCC_CCBIN ?= pgc++

  # In pgi.make, we use gcc_major_version to handle c++11/c++14 flag.
  CXXFLAGS_FROM_HOST := -ccbin=$(NVCC_CCBIN) -Xcompiler='$(CXXFLAGS)' --std=$(CXXSTD)
  CFLAGS_FROM_HOST := $(CXXFLAGS_FROM_HOST)
else
  ifdef CXXSTD
    CXXSTD := $(strip $(CXXSTD))
  else
    CXXSTD := c++11
  endif

  NVCC_CCBIN ?= $(CXX)

  CXXFLAGS_FROM_HOST := -ccbin=$(NVCC_CCBIN) -Xcompiler='$(CXXFLAGS)' --std=$(CXXSTD)
  CFLAGS_FROM_HOST := $(CXXFLAGS_FROM_HOST)
endif

NVCC_FLAGS = -Wno-deprecated-gpu-targets -m64 -arch=compute_$(CUDA_ARCH) -code=sm_$(CUDA_ARCH) -maxrregcount=$(CUDA_MAXREGCOUNT) --expt-relaxed-constexpr --expt-extended-lambda
# Unfortunately, on cori with cuda 10.0 this fails in thrust code
# NVCC_FLAGS += --Werror=cross-execution-space-call

ifeq ($(DEBUG),TRUE)
  NVCC_FLAGS += -g -G
else
  NVCC_FLAGS += -lineinfo --ptxas-options=-O3
endif

ifeq ($(CUDA_VERBOSE),TRUE)
  NVCC_FLAGS += --ptxas-options=-v
endif

ifeq ($(USE_CUPTI),TRUE)
  INCLUDE_LOCATIONS += $(MAKE_CUDA_PATH)/extras/CUPTI/include
  LIBRARY_LOCATIONS += ${MAKE_CUDA_PATH}/extras/CUPTI/lib64
  LIBRARIES += -Wl,-rpath,${MAKE_CUDA_PATH}/extras/CUPTI/lib64 -lcupti
endif

ifneq ($(USE_CUDA_FAST_MATH),FALSE)
  NVCC_FLAGS += --use_fast_math
endif

NVCC_FLAGS += $(XTRA_NVCC_FLAGS)

ifeq ($(nvcc_forward_unknowns),1)
  NVCC_FLAGS += --forward-unknown-to-host-compiler
endif

ifeq ($(shell expr $(nvcc_major_version) \>= 11),1)
ifeq ($(GPU_ERROR_CAPTURE_THIS),TRUE)
  NVCC_FLAGS += --Werror ext-lambda-captures-this
else
ifeq ($(GPU_WARN_CAPTURE_THIS),TRUE)
  NVCC_FLAGS += --Wext-lambda-captures-this
endif
endif
endif

CXXFLAGS = $(CXXFLAGS_FROM_HOST) $(NVCC_FLAGS) -dc -x cu
CFLAGS   =   $(CFLAGS_FROM_HOST) $(NVCC_FLAGS) -dc -x cu

ifeq ($(nvcc_version),9.2)
  # relaxed constexpr not supported
  ifeq ($(USE_EB),TRUE)
    $(error Cuda 9.2 is not supported with USE_EB=TRUE. Use 9.1 or 10.0 instead.)
  endif
endif

CXX = nvcc
CC  = nvcc
