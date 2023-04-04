# Store the CUDA toolkit version.

ifneq ($(NO_CONFIG_CHECKING),TRUE)
  nvcc_version       := $(shell nvcc --version | grep "release" | awk 'BEGIN {FS = ","} {print $$2}' | awk '{print $$2}')
  nvcc_major_version := $(shell nvcc --version | grep "release" | awk 'BEGIN {FS = ","} {print $$2}' | awk '{print $$2}' | awk 'BEGIN {FS = "."} {print $$1}')
  nvcc_minor_version := $(shell nvcc --version | grep "release" | awk 'BEGIN {FS = ","} {print $$2}' | awk '{print $$2}' | awk 'BEGIN {FS = "."} {print $$2}')
else
  nvcc_version       := 99.9
  nvcc_major_version := 99
  nvcc_minor_version := 9
endif

# Disallow CUDA toolkit versions < 11

nvcc_major_lt_11 = $(shell expr $(nvcc_major_version) \< 11)
ifeq ($(nvcc_major_lt_11),1)
  $(error Your nvcc version is $(nvcc_version). This is unsupported. Please use CUDA toolkit version 11.0 or newer.)
endif

ifeq ($(shell expr $(nvcc_major_version) \= 11),1)
ifeq ($(shell expr $(nvcc_minor_version) \= 0),1)
  # -MP not supported in 11.0
  DEPFLAGS = -MMD
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

  ifeq ($(shell expr $(gcc_major_version) \< 8),1)
    $(error GCC >= 8 required.)
  endif

  ifdef CXXSTD
    CXXSTD := $(strip $(CXXSTD))
  else
    CXXSTD = c++17
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
  else
    CXXSTD := c++17
  endif

  CXXFLAGS += -std=$(CXXSTD)

  NVCC_CCBIN ?= pgc++

  # In pgi.make, we use gcc_major_version to handle c++17 flag.
  CXXFLAGS_FROM_HOST := -ccbin=$(NVCC_CCBIN) -Xcompiler='$(CXXFLAGS)' --std=$(CXXSTD)
  CFLAGS_FROM_HOST := $(CXXFLAGS_FROM_HOST)
else
  ifdef CXXSTD
    CXXSTD := $(strip $(CXXSTD))
  else
    CXXSTD := c++17
  endif

  NVCC_CCBIN ?= $(CXX)

  CXXFLAGS_FROM_HOST := -ccbin=$(NVCC_CCBIN) -Xcompiler='$(CXXFLAGS)' --std=$(CXXSTD)
  CFLAGS_FROM_HOST := $(CXXFLAGS_FROM_HOST)
endif

NVCC_FLAGS = -Wno-deprecated-gpu-targets -m64 -arch=compute_$(CUDA_ARCH) -code=sm_$(CUDA_ARCH) -maxrregcount=$(CUDA_MAXREGCOUNT) --expt-relaxed-constexpr --expt-extended-lambda --forward-unknown-to-host-compiler
# This is to work around a bug with nvcc, see: https://github.com/kokkos/kokkos/issues/1473
NVCC_FLAGS += -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored

ifeq ($(GPU_ERROR_CROSS_EXECUTION_SPACE_CALL),TRUE)
  NVCC_FLAGS += --Werror cross-execution-space-call
endif

ifeq ($(DEBUG),TRUE)
  NVCC_FLAGS += -g -G
else
  NVCC_FLAGS += -lineinfo --ptxas-options=-O3
endif

ifeq ($(CUDA_VERBOSE),TRUE)
  NVCC_FLAGS += --ptxas-options=-v
endif

ifeq ($(USE_CUPTI),TRUE)
  SYSTEM_INCLUDE_LOCATIONS += $(MAKE_CUDA_PATH)/extras/CUPTI/include
  LIBRARY_LOCATIONS += ${MAKE_CUDA_PATH}/extras/CUPTI/lib64
  LIBRARIES += -Wl,-rpath,${MAKE_CUDA_PATH}/extras/CUPTI/lib64 -lcupti
endif

ifneq ($(USE_CUDA_FAST_MATH),FALSE)
  NVCC_FLAGS += --use_fast_math
endif

NVCC_FLAGS += $(XTRA_NVCC_FLAGS)

ifeq ($(GPU_ERROR_CAPTURE_THIS),TRUE)
  NVCC_FLAGS += --Werror ext-lambda-captures-this
else
ifeq ($(GPU_WARN_CAPTURE_THIS),TRUE)
  NVCC_FLAGS += --Wext-lambda-captures-this
endif
endif

nvcc_diag_error = 0
ifeq ($(shell expr $(nvcc_major_version) \>= 12),1)
  nvcc_diag_error = 1
else
ifeq ($(shell expr $(nvcc_major_version) \= 11),1)
ifeq ($(shell expr $(nvcc_minor_version) \>= 2),1)
  nvcc_diag_error = 1
endif
endif
endif
# warning #20092-D: a __device__ variable cannot be directly written in a host function
ifeq ($(nvcc_diag_error),1)
  NVCC_FLAGS += --display-error-number --diag-error 20092
endif

CXXFLAGS = $(CXXFLAGS_FROM_HOST) $(NVCC_FLAGS) -x cu
CFLAGS   =   $(CFLAGS_FROM_HOST) $(NVCC_FLAGS) -x cu

ifeq ($(USE_GPU_RDC),TRUE)
  CXXFLAGS += -dc
  CFLAGS   += -dc
else
  CXXFLAGS += -c
  CFLAGS   += -c
endif

CXX = nvcc
CC  = nvcc
