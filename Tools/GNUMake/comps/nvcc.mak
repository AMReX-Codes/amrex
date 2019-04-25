# Store the CUDA toolkit version.

nvcc_version       := $(shell nvcc --version | tail -1 | awk 'BEGIN {FS = ","} {print $$2}' | awk '{print $$2}')
nvcc_major_version := $(shell nvcc --version | tail -1 | awk 'BEGIN {FS = ","} {print $$2}' | awk '{print $$2}' | awk 'BEGIN {FS = "."} {print $$1}')
nvcc_minor_version := $(shell nvcc --version | tail -1 | awk 'BEGIN {FS = ","} {print $$2}' | awk '{print $$2}' | awk 'BEGIN {FS = "."} {print $$2}')

DEFINES += -DAMREX_NVCC_VERSION=$(nvcc_version)
DEFINES += -DAMREX_NVCC_MAJOR_VERSION=$(nvcc_major_version)
DEFINES += -DAMREX_NVCC_MINOR_VERSION=$(nvcc_minor_version)

# Disallow CUDA toolkit versions < 8.0.

nvcc_major_lt_8 = $(shell expr $(nvcc_major_version) \< 8)
ifeq ($(nvcc_major_lt_8),1)
  $(error Your nvcc version is $(nvcc_version). This is unsupported. Please use CUDA toolkit version 8.0 or newer.)
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
  ifeq ($(gcc_major_version),4)
    CXXFLAGS_FROM_HOST := -ccbin=g++ -Xcompiler='$(CXXFLAGS) --std=c++11' --std=c++11
  else
    CXXFLAGS_FROM_HOST := -ccbin=g++ -Xcompiler='$(CXXFLAGS) --std=c++14' --std=c++14
  endif
  CFLAGS_FROM_HOST := $(CXXFLAGS_FROM_HOST)
else ifeq ($(lowercase_nvcc_host_comp),pgi)
  CXXFLAGS_FROM_HOST := -ccbin=pgc++ -Xcompiler='$(CXXFLAGS)' --std=c++11
  CFLAGS_FROM_HOST := $(CXXFLAGS_FROM_HOST)
else
  CXXFLAGS_FROM_HOST := -ccbin=$(CXX) -Xcompiler='$(CXXFLAGS)'
  CFLAGS_FROM_HOST := $(CXXFLAGS_FROM_HOST)
endif

NVCC_FLAGS = -Wno-deprecated-gpu-targets -m64 -arch=compute_$(CUDA_ARCH) -code=sm_$(CUDA_ARCH) -maxrregcount=$(CUDA_MAXREGCOUNT)
# Unfortunately, on cori with cuda 10.0 this fails in thrust code
# NVCC_FLAGS += --Werror=cross-execution-space-call

ifeq ($(DEBUG),TRUE)
  NVCC_FLAGS += -g -G
else
  NVCC_FLAGS += -lineinfo --ptxas-options=-O3,-v
endif

ifneq ($(USE_CUDA_FAST_MATH),FALSE)
  NVCC_FLAGS += --use_fast_math
endif

CXXFLAGS = $(CXXFLAGS_FROM_HOST) $(NVCC_FLAGS) -dc -x cu
CFLAGS   =   $(CFLAGS_FROM_HOST) $(NVCC_FLAGS) -dc -x cu

ifeq ($(nvcc_version),9.2)
  # relaxed constexpr not supported
  ifeq ($(USE_EB),TRUE)
    $(error Cuda 9.2 is not supported with USE_EB=TRUE. Use 9.1 or 10.0 instead.)
  endif
endif

CXXFLAGS += --expt-relaxed-constexpr --expt-extended-lambda

CXX = nvcc
CC  = nvcc
