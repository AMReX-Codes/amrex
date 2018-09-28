#
# If we are using CUDA, pull in the gcc compiler first
# and override it as necessary. This is done because the
# nvcc compiler driver does not always accept pgc++
# as a host compiler at present. However, if we're using
# OpenACC, then PGI is required, so we cannot do this.
#

ifneq ($(USE_ACC),TRUE)
  include $(AMREX_HOME)/Tools/GNUMake/comps/gnu.mak
endif

# Force immediate expansion of the GCC defines,
# since after this point GCC will no longer be
# the actual compiler defined in CXX.

DEFINES := $(DEFINES)

HOST_CXXFLAGS := $(CXXFLAGS)
HOST_CFLAGS   := $(CFLAGS)

HOST_CXX := $(CXX)
HOST_CC := $(CC)

CXXFLAGS_FROM_HOST := -ccbin=$(CXX) -Xcompiler='$(CXXFLAGS) --std=c++11'
CFLAGS_FROM_HOST := -ccbin=$(CC) -Xcompiler='$(CFLAGS)'

CXXFLAGS = $(CXXFLAGS_FROM_HOST) --std=c++11 -Wno-deprecated-gpu-targets -m64 -dc -x cu -arch=compute_$(CUDA_ARCH) -code=sm_$(CUDA_ARCH)
CFLAGS = $(CFLAGS_FROM_HOST) -Wno-deprecated-gpu-targets -m64 -dc -x c -arch=compute_$(CUDA_ARCH) -code=sm_$(CUDA_ARCH)

ifeq ($(DEBUG),TRUE)
  CXXFLAGS += -g -G
  CFLAGS += -g -G
else
  CXXFLAGS += --use_fast_math -lineinfo
  CFLAGS += --use_fast_math -lineinfo
endif

CXX := nvcc
CC := nvcc

override XTRALIBS :=

# Store the CUDA toolkit version.

nvcc_version       := $(shell $(CXX) --version | tail -1 | awk 'BEGIN {FS = ","} {print $$2}' | awk '{print $$2}')
nvcc_major_version := $(shell $(CXX) --version | tail -1 | awk 'BEGIN {FS = ","} {print $$2}' | awk '{print $$2}' | awk 'BEGIN {FS = "."} {print $$1}')
nvcc_minor_version := $(shell $(CXX) --version | tail -1 | awk 'BEGIN {FS = ","} {print $$2}' | awk '{print $$2}' | awk 'BEGIN {FS = "."} {print $$2}')

DEFINES += -DAMREX_NVCC_VERSION=$(nvcc_version)
DEFINES += -DAMREX_NVCC_MAJOR_VERSION=$(nvcc_major_version)
DEFINES += -DAMREX_NVCC_MINOR_VERSION=$(nvcc_minor_version)

# Disallow CUDA toolkit versions < 8.0.

nvcc_major_lt_8 = $(shell expr $(nvcc_major_version) \< 8)
ifeq ($(nvcc_major_lt_8),1)
  $(error Your nvcc version is $(nvcc_version). This is unsupported. Please use CUDA toolkit version 8.0 or newer.)
endif
