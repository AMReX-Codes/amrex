
# hipcc_version       := $(shell hipcc --version | tail -1 | awk 'BEGIN {FS = ","} {print $$2}' | awk '{print $$2}')
# hipcc_major_version := $(shell hipcc --version | tail -1 | awk 'BEGIN {FS = ","} {print $$2}' | awk '{print $$2}' | awk 'BEGIN {FS = "."} {print $$1}')
# hipcc_minor_version := $(shell hipcc --version | tail -1 | awk 'BEGIN {FS = ","} {print $$2}' | awk '{print $$2}' | awk 'BEGIN {FS = "."} {print $$2}')
# 
# DEFINES += -DAMREX_HIPCC_VERSION=$(hipcc_version)
# DEFINES += -DAMREX_HIPCC_MAJOR_VERSION=$(hipcc_major_version)
# DEFINES += -DAMREX_HIPCC_MINOR_VERSION=$(hipcc_minor_version)

HIPCC_HOST_COMP ?= $(AMREX_CCOMP)

lowercase_hipcc_host_comp = $(shell echo $(HIPCC_HOST_COMP) | tr A-Z a-z)

ifeq ($(lowercase_hipcc_host_comp),$(filter $(lowercase_hipcc_host_comp),gcc gnu g++))
  lowercase_hipcc_host_comp = gnu
  AMREX_CCOMP = gnu
  ifndef GNU_DOT_MAK_INCLUDED
    include $(AMREX_HOME)/Tools/GNUMake/comps/gnu.mak
  endif
endif

#ifeq ($(lowercase_hipcc_host_comp),gnu)
#  ifeq ($(gcc_major_version),4)
#    CXXFLAGS_FROM_HOST := -ccbin=g++ -Xcompiler='$(CXXFLAGS) --std=c++11' --std=c++11
#  else
#    CXXFLAGS_FROM_HOST := -ccbin=g++ -Xcompiler='$(CXXFLAGS) --std=c++14' --std=c++14
#  endif
#  CFLAGS_FROM_HOST := $(CXXFLAGS_FROM_HOST)
#else ifeq ($(lowercase_hipcc_host_comp),pgi)
#  CXXFLAGS_FROM_HOST := -ccbin=pgc++ -Xcompiler='$(CXXFLAGS)' --std=c++11
#  CFLAGS_FROM_HOST := $(CXXFLAGS_FROM_HOST)
#else
#  CXXFLAGS_FROM_HOST := -ccbin=$(CXX) -Xcompiler='$(CXXFLAGS)'
#  CFLAGS_FROM_HOST := $(CXXFLAGS_FROM_HOST)
#endif

CXXFLAGS_FROM_HOST := -ccbin=$(CXX) --std=c++14
CFLAGS_FROM_HOST := -ccbin=$(CXX)

HIPCC_FLAGS = -Wno-deprecated-gpu-targets -m64 -arch=compute_$(CUDA_ARCH) -code=sm_$(CUDA_ARCH) -maxrregcount=$(CUDA_MAXREGCOUNT)

ifeq ($(DEBUG),TRUE)
  HIPCC_FLAGS += -g -G
else
  HIPCC_FLAGS += -lineinfo --ptxas-options=-O3,-v
endif

ifneq ($(USE_CUDA_FAST_MATH),FALSE)
  HIPCC_FLAGS += --use_fast_math
endif

CXXFLAGS = $(CXXFLAGS_FROM_HOST) $(HIPCC_FLAGS) -dc
CFLAGS   =   $(CFLAGS_FROM_HOST) $(HIPCC_FLAGS) -dc

CXXFLAGS += --expt-relaxed-constexpr --expt-extended-lambda

CXX = hipcc
CC  = hipcc
