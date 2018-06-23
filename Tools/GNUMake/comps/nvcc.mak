#
# If we are using CUDA, pull in the gcc compiler first
# and override it as necessary. This is done because the
# nvcc compiler driver does not work optimally with pgc++
# as a host compiler at present.
#

include $(AMREX_HOME)/Tools/GNUMake/comps/gnu.mak

# Force immediate expansion of the GCC defines,
# since after this point GCC will no longer be
# the actual compiler defined in CXX.

DEFINES := $(DEFINES)

CXXFLAGS := -Wno-deprecated-gpu-targets -dc -x cu --std=c++11 -ccbin=$(CXX) -Xcompiler='$(CXXFLAGS)'
CFLAGS := -Wno-deprecated-gpu-targets -dc -x c -ccbin=$(CC) -Xcompiler='$(CFLAGS)'

HOST_CXX := $(CXX)
HOST_CC := $(CC)

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


