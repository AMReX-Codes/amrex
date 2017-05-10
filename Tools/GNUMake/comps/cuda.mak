#
# Setup for compiling the CUDA version of AMReX
# Assumes you have set USE_CUDA=TRUE, and have
# set the variables PGI_PATH to the root PGI
# directory and CUDA_PATH to the root CUDA directory.
#
CXX = nvcc
CC  = nvcc
FC  = pgfortran
F90 = pgfortran

# Note that for .c files we're directing nvcc to treat them
# as pure C files that have no CUDA markup. This is a workaround
# to deal with some of the .c files in the code that are called
# from Fortran, which are currently not called in CUDA code.
# nvcc by default uses C++ name mangling for files treated as
# CUDA code, which doesn't play nicely with the Fortran calling
# routines that expect them to use C naming conventions. A longer
# term solution could be to conditionally include extern "C"
# wrappers in those files if USE_CUDA=TRUE.

CXXFLAGS = -Wno-deprecated-gpu-targets -x cu --std=c++11 -ccbin=g++
CFLAGS   = -Wno-deprecated-gpu-targets -x c -ccbin=gcc
FFLAGS   =
F90FLAGS =

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -Xcompiler='-g -O0 -fno-inline -ggdb -Wall -Wno-sign-compare -ftrapv'
  CFLAGS   += -Xcompiler='-g -O0 -fno-inline -ggdb -Wall -Wno-sign-compare -ftrapv'

  FFLAGS   += -g -O0 -Ktrap=divz,inv -Mchkptr
  F90FLAGS += -g -O0 -Ktrap=divz,inv -Mchkptr

else

  CXXFLAGS += -Xcompiler='-g -O3'
  CFLAGS   += -Xcompiler='-g -O3'
  FFLAGS   += -fast
  F90FLAGS += -fast

endif

CXXFLAGS += 

########################################################################

F90FLAGS += -module $(fmoddir) -I$(fmoddir) -Mdclchk
FFLAGS   += -module $(fmoddir) -I$(fmoddir) -Mextend

########################################################################

GENERIC_COMP_FLAGS =

ifeq ($(USE_OMP),TRUE)
  GENERIC_COMP_FLAGS += -mp=nonuma -Minfo=mp
endif

ifeq ($(USE_ACC),TRUE)
  GENERIC_COMP_FLAGS += -acc -Minfo=acc -ta=nvidia -lcudart -mcmodel=medium
else
  GENERIC_COMP_FLAGS += -noacc
endif

CXXFLAGS += 
CFLAGS   += 
FFLAGS   += $(GENERIC_COMP_FLAGS)
F90FLAGS += $(GENERIC_COMP_FLAGS)

########################################################################

FFLAGS   += -Mcuda=cuda8.0 -Mnomain
F90FLAGS += -Mcuda=cuda8.0 -Mnomain

override XTRALIBS += -lstdc++
