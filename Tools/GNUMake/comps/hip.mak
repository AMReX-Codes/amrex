# Setup for HIP, using hipcc (HCC and clang will use the same compiler name).

HIP_PATH=$(shell hipconfig --path)
ifeq ($(HIP_PATH),)
  $(error hipconfig failed. Is the HIP toolkit available?)
endif

CXX = $(HIP_PATH)/bin/hipcc
CC  = $(HIP_PATH)/bin/hipcc
FC = gfortran
F90 = gfortran

ifdef CXXSTD
  CXXSTD := $(strip $(CXXSTD))
else
  CXXSTD := c++14
endif

#if less than a given version, throw error.

# Generic flags, always used
CXXFLAGS := -std=$(CXXSTD) -m64
CFLAGS   := -std=c99 -m64

FFLAGS   := -ffixed-line-length-none -fno-range-check -fno-second-underscore
F90FLAGS := -ffree-line-length-none -fno-range-check -fno-second-underscore -fimplicit-none

FMODULES =  -J$(fmoddir) -I $(fmoddir)

# =============================================================================================

# Taken straight from gnu 
# ask gfortran the name of the library to link in.  First check for the
# static version.  If it returns only the name w/o a path, then it
# was not found.  In that case, ask for the shared-object version.
gfortran_liba  = $(shell $(F90) -print-file-name=libgfortran.a)
gfortran_libso = $(shell $(F90) -print-file-name=libgfortran.so)

ifneq ($(gfortran_liba),libgfortran.a)  # if found the full path is printed, thus `neq`.
  LIBRARY_LOCATIONS += $(dir $(gfortran_liba))
else
  LIBRARY_LOCATIONS += $(dir $(gfortran_libso))
endif

override XTRALIBS += -lgfortran -lquadmath

# =============================================================================================

# This is designed only for dogora for now.
ifeq ($(HIP_PLATFORM),hcc)

  ifeq ($(DEBUG),TRUE)
    # From llvm
    CXXFLAGS += -g
    CFLAGS   += -g 
    FFLAGS   += -g -O0 -ggdb -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv
    F90FLAGS += -g -O0 -ggdb -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv

  else  # DEBUG=FALSE flags
  endif

  # Generic HIP info
  ROC_PATH=/opt/rocm
  INCLUDE_LOCATIONS += $(HIP_PATH)/include

  # rocRand
  INCLUDE_LOCATIONS += $(ROC_PATH)/rocrand/include $(ROC_PATH)/hiprand/include
  LIBRARY_LOCATIONS += $(ROC_PATH)/rocrand/lib $(ROC_PATH)/hiprand/lib
  LIBRARIES += -Wl,--rpath=$(ROC_PATH)/rocrand/lib -Wl,--rpath=$(ROC_PATH)/hiprand/lib -lhiprand -lrocrand 

  # rocPrim - Header only
  INCLUDE_LOCATIONS += $(ROC_PATH)/rocprim/include

  # rocThrust - Header only
  INCLUDE_LOCATIONS += $(ROC_PATH)/rocthrust/include

# =============================================================================================

# This is Summit. Likely broken.
else ifeq ($(HIP_PLATFORM),nvcc)
  $(error HIP_PLATFORM nvcc is not supported at this time. Use USE_CUDA to compile for NVIDIA platforms.)
#
#  CXXFLAGS_FROM_HOST := -ccbin=$(CXX) --std=c++14
#  CFLAGS_FROM_HOST := -ccbin=$(CXX)
#  HIPCC_FLAGS = -Wno-deprecated-gpu-targets -m64 -arch=compute_$(CUDA_ARCH) -code=sm_$(CUDA_ARCH) -maxrregcount=$(CUDA_MAXREGCOUNT)
#
#  ifeq ($(DEBUG),TRUE)
#    HIPCC_FLAGS += -g -G
#  else
#    HIPCC_FLAGS += -lineinfo --ptxas-options=-O3,-v
#  endif
#
#  ifneq ($(USE_CUDA_FAST_MATH),FALSE)
#    HIPCC_FLAGS += --use_fast_math
#  endif
#
#  CXXFLAGS = $(CXXFLAGS_FROM_HOST) $(HIPCC_FLAGS) -c -dc
#  CFLAGS   =   $(CFLAGS_FROM_HOST) $(HIPCC_FLAGS) -dc
#
#  CXXFLAGS += --expt-relaxed-constexpr --expt-extended-lambda
#
endif

# =============================================================================================
