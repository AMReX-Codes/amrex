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

ifneq ($(BL_NO_FORT),TRUE)

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

endif  # BL_NO_FORT

# =============================================================================================

ifeq ($(HIP_PLATFORM),clang)

  ifeq ($(DEBUG),TRUE)
    CXXFLAGS += -g -O0 -ftrapv
    CFLAGS   += -g -O0 -ftrapv

    FFLAGS   += -g -O0 -ggdb -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv
    F90FLAGS += -g -O0 -ggdb -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv

  else  # DEBUG=FALSE flags

    CXXFLAGS += -g -O3
    CFLAGS   += -g -O3
    FFLAGS   += -g -O3
    F90FLAGS += -g -O3

  endif

  CXXFLAGS += -Wno-pass-failed  # disable this warning

  ifeq ($(WARN_ALL),TRUE)
    warning_flags = -Wall -Wextra -Wno-sign-compare -Wunreachable-code -Wnull-dereference
    warning_flags += -Wfloat-conversion -Wextra-semi

    warning_flags += -Wpedantic

    ifneq ($(WARN_SHADOW),FALSE)
      warning_flags += -Wshadow
    endif

    CXXFLAGS += $(warning_flags) -Woverloaded-virtual
    CFLAGS += $(warning_flags)
  endif

#  ifeq ($(WARN_ERROR),TRUE)
#    CXXFLAGS += -Werror
#    CFLAGS += -Werror
#  endif

  # Generic HIP info
  ROC_PATH=$(realpath $(dir $(HIP_PATH)))
  INCLUDE_LOCATIONS += $(HIP_PATH)/include

  # rocRand
  INCLUDE_LOCATIONS += $(ROC_PATH)/rocrand/include $(ROC_PATH)/hiprand/include
  LIBRARY_LOCATIONS += $(ROC_PATH)/rocrand/lib $(ROC_PATH)/hiprand/lib
  LIBRARIES += -Wl,--rpath=$(ROC_PATH)/rocrand/lib -Wl,--rpath=$(ROC_PATH)/hiprand/lib -lhiprand -lrocrand 

  # rocPrim - Header only
  INCLUDE_LOCATIONS += $(ROC_PATH)/rocprim/include

  # rocThrust - Header only
  # INCLUDE_LOCATIONS += $(ROC_PATH)/rocthrust/include

# =============================================================================================

else ifeq ($(HIP_PLATFORM),nvcc)
  $(error HIP_PLATFORM nvcc is not supported at this time. Use USE_CUDA to compile for NVIDIA platforms.)
endif

# =============================================================================================
