
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

# =============================================================================================

ifeq ($(HIP_PLATFORM),hcc)
  CXXFLAGS_FROM_HOST := 
  CFLAGS_FROM_HOST := 
  HIPCC_FLAGS = -m64
  ROC_PATH=/opt/rocm
  HIP_PATH=$(shell hipconfig --path)

  ifeq ($(DEBUG),TRUE)
    HIPCC_FLAGS += -g
  
    # From llvm
    FFLAGS   += -g -O0 -ggdb -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv
    F90FLAGS += -g -O0 -ggdb -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv

  else  # DEBUG=FALSE flags
  endif

  # Generic HIP
  INCLUDE_LOCATIONS += $(HIP_PATH)/include

  # HIP rand
  INCLUDE_LOCATIONS += $(ROC_PATH)/rocrand/include $(ROC_PATH)/hiprand/include
  LIBRARY_LOCATIONS += -Wl,rpath,$(ROC_PATH)/rocrand/lib -Wl,rpath,$(ROC_PATH)/hiprand/lib
  LIBRARIES += -lhiprand -lrocrand 

  FC = gfortran
  F90 = gfortran

  CXXFLAGS = $(CXXFLAGS_FROM_HOST) $(HIPCC_FLAGS)
  CFLAGS   =   $(CFLAGS_FROM_HOST) $(HIPCC_FLAGS)

# =============================================================================================

else ifeq ($(HIP_PLATFORM),nvcc)

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

  CXXFLAGS = $(CXXFLAGS_FROM_HOST) $(HIPCC_FLAGS) -c -dc
  CFLAGS   =   $(CFLAGS_FROM_HOST) $(HIPCC_FLAGS) -dc

  CXXFLAGS += --expt-relaxed-constexpr --expt-extended-lambda

endif

# =============================================================================================

CXXFLAGS += -std=c++14
CFLAGS   += -std=c99

FFLAGS   += -ffixed-line-length-none -fno-range-check -fno-second-underscore
F90FLAGS += -ffree-line-length-none -fno-range-check -fno-second-underscore -fimplicit-none

FMODULES =  -J$(fmoddir) -I $(fmoddir)

# =============================================================================================

# Taken straight from llvm
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


CXX = /opt/rocm/hip/bin/hipcc
CC  = /opt/rocm/hip/bin/hipcc
