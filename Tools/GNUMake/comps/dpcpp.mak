#
# Generic setup for using gcc
#
CXX = icpx
CC  = icx
FC  = ifx
F90 = ifx

CXXFLAGS =
CFLAGS   =
FFLAGS   =
F90FLAGS =

ifeq ($(USE_ONEDPL),TRUE)
# TBB and PSTL are broken in oneAPI 2021.3.0
# https://software.intel.com/content/www/us/en/develop/articles/intel-oneapi-dpcpp-library-release-notes.html#inpage-nav-2-3
  CPPFLAGS += -DAMREX_USE_ONEDPL -D_GLIBCXX_USE_TBB_PAR_BACKEND=0 -DPSTL_USE_PARALLEL_POLICIES=0
endif

########################################################################

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0 #-ftrapv
  CFLAGS   += -g -O0 #-ftrapv

  FFLAGS   += -g -O0 -traceback -check bounds,uninit,pointers
  F90FLAGS += -g -O0 -traceback -check bounds,uninit,pointers

else

  CXXFLAGS += -g1 -O3 # // xxxx SYCL: todo -g in beta6 causes a lot of warning messages
  CFLAGS   += -g1 -O3 #                       and makes linking much slower
  FFLAGS   += -g -O3
  F90FLAGS += -g -O3

endif

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

# disable warning: comparison with infinity always evaluates to false in fast floating point modes [-Wtautological-constant-compare]
#                  return std::isinf(m);
# appeared since 2021.4.0
CXXFLAGS += -Wno-tautological-constant-compare

ifeq ($(WARN_ERROR),TRUE)
  CXXFLAGS += -Werror
  CFLAGS += -Werror
endif

########################################################################

ifdef CXXSTD
  CXXFLAGS += -std=$(strip $(CXXSTD))
else
  CXXFLAGS += -std=c++17
endif

CXXFLAGS += -fsycl
CFLAGS   += -std=c11

ifneq ($(SYCL_SPLIT_KERNEL),FALSE)
  CXXFLAGS += -fsycl-device-code-split=per_kernel
endif

# temporary work-around for oneAPI beta08 bug
#   define "long double" as 64bit for C++ user-defined literals
#   https://github.com/intel/llvm/issues/2187
CXXFLAGS += -mlong-double-64 -Xclang -mlong-double-64

F90FLAGS += -implicitnone

FMODULES = -module $(fmoddir) -I$(fmoddir)

########################################################################

GENERIC_COMP_FLAGS =

ifeq ($(EXPORT_DYNAMIC),TRUE)
  CPPFLAGS += -DAMREX_EXPORT_DYNAMIC
  LIBRARIES += -Xlinker -export_dynamic
  GENERIC_COMP_FLAGS += -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer
endif

ifeq ($(THREAD_SANITIZER),TRUE)
  GENERIC_COMP_FLAGS += -fsanitize=thread
endif
ifeq ($(FSANITIZER),TRUE)
  GENERIC_COMP_FLAGS += -fsanitize=address -fsanitize=undefined
endif

CXXFLAGS += $(GENERIC_COMP_FLAGS) -pthread
CFLAGS   += $(GENERIC_COMP_FLAGS)

ifeq ($(USE_OMP),TRUE)
  CXXFLAGS +=
  CFLAGS   +=
  FFLAGS   += -qopenmp
  F90FLAGS += -qopenmp
endif

########################################################################

ifneq ($(BL_NO_FORT),TRUE)
  override XTRALIBS += -lifcore
  ifeq ($(USE_OMP),TRUE)
    override XTRALIBS += -lifcoremt
  endif
endif

LDFLAGS += -fsycl-device-lib=libc,libm-fp32,libm-fp64

ifeq ($(SYCL_AOT),TRUE)
  ifndef AMREX_INTEL_ARCH
    ifdef INTEL_ARCH
      AMREX_INTEL_ARCH = $(INTEL_ARCH)
    endif
  endif
  ifdef AMREX_INTEL_ARCH
    amrex_intel_gpu_target = $(AMREX_INTEL_ARCH)
  else
    # amrex_intel_gpu_target = *
    $(error Either INTEL_ARCH or AMREX_INTEL_ARCH must be specified when SYCL_AOT is TRUE.)
  endif
  CXXFLAGS += -fsycl-targets=spir64_gen
  amrex_sycl_backend_flags = -device $(amrex_intel_gpu_target)
  SYCL_AOT_GRF_MODE ?= Default
  ifneq ($(SYCL_AOT_GRF_MODE),Default)
    ifeq ($(SYCL_AOT_GRF_MODE),Large)
      amrex_sycl_backend_flags += -internal_options -ze-opt-large-register-file
    else ifeq ($(SYCL_AOT_GRF_MODE),AutoLarge)
      amrex_sycl_backend_flags += -options -ze-intel-enable-auto-large-GRF-mode
    else
      $(error SYCL_AOT_GRF_MODE ($(SYCL_AOT_GRF_MODE)) must be either Default, Large, or AutoLarge)
    endif
  endif
  LDFLAGS += -Xsycl-target-backend '$(amrex_sycl_backend_flags)'
endif

ifeq ($(DEBUG),TRUE)
  # This might be needed for linking device code larger than 2GB.
  LDFLAGS += -fsycl-link-huge-device-code
endif

ifeq ($(FSANITIZER),TRUE)
  override XTRALIBS += -lubsan
endif

AMREX_CCACHE_ENV = CCACHE_DEPEND=1
