#
# Generic setup for using gcc
#
CXX = dpcpp
CC  = dpcpp
FC  = ifx
F90 = ifx

CXXFLAGS =
CFLAGS   =
FFLAGS   =
F90FLAGS =

########################################################################

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0 #-ftrapv
  CFLAGS   += -g -O0 #-ftrapv

  FFLAGS   += -g -O0 -traceback -check bounds,uninit,pointers
  F90FLAGS += -g -O0 -traceback -check bounds,uninit,pointers

else

  CXXFLAGS += -O3 # // xxxx DPCPP: todo -g in beta6 causes a lot of warning messages
  CFLAGS   += -O3 #                       and makes linking much slower
#  CXXFLAGS += -g -O3
#  CFLAGS   += -g -O3
  FFLAGS   += -g -O3
  F90FLAGS += -g -O3

endif

CXXFLAGS += -Wno-pass-failed # disable this warning

ifeq ($(WARN_ALL),TRUE)
  warning_flags = -Wall -Wextra -Wno-sign-compare -Wunreachable-code -Wnull-dereference
  warning_flags += -Wfloat-conversion -Wextra-semi

  ifneq ($(USE_CUDA),TRUE)
    warning_flags += -Wpedantic
  endif

  ifneq ($(WARN_SHADOW),FALSE)
    warning_flags += -Wshadow
  endif

  CXXFLAGS += $(warning_flags) -Woverloaded-virtual
  CFLAGS += $(warning_flags)
endif

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

CXXFLAGS += -Wno-error=sycl-strict -fsycl
CFLAGS   += -std=c99

ifneq ($(DEBUG),TRUE)  # There is currently a bug that DEBUG build will crash.
ifeq ($(DPCPP_AOT),TRUE)
  INTEL_CPU_LONG_NAME = $(shell cat /sys/devices/cpu/caps/pmu_name)
  ifneq ($(INTEL_CPU_LONG_NAME),)
    ifeq ($(INTEL_CPU_LONG_NAME),skylake)
      INTEL_CPU_SHORT_NAME = skl
    else ifeq ($(INTEL_CPU_LONG_NAME),kabylake)
      INTEL_CPU_SHORT_NAME = kbl
    else ifeq ($(INTEL_CPU_LONG_NAME),cascadelake)
      INTEL_CPU_SHORT_NAME = cfl
    else
      $(error AOT TODO: $(INTEL_CPU_LONG_NAME))
    endif
    CXXFLAGS += -fsycl-targets=spir64_gen-unknown-unknown-sycldevice -Xsycl-target-backend '-device $(INTEL_CPU_SHORT_NAME)'
  endif
endif
endif

ifneq ($(DPCPP_AOT),TRUE)
ifneq ($(DPCPP_SPLIT_KERNEL),FALSE)
  CXXFLAGS += -fsycl-device-code-split=per_kernel
endif
endif

# temporary work-around for DPC++ beta08 bug
#   define "long double" as 64bit for C++ user-defined literals
#   https://github.com/intel/llvm/issues/2187
CXXFLAGS += -mlong-double-64 -Xclang -mlong-double-64

# Beta09 has enabled early optimizations by default.  But this causes many
# tests to crash.  So we disable it.
CXXFLAGS += -fno-sycl-early-optimizations

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

ifeq ($(FSANITIZER),TRUE)
  override XTRALIBS += -lubsan
endif
