#
# Generic setup for using gcc
#
CXX = dpcpp
CC  = dpcpp
FC  = gfortran
F90 = gfortran

CXXFLAGS =
CFLAGS   =
FFLAGS   =
F90FLAGS =

########################################################################

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0 #-ftrapv
  CFLAGS   += -g -O0 #-ftrapv

  FFLAGS   += -g -O0 -ggdb -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 #-ftrapv
  F90FLAGS += -g -O0 -ggdb -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 #-ftrapv

else

  CXXFLAGS += -O3 # // xxxx DPCPP: todo -g in beta6 causes a lot of warning messages
  CFLAGS   += -O3 #                       and makes linking much slower
#  CXXFLAGS += -g -O3
#  CFLAGS   += -g -O3
  FFLAGS   += -g -O3
  F90FLAGS += -g -O3

endif

CXXFLAGS += -Wno-pass-failed  # disable this warning

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

FFLAGS   += -ffixed-line-length-none -fno-range-check -fno-second-underscore
F90FLAGS += -ffree-line-length-none -fno-range-check -fno-second-underscore -fimplicit-none

FMODULES =  -J$(fmoddir) -I $(fmoddir)

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

ifeq ($(USE_OMP),TRUE)
  GENERIC_COMP_FLAGS += -fopenmp
endif

CXXFLAGS += $(GENERIC_COMP_FLAGS) -pthread
CFLAGS   += $(GENERIC_COMP_FLAGS)
FFLAGS   += $(GENERIC_COMP_FLAGS)
F90FLAGS += $(GENERIC_COMP_FLAGS)

########################################################################

ifneq ($(BL_NO_FORT),TRUE)

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

endif

DPCPP_DIR = $(shell dpcpp --version | tail -1 | sed -e 's/InstalledDir: //' | sed -e 's/linux\/bin/linux/')
ifeq ($(findstring beta09, $(DPCPP_DIR)), beta09)
  override XTRAOBJS += $(DPCPP_DIR)/lib/libsycl-glibc.o
  LDFLAGS += -device-math-lib=fp32,fp64
else
  LDFLAGS += -fsycl-device-lib=libc,libm-fp32,libm-fp64
endif

ifeq ($(FSANITIZER),TRUE)
  override XTRALIBS += -lubsan
endif
