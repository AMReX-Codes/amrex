GNU_DOT_MAK_INCLUDED = TRUE

########################################################################

ifndef AMREX_CCOMP
  AMREX_CCOMP = gnu
endif

ifndef AMREX_FCOMP
  AMREX_FCOMP = gnu
endif

########################################################################

ifeq ($(USE_CUDA),TRUE)
  ifdef NVCC_CCBIN
    GCC_VERSION_COMP = $(NVCC_CCBIN)
  else
    GCC_VERSION_COMP = g++
  endif
else
  GCC_VERSION_COMP = $(CXX)
endif

gcc_version       = $(shell $(GCC_VERSION_COMP) -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;')
gcc_major_version = $(shell $(GCC_VERSION_COMP) -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;' | sed -e 's;\..*;;')
gcc_minor_version = $(shell $(GCC_VERSION_COMP) -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;' | sed -e 's;[^.]*\.;;' | sed -e 's;\..*;;')

COMP_VERSION = $(gcc_version)

########################################################################

GENERIC_GNU_FLAGS =

ifeq ($(EXPORT_DYNAMIC),TRUE)
  CPPFLAGS += -DAMREX_EXPORT_DYNAMIC
  LIBRARIES += -ldl
  GENERIC_GNU_FLAGS += -rdynamic -fno-omit-frame-pointer
endif

gcc_major_ge_5 = $(shell expr $(gcc_major_version) \>= 5)
gcc_major_ge_6 = $(shell expr $(gcc_major_version) \>= 6)
gcc_major_ge_7 = $(shell expr $(gcc_major_version) \>= 7)
gcc_major_ge_8 = $(shell expr $(gcc_major_version) \>= 8)
gcc_major_ge_9 = $(shell expr $(gcc_major_version) \>= 9)
gcc_major_ge_10 = $(shell expr $(gcc_major_version) \>= 10)

ifeq ($(THREAD_SANITIZER),TRUE)
  GENERIC_GNU_FLAGS += -fsanitize=thread
endif
ifeq ($(FSANITIZER),TRUE)
  GENERIC_GNU_FLAGS += -fsanitize=address -fsanitize=undefined
  ifeq ($(gcc_major_ge_8),1)
    GENERIC_GNU_FLAGS += -fsanitize=pointer-compare -fsanitize=pointer-subtract
    GENERIC_GNU_FLAGS += -fsanitize=builtin -fsanitize=pointer-overflow
  endif
endif

ifeq ($(USE_OMP),TRUE)
  GENERIC_GNU_FLAGS += -fopenmp
endif

########################################################################
########################################################################
########################################################################

ifeq ($(AMREX_CCOMP),gnu)

CXX = g++
CC  = gcc

CXXFLAGS =
CFLAGS   =

########################################################################

CXXFLAGS += -Werror=return-type
CFLAGS   += -Werror=return-type

ifeq ($(DEBUG),TRUE)
  CXXFLAGS += -g -O0 -ggdb -ftrapv
  CFLAGS   += -g -O0 -ggdb -ftrapv
else
  CXXFLAGS += -g -O3
  CFLAGS   += -g -O3
endif

ifeq ($(WARN_ALL),TRUE)
  warning_flags = -Wall -Wextra -Wno-sign-compare -Wunreachable-code

  ifneq ($(USE_CUDA),TRUE)
    # With -Wpedantic I got 650 MB of warnings
    warning_flags += -Wpedantic
  endif

  ifeq ($(gcc_major_ge_6),1)
    warning_flags += -Wnull-dereference
  endif

  ifeq ($(gcc_major_ge_5),1)
    warning_flags += -Wfloat-conversion
  endif

  ifneq ($(WARN_SHADOW),FALSE)
    warning_flags += -Wshadow
  endif

  ifeq ($(gcc_major_version),7)
    warning_flags += -Wno-array-bounds
  endif

  CXXFLAGS += $(warning_flags) -Woverloaded-virtual
  CFLAGS += $(warning_flags)
endif

ifeq ($(WARN_ERROR),TRUE)
  CXXFLAGS += -Werror
  CFLAGS += -Werror
endif

ifeq ($(USE_GPROF),TRUE)
  CXXFLAGS += -pg
  CFLAGS += -pg
endif


ifeq ($(USE_COMPILE_PIC),TRUE)
  CXXFLAGS = -fPIC
  CFLAGS = -fPIC
endif

ifeq ($(ERROR_DEPRECATED),TRUE)
  CXXFLAGS += -Werror=deprecated
  CFLAGS += -Werror=deprecated
endif

########################################################################

ifdef CXXSTD
  CXXSTD := $(strip $(CXXSTD))
  ifeq ($(shell expr $(gcc_major_version) \< 5),1)
    ifeq ($(CXXSTD),c++14)
      $(error C++14 support requires GCC 5 or newer.)
    endif
  endif
  CXXFLAGS += -std=$(CXXSTD)
else
  ifeq ($(gcc_major_version),4)
    CXXFLAGS += -std=c++11
  else ifeq ($(gcc_major_version),5)
    CXXFLAGS += -std=c++14
  endif
endif

CFLAGS   += -std=gnu99

########################################################################

CXXFLAGS += $(GENERIC_GNU_FLAGS) -pthread
CFLAGS   += $(GENERIC_GNU_FLAGS)

endif # AMREX_CCOMP == gnu

########################################################################
########################################################################
########################################################################

ifeq ($(AMREX_FCOMP),gnu)

FC  = gfortran
F90 = gfortran

FFLAGS   =
F90FLAGS =

########################################################################

ifeq ($(DEBUG),TRUE)

  FFLAGS   += -g -O0 -ggdb -fcheck=bounds -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv
  F90FLAGS += -g -O0 -ggdb -fcheck=bounds -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv

else

  FFLAGS   += -g -O3
  F90FLAGS += -g -O3

endif

ifeq ($(USE_GPROF),TRUE)

  FFLAGS += -pg
  F90FLAGS += -pg

endif

ifeq ($(USE_COMPILE_PIC),TRUE)

  FFLAGS = -fPIC
  F90FLAGS = -fPIC

endif

########################################################################

FFLAGS   += -ffixed-line-length-none -fno-range-check -fno-second-underscore
F90FLAGS += -ffree-line-length-none -fno-range-check -fno-second-underscore -fimplicit-none

FMODULES =  -J$(fmoddir) -I $(fmoddir)

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

FFLAGS   += $(GENERIC_GNU_FLAGS)
F90FLAGS += $(GENERIC_GNU_FLAGS)

endif  # BL_NO_FORT

endif # AMREX_FCOMP == gnu
