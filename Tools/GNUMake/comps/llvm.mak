#
# Generic setup for using gcc
#
CXX = clang++
CC  = clang
FC  = gfortran
F90 = gfortran

CXXFLAGS =
CFLAGS   =
FFLAGS   =
F90FLAGS =

########################################################################

clang_version       = $(shell $(CXX) --version | head -1 | sed -e 's/.*version.*\([0-9]\+\.[0-9]\+\.[0-9]\+\).*/\1/')
clang_major_version = $(shell $(CXX) --version | head -1 | sed -e 's/.*version.*\([0-9]\+\.[0-9]\+\.[0-9]\+\).*/\1/' | sed -e 's;\..*;;')
clang_minor_version = $(shell $(CXX) --version | head -1 | sed -e 's/.*version.*\([0-9]\+\.[0-9]\+\.[0-9]\+\).*/\1/' | sed -e 's;[^.]*\.;;' | sed -e 's;\..*;;')

COMP_VERSION = $(clang_version)

########################################################################

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0 -ftrapv
  CFLAGS   += -g -O0 -ftrapv

  FFLAGS   += -g -O0 -ggdb -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv
  F90FLAGS += -g -O0 -ggdb -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv

else

  CXXFLAGS += -g -O3
  CFLAGS   += -g -O3
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
  CXXSTD := $(strip $(CXXSTD))
else
  CXXSTD := c++14
endif

CXXFLAGS += -std=$(CXXSTD)
CFLAGS   += -std=c99

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

ifeq ($(FSANITIZER),TRUE)
  override XTRALIBS += -lubsan
endif
