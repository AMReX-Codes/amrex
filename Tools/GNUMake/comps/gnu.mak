#
# Generic setup for using gcc
#
CXX = g++
CC  = gcc
FC  = gfortran
F90 = gfortran

CXXFLAGS =
CFLAGS   =
FFLAGS   =
F90FLAGS =

########################################################################

gcc_version       = $(shell $(CXX) -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;')
gcc_major_version = $(shell $(CXX) -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;' | sed -e 's;\..*;;')
gcc_minor_version = $(shell $(CXX) -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;' | sed -e 's;[^.]*\.;;' | sed -e 's;\..*;;')

COMP_VERSION = $(gcc_version)

DEFINES += -DBL_GCC_VERSION=$(gcc_version)
DEFINES += -DBL_GCC_MAJOR_VERSION=$(gcc_major_version)
DEFINES += -DBL_GCC_MINOR_VERSION=$(gcc_minor_version)

########################################################################

gcc_major_le_4 = $(shell expr $(gcc_major_version) \<= 4)
gcc_minor_lt_8 = $(shell expr $(gcc_minor_version) \< 8)
ifeq ($(gcc_major_le_4),1)
  ifeq ($(gcc_minor_lt_8),1)
    $(warning Your default GCC is version $(gcc_version). This might break during build. We therefore recommend that you specify a GCC >= 4.8 in your Make.local. The the docs on building AMReX for an example.)
  endif
endif

########################################################################

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0 -fno-inline -ggdb -Wshadow -Wall -Wno-sign-compare -ftrapv -Wno-unused-but-set-variable
  CFLAGS   += -g -O0 -fno-inline -ggdb -Wshadow -Wall -Wno-sign-compare -ftrapv

  FFLAGS   += -g -O0 -ggdb -fcheck=bounds -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv
  F90FLAGS += -g -O0 -ggdb -fcheck=bounds -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv

  ifneq ($(gcc_major_version),$(filter $(gcc_major_version),4 5))
    CXXFLAGS += -Wnull-dereference
    CFLAGS += -Wnull-dereference
  endif

else

  CXXFLAGS += -g -O3
  CFLAGS   += -g -O3
  FFLAGS   += -g -O3
  F90FLAGS += -g -O3

endif


ifeq ($(USE_GPROF),TRUE)

  CXXFLAGS += -pg
  CFLAGS += -pg
  FFLAGS += -pg
  F90FLAGS += -pg

endif

########################################################################

ifeq ($(gcc_major_version),4)
  CXXFLAGS += -std=c++11
else ifeq ($(gcc_major_version),5)
  CXXFLAGS += -std=c++14
endif
CFLAGS     += -std=gnu99

FFLAGS   += -ffixed-line-length-none -fno-range-check -fno-second-underscore -J$(fmoddir) -I $(fmoddir)
F90FLAGS += -ffree-line-length-none -fno-range-check -fno-second-underscore -J$(fmoddir) -I $(fmoddir) -fimplicit-none

########################################################################

GENERIC_GNU_FLAGS =

ifeq ($(THREAD_SANITIZER),TRUE)
  GENERIC_GNU_FLAGS += -fsanitize=thread
endif
ifeq ($(FSANITIZER),TRUE)
  GENERIC_GNU_FLAGS += -fsanitize=address -fsanitize=undefined
endif

ifeq ($(USE_OMP),TRUE)
  GENERIC_GNU_FLAGS += -fopenmp
endif

CXXFLAGS += $(GENERIC_GNU_FLAGS)
CFLAGS   += $(GENERIC_GNU_FLAGS)
FFLAGS   += $(GENERIC_GNU_FLAGS)
F90FLAGS += $(GENERIC_GNU_FLAGS)

########################################################################

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
