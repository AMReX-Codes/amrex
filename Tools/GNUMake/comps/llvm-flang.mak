#
# Setup for using clang/flang
#
CXX = clang++
CC  = clang
FC  = flang
F90 = flang

CXXFLAGS =
CFLAGS   =
FFLAGS   =
F90FLAGS =

########################################################################

clang_version       = $(shell $(CXX) --version | head -1 | sed -e 's/.*version.*\([0-9]\+\.[0-9]\+\.[0-9]\+\).*/\1/')
clang_major_version = $(shell $(CXX) --version | head -1 | sed -e 's/.*version.*\([0-9]\+\.[0-9]\+\.[0-9]\+\).*/\1/' | sed -e 's;\..*;;')
clang_minor_version = $(shell $(CXX) --version | head -1 | sed -e 's/.*version.*\([0-9]\+\.[0-9]\+\.[0-9]\+\).*/\1/' | sed -e 's;[^.]*\.;;' | sed -e 's;\..*;;')

COMP_VERSION = $(clang_version)

DEFINES += -DBL_CLANG_VERSION='$(clang_version)'
DEFINES += -DBL_CLANG_MAJOR_VERSION='$(clang_major_version)'
DEFINES += -DBL_CLANG_MINOR_VERSION='$(clang_minor_version)'

########################################################################

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0 -Wall -Wextra -Wno-sign-compare -Wno-unused-parameter -Wno-unused-variable -ftrapv
  CFLAGS   += -g -O0 -Wall -Wextra -Wno-sign-compare -Wno-unused-parameter -Wno-unused-variable -ftrapv

  FFLAGS   += -g -O0 -ggdb -Wuninitialized -Wunused -ftrapv
  F90FLAGS += -g -O0 -ggdb -Wuninitialized -Wunused -ftrapv

else

  CXXFLAGS += -g -O3
  CFLAGS   += -g -O3
  FFLAGS   += -g -O3
  F90FLAGS += -g -O3

endif

########################################################################

CXXFLAGS += -std=c++11
CFLAGS   += -std=c99

FMODULES = -J$(fmoddir) -I $(fmoddir)

########################################################################

GENERIC_COMP_FLAGS =

ifeq ($(THREAD_SANITIZER),TRUE)
  GENERIC_COMP_FLAGS += -fsanitize=thread
endif
ifeq ($(FSANITIZER),TRUE)
  GENERIC_COMP_FLAGS += -fsanitize=address -fsanitize=undefined
endif

ifeq ($(USE_OMP),TRUE)
  GENERIC_COMP_FLAGS += -fopenmp
endif

CXXFLAGS += $(GENERIC_COMP_FLAGS)
CFLAGS   += $(GENERIC_COMP_FLAGS)
FFLAGS   += $(GENERIC_COMP_FLAGS)
F90FLAGS += $(GENERIC_COMP_FLAGS)

########################################################################

# libflangrti.so is needed when using OpenMP. It is also needed in
# Castro_util.o to provide the symbol "__mth_i_idnint" in non-OpenMP builds
override XTRALIBS += -lflangrti -lflang -lpgmath

ifeq ($(FSANITIZER),TRUE)
  override XTRALIBS += -lubsan
endif

# Follow the method in pgi.mak and link with the Fortran compiler
# when we are using a Fortran main.

override XTRALIBS += -lstdc++
LINK_WITH_FORTRAN_COMPILER ?= $(USE_F_INTERFACES)
