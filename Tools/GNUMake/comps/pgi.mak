
ifndef AMREX_CCOMP
  AMREX_CCOMP = pgi
endif

ifndef AMREX_FCOMP
  AMREX_FCOMP = pgi
endif

########################################################################

pgi_version = $(shell $(CXX) -V 2>&1 | grep 'target' | sed 's|.*$(CXX) \([0-9\.]*\).*|\1|')
pgi_major_version = $(shell echo $(pgi_version) | cut -f1 -d.)
pgi_minor_version = $(shell echo $(pgi_version) | cut -f2 -d.)

gcc_version       = $(shell g++ -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;')
gcc_major_version = $(shell g++ -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;' | sed -e 's;\..*;;')
gcc_minor_version = $(shell g++ -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;' | sed -e 's;[^.]*\.;;' | sed -e 's;\..*;;')

COMP_VERSION = $(pgi_version)

########################################################################

GENERIC_PGI_FLAGS =

ifeq ($(USE_OMP),TRUE)
  GENERIC_PGI_FLAGS += -mp=nonuma -Minfo=mp
endif

ifeq ($(USE_ACC),TRUE)
  GENERIC_PGI_FLAGS += -acc -Minfo=accel -mcmodel=medium
  ifneq ($(CUDA_ARCH),)
    GENERIC_PGI_FLAGS += -ta=tesla:cc$(CUDA_ARCH)
  else
    GENERIC_PGI_FLAGS += -ta=tesla
  endif
else
  GENERIC_PGI_FLAGS += -noacc
endif

# Note that -O2 is the default optimization level for PGI

PGI_OPT := -O2 -fast

########################################################################
########################################################################
########################################################################

ifeq ($(AMREX_CCOMP),pgi)

CXX = pgc++
CC  = pgcc

########################################################################

CXXFLAGS =
CFLAGS   =

# Allow -gopt to be disabled to work around a compiler bug on P9.

PGI_GOPT ?= TRUE

ifeq ($(DEBUG),TRUE)

  # 2016-12-02: pgi 16.10 doesn't appear to like -traceback together with c++11

  CXXFLAGS += -g -O0 -Mbounds
  CFLAGS   += -g -O0 -Mbounds

else

  CXXFLAGS += $(PGI_OPT)
  CFLAGS   += $(PGI_OPT)

  ifeq ($(PGI_GOPT),TRUE)

    CXXFLAGS += -gopt
    CFLAGS   += -gopt

  endif

endif

# The logic here should be consistent with what's in nvcc.mak
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

CFLAGS   += -c99

CXXFLAGS += $(GENERIC_PGI_FLAGS)
CFLAGS   += $(GENERIC_PGI_FLAGS)

else # AMREX_CCOMP == pgi

# If we're using OpenACC but also CUDA, then nvcc will be the C++ compiler. If
# we want to call the OpenACC API from C++ then we need to make sure we have
# the includes for it, because PGI may not be the host compiler for nvcc.

ifeq ($(USE_ACC),TRUE)
  PGI_BIN_LOCATION := $(shell pgc++ -show 2>&1 | grep CPPCOMPDIR | awk '{print $$7}' | cut -c2-)
  PGI_LOCATION := $(shell dirname $(PGI_BIN_LOCATION))
  INCLUDE_LOCATIONS += $(PGI_LOCATION)/etc/include_acc
endif

endif # AMREX_CCOMP == pgi

########################################################################
########################################################################
########################################################################

ifeq ($(AMREX_FCOMP),pgi)

#
# Now set the Fortran flags. Since this is done after the GNU include
# in the CUDA version, all of the GNU specific options are overriden.
#

FC  = pgfortran
F90 = pgfortran

FFLAGS   =
F90FLAGS =

ifeq ($(DEBUG),TRUE)

  # 2016-12-02: pgi 16.10 doesn't appear to like -traceback together with c++11

  FFLAGS   += -g -O0 -Mbounds -Ktrap=divz,inv -Mchkptr
  F90FLAGS += -g -O0 -Mbounds -Ktrap=divz,inv -Mchkptr

else

  FFLAGS   += $(PGI_OPT)
  F90FLAGS += $(PGI_OPT)

  ifeq ($(PGI_GOPT),TRUE)

    FFLAGS   += -gopt
    F90FLAGS += -gopt

  endif

endif

# Note that we do not have a Fortran main

ifneq ($(USE_F_INTERFACES),TRUE)
  F90FLAGS += -Mnomain
  FFLAGS   += -Mnomain
endif

ifeq ($(USE_CUDA),TRUE)

  F90FLAGS += -Mcuda=cc$(CUDA_ARCH),fastmath,charstring
  FFLAGS   += -Mcuda=cc$(CUDA_ARCH),fastmath,charstring

  ifeq ($(DEBUG),TRUE)
    F90FLAGS += -Mcuda=debug
    FFLAGS   += -Mcuda=debug
  else
    F90FLAGS += -Mcuda=lineinfo
    FFLAGS   += -Mcuda=lineinfo
  endif

  ifeq ($(CUDA_VERBOSE),TRUE)
    F90FLAGS += -Mcuda=ptxinfo
    FFLAGS   += -Mcuda=ptxinfo
  endif

  F90FLAGS += CUDA_HOME=$(COMPILE_CUDA_PATH)
  FFLAGS   += CUDA_HOME=$(COMPILE_CUDA_PATH)

  ifdef CUDA_MAXREGCOUNT
    F90FLAGS += -Mcuda=maxregcount:$(CUDA_MAXREGCOUNT)
    FFLAGS   += -Mcuda=maxregcount:$(CUDA_MAXREGCOUNT)
  endif

  DEFINES += -DAMREX_USE_CUDA_FORTRAN

  LINK_WITH_FORTRAN_COMPILER = TRUE

  ifeq ($(USE_MPI),TRUE)
  ifneq ($(findstring Open MPI, $(shell mpicxx -showme:version 2>&1)),)
    OMPI_FCFLAGS_ORIG = $(shell mpif90 -showme:compile)
    export OMPI_FCFLAGS := $(subst -pthread,-lpthread,$(OMPI_FCFLAGS_ORIG))
  endif
  endif

endif


########################################################################

F90FLAGS += -Mdclchk
FFLAGS   += -Mextend

FMODULES = -module $(fmoddir) -I$(fmoddir)

########################################################################

FFLAGS   += $(GENERIC_PGI_FLAGS)
F90FLAGS += $(GENERIC_PGI_FLAGS)

########################################################################

override XTRALIBS += -lstdc++ -pgf90libs -latomic

LINK_WITH_FORTRAN_COMPILER ?= $(USE_F_INTERFACES)

endif # AMREX_FCOMP == pgi
