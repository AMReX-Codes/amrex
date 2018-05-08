#
# Generic setup for using PGI
#

# Due to how we are handling CUDA C compiling,
# we must handle all of the C/C++ options first,
# and only handle Fortran afterward.

CXX = pgc++
CC  = pgcc

########################################################################

pgi_version = $(shell $(CXX) -V 2>&1 | grep 'target')

COMP_VERSION = $(pgi_version)

########################################################################

CXXFLAGS =
CFLAGS   =

# 2017-06-23: PGI 17.4 causes the CUDA code to break when using -O2/-fast

PGI_OPT := -fast

ifeq ($(USE_CUDA),TRUE)
  PGI_OPT := -O
endif

ifeq ($(DEBUG),TRUE)

  # 2016-12-02: pgi 16.10 doesn't appear to like -traceback together with c++11

  CXXFLAGS += -g -O0 -Mbounds
  CFLAGS   += -g -O0 -Mbounds

else

  CXXFLAGS += -gopt $(PGI_OPT)
  CFLAGS   += -gopt $(PGI_OPT)

endif

CXXFLAGS += --c++11
CFLAGS   += -c99



GENERIC_PGI_FLAGS =

ifeq ($(USE_OMP),TRUE)
  GENERIC_PGI_FLAGS += -mp=nonuma -Minfo=mp
endif

ifeq ($(USE_ACC),TRUE)
  GENERIC_PGI_FLAGS += -acc -Minfo=acc -ta=nvidia -mcmodel=medium
else
  GENERIC_PGI_FLAGS += -noacc
endif

CXXFLAGS += $(GENERIC_PGI_FLAGS)
CFLAGS   += $(GENERIC_PGI_FLAGS)

ifeq ($(USE_CUDA),TRUE)
  include $(AMREX_HOME)/Tools/GNUMake/comps/nvcc.mak
endif

override XTRALIBS += -lstdc++



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

  FFLAGS   += -gopt $(PGI_OPT)
  F90FLAGS += -gopt $(PGI_OPT)

endif

# Note that we do not have a Fortran main

F90FLAGS += -Mnomain
FFLAGS   += -Mnomain

ifeq ($(USE_CUDA),TRUE)

  F90FLAGS += -Mcuda=$(CUDA_VERSION)
  FFLAGS   += -Mcuda=$(CUDA_VERSION)

  F90FLAGS += CUDA_HOME=$(COMPILE_CUDA_PATH)
  FFLAGS   += CUDA_HOME=$(COMPILE_CUDA_PATH)

  ifdef CUDA_MAXREGCOUNT
    F90FLAGS += -Mcuda=maxregcount:$(CUDA_MAXREGCOUNT)
    FFLAGS   += -Mcuda=maxregcount:$(CUDA_MAXREGCOUNT)
  endif

  LINK_WITH_FORTRAN_COMPILER = TRUE

endif


########################################################################

F90FLAGS += -module $(fmoddir) -I$(fmoddir) -Mdclchk
FFLAGS   += -module $(fmoddir) -I$(fmoddir) -Mextend

########################################################################

FFLAGS   += $(GENERIC_PGI_FLAGS)
F90FLAGS += $(GENERIC_PGI_FLAGS)

########################################################################

override XTRALIBS += -lstdc++ -pgf90libs -latomic

LINK_WITH_FORTRAN_COMPILER ?= $(USE_F_INTERFACES)
