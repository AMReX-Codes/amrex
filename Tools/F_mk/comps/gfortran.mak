# to compile mt19937ar.f90, we need -fno-range-check, since that
# routine relies on overflows when doing initializations

  FCOMP_VERSION = $(shell $(COMP) -v 2>&1 | grep 'version')

  FC  = $(COMP)
  F90 = $(COMP)
  ifdef CCOMP
    CC  = $(CCOMP)
    CXX = $(CCOMP)
  else
    CC  = gcc
    CXX = g++
  endif

  # Enforce gcc minimum version of 4.8.

  gcc_version       = $(shell $(CXX) -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;')
  gcc_major_version = $(shell $(CXX) -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;' | sed -e 's;\..*;;')
  gcc_minor_version = $(shell $(CXX) -dumpfullversion -dumpversion | head -1 | sed -e 's;.*  *;;' | sed -e 's;[^.]*\.;;' | sed -e 's;\..*;;')

  gcc_major_le_4 = $(shell expr $(gcc_major_version) \<= 4)
  gcc_minor_lt_8 = $(shell expr $(gcc_minor_version) \< 8)
  ifeq ($(gcc_major_le_4),1)
    ifeq ($(gcc_minor_lt_8),1)
      $(warning Your default GCC is version $(gcc_version). This might break during build. We therefore recommend that you specify a GCC >= 4.8 in your Make.local. The the docs on building AMReX for an example.)
    endif
  endif

  F90FLAGS += -J$(mdir) -I$(mdir)
  FFLAGS   += -J$(mdir) -I$(mdir)
  CFLAGS   += -std=c99 -Wall
  CXXFLAGS += -std=c++11 -Wall

  ifdef NDEBUG
    F90FLAGS += -g -O2 -ftree-vectorize -fno-range-check
    FFLAGS   += -g -O2 -ftree-vectorize -fno-range-check
    CFLAGS   += -g -O2 -ftree-vectorize
    CXXFLAGS += -g -O2 -ftree-vectorize
  else
    F90FLAGS += -g -fno-range-check -O1 -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid -finit-real=snan
    FFLAGS   += -g -fno-range-check -O1 -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid -finit-real=snan
    CFLAGS   += -g -O1
    CXXFLAGS += -g -O1
  endif

  ifdef OMP
    F90FLAGS += -fopenmp
    FFLAGS   += -fopenmp
    CFLAGS   += -fopenmp
    CXXFLAGS += -fopenmp
  endif

  ifdef GPROF
    F90FLAGS += -pg
    FFLAGS   += -pg
    CFLAGS   += -pg
    CXXFLAGS += -pg
  endif

  ifdef ROSE
    F90FLAGS += -ffree-line-length-none
    FFLAGS   += -ffree-line-length-none -ffixed-line-length-none

    ROSEVERBOSE ?= 0
    ROSEFLAGS := -rose:verbose $(ROSEVERBOSE)

    ifdef OMP
      ROSEFLAGS += -rose:openmp
    endif
  endif

  ifdef FSANITIZER
    F90FLAGS += -fsanitize=address -fsanitize=undefined
    FFLAGS   += -fsanitize=address -fsanitize=undefined
    CFLAGS   += -fsanitize=address -fsanitize=undefined
    CXXFLAGS += -fsanitize=address -fsanitize=undefined
  endif

  ifdef THREAD_SANITIZER
    F90FLAGS += -fsanitize=thread
    FFLAGS   += -fsanitize=thread
    CFLAGS   += -fsanitize=thread
    CXXFLAGS += -fsanitize=thread
  endif

  xtr_libraries += -lstdc++

  ifeq ($(ARCH), Darwin)
       xtr_libraries += -lc++
  endif
