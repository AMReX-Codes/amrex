PBOXLIB_HOME := ../..
FBOXLIB_HOME = ../../../fParallel
FBOXLIB_HOME =

PRECISION = DOUBLE
PROFILE = TRUE
PROFILE = FALSE
DEBUG	= TRUE
DEBUG	= FALSE
DIM	= 2
DIM	= 3

COMP=CC
FCOMP=ftn

USE_MPI=FALSE
USE_MPI=TRUE
USE_THREADS=TRUE
USE_THREADS=FALSE

USE_ARRAYVIEW = TRUE
USE_ARRAYVIEW = FALSE

#USE_XT3 = TRUE

EBASE = main

include $(PBOXLIB_HOME)/mk/Make.defs

# multigrid test sources
CEXE_sources += $(EBASE).cpp

include Make.package
include $(PBOXLIB_HOME)/bndrylib/Make.package
include $(PBOXLIB_HOME)/mglib/Make.package
include $(PBOXLIB_HOME)/BoxLib/Make.package

INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/bndrylib
vpathdir += $(PBOXLIB_HOME)/bndrylib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/BoxLib
vpathdir += $(PBOXLIB_HOME)/BoxLib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/iamrlib
vpathdir += $(PBOXLIB_HOME)/iamrlib
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/mglib
vpathdir += $(PBOXLIB_HOME)/mglib

ifeq ($(USE_XT3),TRUE)
  LIBRARIES += -lpgf90 -lpgf90_rpm1 -lpgf902 -lpgf90rtl -lpgftnrtl -lpghpf -lrt
  DEFINES += -DBL_PGI -DMPICH_SKIP_MPICXX -DBL_XT3
ifeq ($(USE_MPI),TRUE)
  CXX := CC -target=catamount
  FC  := ftn -target=catamount
  fC  := ftn -target=catamount
else
  CXX := pgCC
  FC  := pgf90
  fC  := pgf90
endif
  CXXOPTF += -O -Minline=levels:10 --no_exceptions
  FOPTF   += -O -fastsse -Minline=levels:10
  fOPTF   += -O -fastsse -Minline=levels:10
endif

ifdef FBOXLIB_HOME
  DEFINES += -DMG_USE_FBOXLIB

  include $(PBOXLIB_HOME)/MGT_Solver/Make.package
  INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/MGT_Solver
  vpathdir += $(PBOXLIB_HOME)/MGT_Solver

  include $(FBOXLIB_HOME)/mg/FParallelMG.mak
  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/mg/
  Fdirs   := boxlib mg extern/SPARSKIT extern/LAPACK
  Flocs   := $(foreach dir, $(Fdirs), $(FBOXLIB_HOME)/$(dir))
  vpathdir += $(Flocs)
endif

# We need these from amrlib 

CEXE_sources += BCRec.cpp 
CEXE_headers += BCRec.H BC_TYPES.H
INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/amrlib
vpathdir += $(PBOXLIB_HOME)/amrlib

# Needed from AmrVis
#CEXE_headers += AmrData.H
#INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/pAmrvis
#vpathdir += $(PBOXLIB_HOME)/pAmrvis

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(PBOXLIB_HOME)/mk/Make.rules
