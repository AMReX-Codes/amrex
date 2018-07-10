# Makefile for AmrDeriveSpectrum using AMReX tools and FFTW.
#
# AmrDeriveSpectrum has the following requirements:
#
# - MPI and its headers should be installed. Tested to work with mpich.
#   The path to the MPI library is supplied via MPI_HOME.
#
# - FFTW2 and its headers should be installed.
#   FFTW2 should have been compiled with MPI support with double precision.
#   The path to the FFTW2 installation directory is supplied via FFTW2_HOME.
#
# - The amrex source code should be available, supplied via AMREX_HOME.
#

#############################################

AMREX_HOME = ../../..
MPI_HOME = /usr/lib/mpich
USRLIB = /usr/lib

#############################################

BOXLIB_HOME = $(AMREX_HOME)

BL_NOFAST=TRUE

#
# Variables for the user to set ...
#
########## define EBASE and BL_HAS_FORT here

EBASE = AmrDeriveSpectrum

#BL_HAS_FORT = TRUE

DEBUG	      = FALSE
PROFILE       = FALSE
USE_THREADS   = TRUE
DIM           = 3
USE_MPI       = TRUE
AIXLARGEMEM   = FALSE
COMP          = g++
FCOMP         = gfortran
USE_ARRAYVIEW = FALSE
BUILDCHEM     = FALSE

HERE = .

ifndef MPI_HOME
 $(error MPI_HOME is not defined.)
endif

ifndef AMREX_HOME
 $(error AMREX_HOME is not defined.)
endif


ifeq ($(EBASE), AmrDeriveSpectrum)
 INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Amr
 INCLUDE_LOCATIONS += $(MPI_HOME)
 LIBRARY_LOCATIONS += $(USRLIB)/lib
 include $(BOXLIB_HOME)/Src/Boundary/Make.package
 include $(BOXLIB_HOME)/Src/Extern/amrdata/Make.package
 LIBRARIES += -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
endif

include $(BOXLIB_HOME)/Tools/GNUMake/Make.defs

#DEFINES           += -DBL_DERIVE_IAMR
#DEFINES           += -DBL_DERIVE_HCHEM
DEFINES           += -DBL_NOLINEVALUES
#DEFINES           += -static

###### arrayview
USE_ARRAYVIEW = FALSE
ifeq (USE_ARRAYVIEW, TRUE)
DEFINES += -DBL_USE_ARRAYVIEW
ARRAYVIEWDIR = .
INCLUDE_LOCATIONS += $(ARRAYVIEWDIR)
endif

CEXE_sources += $(EBASE).cpp

include ./Make.package

DEFINES += -DBL_PARALLEL_IO
ifeq ($(COMP),KCC)
  ifeq ($(KCC_VERSION),3.3)
    CXXFLAGS += --diag_suppress 837
    CXXOPTF  += -Olimit 2400
    #CXXFLAGS += -fprm n -fptm u
  endif
endif

include $(BOXLIB_HOME)/Src/Base/Make.package

VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/Extern/amrdata
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/Extern/ProfParser

INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Base
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Amr
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Extern/amrdata
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Extern/ProfParser
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Boundary

ifeq ($(USE_ARRAYVIEW),TRUE)
  DEFINES += -DBL_USE_ARRAYVIEW
  DEFINES += -DBL_ARRAYVIEW_TAGBOX
endif

ifeq ($(MACHINE),OSF1)
#
# Some additional stuff for our preferred development/debugging environment.
#
  ifeq ($(PRECISION),DOUBLE)
    FFLAGS += -real_size 64
  endif
  FDEBF += -C
  FDEBF += -warn argument_checking
  FDEBF += -warn declarations
  ifneq ($(FC), f90)
    FDEBF += -warn truncated_source
    FDEBF += -warn unused
  endif
endif

include $(BOXLIB_HOME)/Src/Boundary/Make.package
INCLUDES += $(BOXLIB_HOME)/Src/Boundary

PATHDIRS  = $(HERE)
PATHDIRS += $(BOXLIB_HOME)/Src/Base
PATHDIRS += $(BOXLIB_HOME)/Src/Amr
PATHDIRS += $(BOXLIB_HOME)/Src/Extern/amrdata
PATHDIRS += $(BOXLIB_HOME)/Src/Boundary

ifeq (${BL_HAS_FORT}, TRUE)
FEXE_sources += ${EBASE}_F.F
endif

vpath %.f   $(PATHDIRS)
vpath %.F   $(PATHDIRS)
vpath %.f90 $(PATHDIRS)
vpath %.h   $(PATHDIRS)
vpath %.H   $(PATHDIRS)
vpath %.cpp $(PATHDIRS)

all: $(executable)

include $(BOXLIB_HOME)/Tools/GNUMake/Make.rules
