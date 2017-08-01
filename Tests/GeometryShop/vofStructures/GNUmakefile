#
# Set these to the appropriate value.
#
DIM          = 3
DIM          = 2

COMP         = g++

DEBUG        = FALSE
DEBUG        = TRUE

USE_MPI      = TRUE
USE_MPI      = FALSE

USE_OMP      = TRUE
USE_OMP      = FALSE

PROFILE       = FALSE
COMM_PROFILE  = FALSE
TRACE_PROFILE = FALSE
#DEFINES += -DBL_PROFILING_SPECIAL

AMREX_HOME = ../../..
include $(AMREX_HOME)/Tools/GNUMake/Make.defs

#EBASE  := umapTest
EBASE  := nbrsTest

DEFINES+=-DNCELLMAX=4
DEFINES+=-DNFACEMAX=5
ifeq ($(EBASE), umapTest)
  CEXE_headers += AMReX_EBStruct.H umapTest.H
  f90EXE_sources += umapTest_nd.f90 umapTest_mod.f90
  F90EXE_sources += AMReX_ebstruct_mod.F90
endif

ifeq ($(EBASE), nbrsTest)
  f90EXE_sources += nbrsTest_nd.f90
  F90EXE_sources += nbrsTest_mod.F90
endif
CEXE_sources += $(addsuffix .cpp, $(EBASE))

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/GeometryShop/Make.package
include $(AMREX_HOME)/Src/EBAMRTools/Make.package

all: $(addsuffix $(optionsSuffix).ex, $(EBASE))

clean::
	$(RM) bl_prof
	$(RM) *.ex

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

