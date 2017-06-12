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

#
# Base name of each of the executables we want to build.
# I'm assuming that each of these is a stand-alone program,
# that simply needs to link against BoxLib.
#
_progs  := vofStructures

ifeq ($(_progs), vofStructures)
  DEFINES+=-DNCELLMAX=4
  DEFINES+=-DNFACEMAX=5
  F90EXE_sources += vofStructures_$(DIM)d.F90
  F90EXE_sources += AMReX_ebstruct_mod.F90
  CEXE_headers += Node.H
endif

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/GeometryShop/Make.package

all: $(addsuffix $(optionsSuffix).ex, $(_progs))


$(addsuffix $(optionsSuffix).ex, $(_progs)) \
   : %$(optionsSuffix).ex : %.cpp $(objForExecs)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(includes) $(LDFLAGS) -o $@ $< $(objForExecs) $(libraries)
	$(RM) $@.o

clean::
	$(RM) bl3_prof bl3_prof.m
	$(RM) *.ex *.o

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

