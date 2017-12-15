#
# Set these to the appropriate value.
#
DIM          = 3
DIM          = 2
USE_EB = FALSE

COMP         = g++

DEBUG        = FALSE
DEBUG        = TRUE

USE_MPI      = TRUE
USE_MPI      = FALSE

PROFILE       = FALSE
COMM_PROFILE  = FALSE
TRACE_PROFILE = FALSE
#DEFINES += -DBL_PROFILING_SPECIAL

AMREX_HOME = ../../../../..
include $(AMREX_HOME)/Tools/GNUMake/Make.defs

#
# Base name of each of the executables we want to build.
# I'm assuming that each of these is a stand-alone program,
# that simply needs to link against BoxLib.
#
#_progs  := tFB
#_progs  := tRABcast.cpp
_progs  += fourthOrderCFInterpTest



include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
#include $(AMREX_HOME)/Src/Amr/Make.package


INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Boundary
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/AmrCore
#INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Amr
VPATH_LOCATIONS += $(AMREX_HOME)/Src/Base
VPATH_LOCATIONS += $(AMREX_HOME)/Src/Boundary
VPATH_LOCATIONS += $(AMREX_HOME)/Src/AmrCore
#VPATH_LOCATIONS += $(AMREX_HOME)/Src/Amr

all: $(addsuffix $(optionsSuffix).ex, $(_progs))


$(addsuffix $(optionsSuffix).ex, $(_progs)) \
   : %$(optionsSuffix).ex : %.cpp $(objForExecs)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(includes) $(LDFLAGS) -o $@ $< $(objForExecs) $(libraries)
	$(RM) $@.o

clean::
	$(RM) bl3_prof bl3_prof.m
	$(RM) *.ex *.o

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

