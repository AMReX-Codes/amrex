#
# Set these to the appropriate value.
#
DIM          = 3
DIM          = 2

COMP         = g++

DEBUG        = TRUE
DEBUG        = FALSE

USE_MPI      = TRUE
USE_MPI      = FALSE

USE_OMP      = FALSE
USE_OMP      = TRUE

PROFILE       = FALSE
COMM_PROFILE  = FALSE
TRACE_PROFILE = FALSE
TINY_PROFILE  = TRUE
#DEFINES += -DBL_PROFILING_SPECIAL

AMREX_HOME = ../../../..
include $(AMREX_HOME)/Tools/GNUMake/Make.defs

#_progs  := stencilTest
#_progs  := stencilTestMSD
_progs  := dirichletTest

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/EB/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/GeometryShop/Make.package
include ../src/Make.package
VPATH_LOCATIONS   += ../src
INCLUDE_LOCATIONS += ../src
all: $(addsuffix $(optionsSuffix).ex, $(_progs))


$(addsuffix $(optionsSuffix).ex, $(_progs)) \
   : %$(optionsSuffix).ex : %.cpp $(objForExecs)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(includes) $(LDFLAGS) -o $@ $< $(objForExecs) $(libraries)
	$(RM) $@.o

clean::
	$(RM) bl3_prof bl3_prof.m
	$(RM) *.ex *.o

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

