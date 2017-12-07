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

PROFILE       = FALSE
COMM_PROFILE  = FALSE
TRACE_PROFILE = FALSE
#DEFINES += -DBL_PROFILING_SPECIAL

AMREX_HOME = ../../../..
include $(AMREX_HOME)/Tools/GNUMake/Make.defs

#
# Base name of each of the executables we want to build.
# I'm assuming that each of these is a stand-alone program,
# that simply needs to link against BoxLib.
#
_progs  += sphere
_progs  += interiorBox
_progs  += surfaceOfRevolution
_progs  += parabolaWithSphere

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/EB/Make.package
include $(AMREX_HOME)/Src/GeometryShop/Make.package
include $(AMREX_HOME)/Src/Extern/amrdata/Make.package
include $(AMREX_HOME)/Src/Extern/ProfParser/Make.package
include $(AMREX_HOME)/Tools/C_util/Make.package
include $(AMREX_HOME)/Tutorials/EB/GeometryGeneration/src/Make.package

INCLUDE_LOCATIONS += .
vpathdir += .

INCLUDE_LOCATIONS += $(AMREX_HOME)/Tutorials/EB/GeometryGeneration/src
vpathdir += $(AMREX_HOME)/Tutorials/EB/GeometryGeneration/src


INCLUDE_LOCATIONS += ../src
vpathdir += ../src

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
vpathdir += $(AMREX_HOME)/Src/Base


INCLUDE_LOCATIONS += $(AMREX_HOME)/Tools/C_util
vpathdir += $(AMREX_HOME)/Tools/C_util

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Boundary
vpathdir += $(AMREX_HOME)/Src/Boundary

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/EB
vpathdir += $(AMREX_HOME)/Src/EB

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/GeometryShop
vpathdir += $(AMREX_HOME)/Src/GeometryShop

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/AmrCore
vpathdir += $(AMREX_HOME)/Src/AmrCore


INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata
vpathdir += $(AMREX_HOME)/Src/Extern/amrdata

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/ProfParser
vpathdir += $(AMREX_HOME)/Src/Extern/Extern/ProfParser


vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(addsuffix $(optionsSuffix).ex, $(_progs))


$(addsuffix $(optionsSuffix).ex, $(_progs)) \
   : %$(optionsSuffix).ex : %.cpp $(objForExecs)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(includes) $(LDFLAGS) -o $@ $< $(objForExecs) $(libraries)
	$(RM) $@.o

clean::
	$(RM) bl3_prof bl3_prof.m
	$(RM) *.ex *.o

include $(AMREX_HOME)/Tools/GNUMake/Make.rules






