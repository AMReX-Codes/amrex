#
# Set these to the appropriate value.
#
DIM          = 3
DIM          = 2

COMP         = g++
FCOMP        = gfortran

DEBUG        = TRUE
DEBUG        = FALSE

USE_MPI      = TRUE
USE_MPI      = FALSE

PROFILE       = TRUE
COMM_PROFILE  = FALSE
TRACE_PROFILE = FALSE


BOXLIB_HOME = ../..
include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

CXXFLAGS += -std=c++0x

#
# Base name of each of the executables we want to build.
# I'm assuming that each of these is a stand-alone program,
# that simply needs to link against BoxLib.
#
#_progs += tVisMF tDir t8BIT tFB tFAC tCArena
#_progs += tRan
#_progs  := tread
#_progs  := tParmParse
#_progs  := tCArena
#_progs  := tBA
#_progs  := tDM
#_progs  := tFillFab
#_progs  := tMF
#_progs  := tFB
#_progs  := tMFcopy
_progs  := tProfiler

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BoundaryLib

ifeq ($(COMP),g++)
  ifeq ($(__gcc_major_version),2)
    INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib/std
  endif
endif

ifeq ($(_progs),tProfiler)
  FEXE_sources += TPROFILER.F
endif

ifeq ($(_progs),tFillFab)
  fEXE_sources += fillfab.f
endif

VPATH += $(BOXLIB_HOME)/Src/C_BaseLib

include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

all: $(addsuffix $(optionsSuffix).ex, $(_progs))


$(addsuffix $(optionsSuffix).ex, $(_progs)) \
   : %$(optionsSuffix).ex : %.cpp $(objForExecs)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(objForExecs) $(libraries)
	$(RM) $@.o

clean::
	$(RM) bl3_prof bl3_prof.m
	$(RM) *.ex *.o

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
