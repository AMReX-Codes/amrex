#
# Set these to the appropriate value.
#
DIM          = 2
DIM          = 3

COMP         = g++
FCOMP        = gfortran

DEBUG        = TRUE
DEBUG        = FALSE

USE_MPI      = TRUE
USE_MPI      = FALSE

BOXLIB_HOME = ../..
include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

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
#_progs  := tMF
#_progs  := tFB
_progs  := tMFcopy

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib

ifeq ($(COMP),g++)
  ifeq ($(__gcc_major_version),2)
    INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib/std
  endif
endif

VPATH = $(BOXLIB_HOME)/Src/C_BaseLib
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
