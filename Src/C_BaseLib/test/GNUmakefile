#
# $Id: GNUmakefile,v 1.4 2001-07-22 18:11:03 car Exp $
#

#
# Set these to the appropriate value.
#
COMP         = KCC
COMP         = g++
DIM          = 2
PRECISION    = DOUBLE
DEBUG        = TRUE
USE_MPI      = TRUE
USE_MPI      = FALSE
USE_THREADS = TRUE

PBOXLIB_HOME = ../..
include $(PBOXLIB_HOME)/mk/Make.defs

DEFINES += -DBL_PTHREADS
DEFINES += -DBL_PROFILING

#
# Base name of each of the executables we want to build.
# I'm assuming that each of these is a stand-alone program,
# that simply needs to link against BoxLib.
#
#_progs = tVisMF tDir t8BIT tFB tFAC
#_progs = tProfiler tThread tWorkQueue
_progs = tRan

INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/BoxLib
ifeq ($(COMP),g++)
  INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/BoxLib/std
endif
VPATH = $(PBOXLIB_HOME)/BoxLib
include $(PBOXLIB_HOME)/BoxLib/Make.package

all: $(addsuffix $(optionsSuffix).ex, $(_progs))

$(addsuffix $(optionsSuffix).ex, $(_progs)) \
   : %$(optionsSuffix).ex : %.cpp $(objForExecs)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(objForExecs) $(libraries)
	$(RM) $@.o

#
# Run the tests ...
#
runtests : $(addsuffix $(optionsSuffix).ex, $(_progs))
        #
        # 
        #
	for f in $(filter-out tFABiofilt%,                 \
          $(addsuffix $(optionsSuffix).ex, $(_progs))); do \
          echo "====> Testing $$f ..."; echo; ./$$f; echo; \
        done

clean::
	$(RM) bl3_prof bl3_prof.m

include $(PBOXLIB_HOME)/mk/Make.rules
