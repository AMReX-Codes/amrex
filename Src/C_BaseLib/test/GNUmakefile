#
# $Id: GNUmakefile,v 1.1 2001-07-20 19:37:36 car Exp $
#

#
# Set these to the appropriate value.
#
COMP         = KCC
COMP         = g++
DIM          = 2
PRECISION    = DOUBLE
DEBUG        = TRUE
USE_MPI      = FALSE

PBOXLIB_HOME = ../..
include $(PBOXLIB_HOME)/mk/Make.defs

MPI_HOME =

#
# Base name of each of the executables we want to build.
# I'm assuming that each of these is a stand-alone program,
# that simply needs to link against BoxLib.
#
#_progs = tVisMF tDir t8BIT tFB tFAC
_progs = tProfiler

all: $(addsuffix $(optionsSuffix).ex, $(_progs))

INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/BoxLib
VPATH = $(PBOXLIB_HOME)/BoxLib

$(addsuffix $(optionsSuffix).ex, $(_progs)) \
   : %$(optionsSuffix).ex : %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(libraries)
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

include $(PBOXLIB_HOME)/mk/Make.rules
