#
# Set these to the appropriate value.
#
DIM          = 2
DIM          = 3
PRECISION    = DOUBLE
COMP = PathScale
FCOMP = PathScale
COMP  = Intel
FCOMP = Intel
DEBUG        = TRUE
DEBUG        = FALSE
USE_MPI      = TRUE
USE_MPI      = FALSE
PROFILE      = TRUE
PROFILE      = FALSE

#DEFINES += -DBL_COALESCE_FABS

BOXLIB_HOME = ../..
include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

#ifeq ($(COMP),Intel)
#  COPTF   += -fast
#  CXXOPTF += -fast
#endif

#ifeq ($(FCOMP),Intel)
#  fOPTF += -fast
#  FOPTF += -fast
#endif

#FDEBF += -fpe0
#fDEBF += -fpe0

#FOPTF += -pg
#fOPTF += -pg
#CXXOPTF += -pg
#LDFLAGS += -pg
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
_progs  := tMF
#_progs  := tFB

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib
ifeq ($(COMP),g++)
  ifeq ($(__gcc_major_version),2)
    INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BaseLib/std
  endif
endif
VPATH = $(BOXLIB_HOME)/Src/C_BaseLib
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

FEXE_sources += util.F

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
	$(RM) *.ex *.o

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
