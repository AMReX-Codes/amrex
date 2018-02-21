#
# Generic setup for using gcc
#
ifeq ($(USE_OMP),TRUE)
  CXX = xlC_r
  CC  = xlc_r
  FC  = xlf_r
  F90 = xlf_r
else
  CXX = xlC
  CC  = xlc
  FC  = xlf
  F90 = xlf
endif

CXXFLAGS =
CFLAGS   =
FFLAGS   =
F90FLAGS =

########################################################################

ibm_version  = $(shell $(CXX) --version | head -1)

COMP_VERSION = $(ibm_version)


########################################################################

ifeq ($(DEBUG),TRUE)

  CXXFLAGS += -g -O0 
  CFLAGS   += -g -O0 

  FFLAGS   += -g -O0 
  F90FLAGS += -g -O0 

else

  CXXFLAGS += -g -O2
  CFLAGS   += -g -O2
  FFLAGS   += -g -O2
  F90FLAGS += -g -O2

endif


########################################################################

CXXFLAGS += -std=c++1y
CFLAGS   += -std=gnu99
F90FLAGS += -qlanglvl=extended -qxlf2003=polymorphic

FFLAGS   += -qmoddir=$(fmoddir) -I $(fmoddir) -WF,-C!
F90FLAGS += -qmoddir=$(fmoddir) -I $(fmoddir) -WF,-C!

########################################################################

GENERIC_COMP_FLAGS =


#ifeq ($(USE_OMP),TRUE)
#  GENERIC_COMP_FLAGS += -fopenmp
#endif

CXXFLAGS += $(GENERIC_COMP_FLAGS)
CFLAGS   += $(GENERIC_COMP_FLAGS)
FFLAGS   += $(GENERIC_COMP_FLAGS)
F90FLAGS += $(GENERIC_COMP_FLAGS)

CPP_PREFIX = -WF,

#override XTRALIBS += 
override XTRALIBS = $(shell mpifort -showme:link) -L $(OLCF_XLF_ROOT)/lib -lxlf90_r -lm  -lxlfmath

FORTLINK = LOWERCASE
