PRECISION = DOUBLE
DEBUG     = FALSE
DEBUG     = TRUE
PROFILE   = FALSE
DIM       = 2
DIM       = 3
COMP      = cxx
COMP      = KCC


EBASE=proj
LBASE=

GUTHAMR_HOME = ..
PBOXLIB_HOME = ..
TOP = $(GUTHAMR_HOME)

INCLUDE_LOCATIONS += ../pBoxLib_2

include $(TOP)/mk/Make.defs

CPPFLAGS += -DBL_USE_NEW_HFILES

include ./Make.package
include ../pBoxLib_2/Make.package

vpath %.cpp :	. ../pBoxLib_2
vpath %.F   :	. ../pBoxLib_2

ifeq ($(LBASE),proj)
  optionsSuffix      = $(DIM)d.$(machineSuffix)
  ifeq ($(COMP),KCC)
    LibsToCloseAgainst += ../pBoxLib_2/libbox$(DIM)d.$(machineSuffix).a
  endif
  all: $(instLib)
endif

ifeq ($(EBASE),proj)
  all: $(executable)
endif

ifeq ($(MACHINE),Linux)
  LIBRARIES += $(LINUXLIBS)
endif

#3RD = 1
ifdef 3RD
# FOR RUNNING 3RD ONLY
CXXFLAGS += --link_command_prefix 3rd
CXXDEBF = +K0 --link_command_prefix 3rd
LIBRARIES += -ldnet_stub
FDEBF += -automatic
# FOR RUNNING 3RD ONLY
endif
FDEBF += -u

include $(TOP)/mk/Make.rules
