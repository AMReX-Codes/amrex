#
# $Id: GNUmakefile,v 1.39 1999-08-06 23:10:05 car Exp $
#
PBOXLIB_HOME = ..

TOP = $(PBOXLIB_HOME)

PRECISION = DOUBLE
DEBUG     = FALSE
DEBUG     = TRUE
PROFILE   = FALSE
DIM       = 2
DIM       = 3
COMP      = cxx
COMP	  = g++
COMP      = KCC
USE_MPI   = FALSE
USE_MPI   = TRUE

#
# To build the executable: EBASE=proj LBASE=
#
EBASE = 
EBASE = proj

ifneq ($(EBASE),proj)
  override EBASE = 
  override LBASE = proj
endif
#
# Where libraries and include files will be installed.
#
INSTALL_ROOT = $(PBOXLIB_HOME)

include $(TOP)/mk/Make.defs
include ./Make.package

ifeq ($(LBASE),proj)
  ifeq ($(KCC_VERSION),3.3)
    CXXFLAGS += --one_instantiation_per_object
  endif
  INCLUDE_LOCATIONS += $(TOP)/include

  ifeq ($(COMP),KCC)
    LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libbox$(DIM)d.a
  endif
  all: $(optionsLib)
endif

ifeq ($(EBASE),proj)
  include $(TOP)/bndrylib/Make.package
  include $(TOP)/pBoxLib_2/Make.package

  INCLUDE_LOCATIONS += $(TOP)/bndrylib
  INCLUDE_LOCATIONS += $(TOP)/pBoxLib_2 

  ifeq ($(DEBUG),TRUE)
    CPPFLAGS          += -DHG_DEBUG
  endif
  ifeq ($(COMP),KCC)
    CXXFLAGS += --strict
    CXXFLAGS += --display_error_number
    CXXFLAGS += --diag_suppress 450
  endif

  LIBRARIES         +=
  LIBRARY_LOCATIONS += 
  CEXE_sources      += proj.cpp

  vpath %.cpp	. $(TOP)/pBoxLib_2 $(TOP)/bndrylib
  vpath %.H     . $(TOP)/pBoxLib_2 $(TOP)/bndrylib
  vpath %.F	.                  $(TOP)/bndrylib
  all: $(executable)
endif

include $(TOP)/mk/Make.rules

#LDFLAGS += -v

#
# For 3rd ...
#
#3RD = 0
ifdef 3RD
  CXXFLAGS  += --link_command_prefix 3rd
  LIBRARIES += -ldnet_stub
  FDEBF     += -automatic
endif

ifeq ($(MACHINE),T3E)
  #
  # Turn off warning about double precision constants.
  #
  FOPTF += -M 1110
  FDEBF += -M 1110
endif

ifeq ($(MACHINE),OSF1)
  #
  # Warn about unused variables.
  #
  FDEBF += -u
endif

#
# Lower optimization level for this big file
#
o/3d.OSF1.KCC/hg_multi3d_terrain.o: hg_multi3d_terrain.F
	@if [ ! -d $(objTempDir) ]; then mkdir -p $(objTempDir); fi
	@if [ ! -d $(f77TempDir) ]; then mkdir -p $(f77TempDir); fi
	$(CPP) -DBL_LANG_FORT $(CPPFLAGS) $< | $(FORTPREP) > f/3d.OSF1.KCC/hg_multi3d_terrain.f
	$(COMPILE.f) $(FEXEFLAGS) -O f/3d.OSF1.KCC/hg_multi3d_terrain.f $(FORT_OUTPUT_OPTION)

libs:
	cd $(TOP)/pBoxLib_2; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install

