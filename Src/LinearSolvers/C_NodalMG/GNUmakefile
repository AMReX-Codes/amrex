#
# $Id: GNUmakefile,v 1.19 1998-11-10 02:21:59 car Exp $
#
PRECISION = DOUBLE
DEBUG     = FALSE
DEBUG     = TRUE
PROFILE   = FALSE
DIM       = 2
DIM       = 3
COMP	  = g++
COMP      = cxx
COMP      = KCC
USE_MPI   = FALSE
USE_MPI   = TRUE
#
# To build the executable: EBASE=proj LBASE=
#
EBASE = 

ifneq ($(EBASE),proj)
override EBASE = 
override LBASE = proj
endif

PBOXLIB_HOME = ..

TOP = $(PBOXLIB_HOME)
#
# Where libraries and include files will be installed.
#
INSTALL_ROOT = $(PBOXLIB_HOME)

include $(TOP)/mk/Make.defs ./Make.package

ifeq ($(LBASE),proj)
ifeq ($(KCC_VERSION),3.3)
CXXFLAGS += --one_instantiation_per_object
endif
endif

INCLUDE_LOCATIONS += $(TOP)/include

ifeq ($(LBASE),proj)
ifeq ($(COMP),KCC)
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libbox$(DIM)d.a
endif
all: $(optionsLib)
endif

ifeq ($(EBASE),proj)
CPPFLAGS          += -DHG_DEBUG
LIBRARIES         += -lproj$(optionsSuffix) -lbox$(DIM)d
LIBRARY_LOCATIONS += . $(TOP)/lib/$(machineSuffix)
CEXE_sources      += proj.cpp
all: $(executable)
endif

include $(TOP)/mk/Make.rules

ifeq ($(USE_MPI), TRUE)
DEFINES += -DBL_USE_MPI
LIBRARIES += -lmpi
ifeq ($(MACHINE), OSF1)
INCLUDE_LOCATIONS += /usr/local/mpi/include
LIBRARY_LOCATIONS += /usr/local/mpi/lib/alpha/ch_p4
endif
endif
#
# For 3rd ...
#
#3RD = 1
ifdef 3RD
CXXFLAGS += --link_command_prefix 3rd
CXXDEBF = +K0 --link_command_prefix 3rd
LIBRARIES += -ldnet_stub
FDEBF += -automatic
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
