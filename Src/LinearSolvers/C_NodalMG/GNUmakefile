#
# $Id: GNUmakefile,v 1.15 1998-10-28 22:47:55 lijewski Exp $
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
# To build the closed library: EBASE= LBASE=proj
#
EBASE=proj
LBASE=

GUTHAMR_HOME = ..
PBOXLIB_HOME = ..

TOP = $(GUTHAMR_HOME)

include $(TOP)/mk/Make.defs ./Make.package

CPPFLAGS += -DHG_DEBUG

ifeq ($(KCC_VERSION),3.3)
CXXFLAGS += --one_instantiation_per_object
endif

INCLUDE_LOCATIONS += $(TOP)/include

ifeq ($(LBASE),proj)
ifeq ($(COMP),KCC)
LibsToCloseAgainst += ../lib/$(machineSuffix)/libbox$(DIM)d.a
endif
all: $(instLib)
endif

ifeq ($(EBASE),proj)
  all: $(executable)
endif

ifeq ($(USE_MPI), TRUE)
DEFINES += -DBL_USE_MPI
LIBRARIES += -lmpi
ifeq ($(MACHINE), OSF1)
MPI_HOME = /usr/local/mpi
INCLUDE_LOCATIONS += $(MPI_HOME)/include
LIBRARY_LOCATIONS += $(MPI_HOME)/lib/alpha/ch_p4
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

include $(TOP)/mk/Make.rules
