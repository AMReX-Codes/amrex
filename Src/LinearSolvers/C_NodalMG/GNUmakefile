PRECISION = DOUBLE
DEBUG     = FALSE
DEBUG     = TRUE
PROFILE   = FALSE
DIM       = 2
DIM       = 3
COMP	  = g++
COMP      = cxx
COMP      = KCC

USE_MPI=FALSE
USE_MPI=TRUE

EBASE=proj
LBASE=

GUTHAMR_HOME = ..
PBOXLIB_HOME = ..
TOP = $(GUTHAMR_HOME)

INCLUDE_LOCATIONS += ../pBoxLib_2

include $(TOP)/mk/Make.defs

CPPFLAGS += -DBL_USE_NEW_HFILES
CPPFLAGS += -DHG_DEBUG

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


ifeq ($(USE_MPI), TRUE)
DEFINES += -DBL_USE_MPI
MPI_HOME = /usr/local/mpi
INCLUDE_LOCATIONS += $(MPI_HOME)/include
LIBRARY_LOCATIONS += $(MPI_HOME)/lib/alpha/ch_p4
LIBRARIES += -lmpi
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
FDEBF += -u -C

# Lower optimization level for this big file
o/3d.OSF1.KCC/hg_multi3d_terrain.o: hg_multi3d_terrain.F
	@if [ ! -d $(objTempDir) ]; then mkdir -p $(objTempDir); fi
	@if [ ! -d $(f77TempDir) ]; then mkdir -p $(f77TempDir); fi
	$(CPP) -DBL_LANG_FORT $(CPPFLAGS) $< | $(FORTPREP) > f/3d.OSF1.KCC/hg_multi3d_terrain.f
	$(COMPILE.f) $(FEXEFLAGS) -O f/3d.OSF1.KCC/hg_multi3d_terrain.f $(FORT_OUTPUT_OPTION)

include $(TOP)/mk/Make.rules
