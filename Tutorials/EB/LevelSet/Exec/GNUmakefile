AMREX_HOME := ../../../..
TOP := $(AMREX_HOME)/Tutorials/EB/LevelSet

EBASE = EBL

DIM          = 3
COMP         = gnu

DEBUG        = FALSE

USE_MPI      = TRUE

PROFILE       = FALSE
COMM_PROFILE  = FALSE
TRACE_PROFILE = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/EB/Make.package
include $(AMREX_HOME)/Src/GeometryShop/Make.package
include $(AMREX_HOME)/Src/Extern/amrdata/Make.package
include $(AMREX_HOME)/Tools/C_util/Make.package
include $(AMREX_HOME)/Tutorials/EB/LevelSet/Src/Make.package

INCLUDE_LOCATIONS += .
vpathdir += .

INCLUDE_LOCATIONS += ../Src
vpathdir += ../Src

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
vpathdir += $(AMREX_HOME)/Src/Base

INCLUDE_LOCATIONS += $(AMREX_HOME)/Tools/C_util
vpathdir += $(AMREX_HOME)/Tools/C_util

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Boundary
vpathdir += $(AMREX_HOME)/Src/Boundary

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/EB
vpathdir += $(AMREX_HOME)/Src/EB

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/GeometryShop
vpathdir += $(AMREX_HOME)/Src/GeometryShop

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/AmrCore
vpathdir += $(AMREX_HOME)/Src/AmrCore

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata
vpathdir += $(AMREX_HOME)/Src/Extern/amrdata

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

all: $(executable)
	$(SILENT) $(RM) AMReX_buildInfo.cpp
	@echo SUCCESS

clean::
	$(RM) bl3_prof bl3_prof.m
	$(RM) *.ex *.o
