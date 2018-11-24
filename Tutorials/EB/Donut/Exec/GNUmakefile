AMREX_HOME := ../../../..

EBASE = EBL

DIM          = 3
COMP         = gnu

DEBUG        = FALSE

USE_MPI      = TRUE
USE_OMP      = FALSE

USE_EB       = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/EB/Make.package

include $(AMREX_HOME)/Tutorials/EB/LevelSet/Src/Make.package
VPATH_LOCATIONS += $(AMREX_HOME)/Tutorials/EB/LevelSet/Src
INCLUDE_LOCATIONS += $(AMREX_HOME)/Tutorials/EB/LevelSet/Src

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
