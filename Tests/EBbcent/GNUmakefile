AMREX_HOME ?= ../../

#DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 2

COMP    = gnu

USE_MPI   = FALSE
USE_OMP   = FALSE
USE_EB    = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/GeometryShop/Make.package
include $(AMREX_HOME)/Src/EB/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
