#
# Set these to the appropriate value.
#
DIM          = 3

COMP         = g++

DEBUG        = TRUE

USE_MPI      = FALSE
USE_MPI      = TRUE

USE_OMP      = FALSE

AMREX_HOME = ../../..
include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/GeometryShop/Make.package
include $(AMREX_HOME)/Src/EB/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

