AMREX_HOME := ../../..

PRECISION = DOUBLE

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 3
DIM	= 2

COMP =g++

USE_MPI=FALSE

EBASE = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += $(EBASE).cpp

include ./Make.package

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/GeometryShop/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
