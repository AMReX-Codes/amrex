AMREX_HOME ?= ../../..

DEBUG ?= FALSE
DIM ?= 3
COMP ?= gnu
USE_MPI ?= TRUE
USE_OMP ?= FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include Make.package
Pdirs := Base
Ppack += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)
include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

