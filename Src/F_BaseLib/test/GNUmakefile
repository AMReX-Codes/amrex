NDEBUG :=
MPI    :=
OMP    :=

COMP := g95
COMP := Lahey
COMP := Intel

FPARALLEL := ../..

include $(FPARALLEL)/mk/GMakedefs.mak
f90sources += t_main.f90
f90sources += t_knapsack.f90

programs += main.f90

all: $(pnames)

include $(FPARALLEL)/boxlib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib
INCLUDE_LOCATIONS += $(FPARALLEL)/boxlib

ifdef USE_METIS
include $(FPARALLEL)/extern/metis-4.0/Lib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/metis-4.0/Lib
INCLUDE_LOCATIONS += $(FPARALLEL)/extern/metis-4.0/Lib
endif

include $(FPARALLEL)/mk/GMakerules.mak
