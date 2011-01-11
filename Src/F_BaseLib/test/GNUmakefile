NDEBUG :=
MPI    :=
OMP    :=

COMP := Lahey
COMP := g95
COMP := Intel

FPARALLEL := ../..

include $(FPARALLEL)/mk/GMakedefs.mak

f90sources += t_main.f90
f90sources += t_knapsack.f90
f90sources += t_bx.f90
f90sources += t_cls.f90
f90sources += t_bxasc.f90
f90sources += t_bl_prof.f90
f90sources += t_particles.f90

programs += main.f90

all: $(pnames)

include $(FPARALLEL)/boxlib/GPackage.mak
include $(FPARALLEL)/extern/random/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib
VPATH_LOCATIONS += $(FPARALLEL)/extern/random
INCLUDE_LOCATIONS += $(FPARALLEL)/boxlib


include $(FPARALLEL)/mk/GMakerules.mak
