NDEBUG :=
MPI    :=
OMP    :=

MKVERBOSE:=t

COMP := g95
COMP := Lahey
COMP := Intel

FPARALLEL := ../..

include $(FPARALLEL)/mk/GMakedefs.mak
f90sources += t_main.f90
f90sources += t_ieee.f90

programs += main.f90

all: $(pnames)

include $(FPARALLEL)/boxlib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib
INCLUDE_LOCATIONS += $(FPARALLEL)/boxlib

include $(FPARALLEL)/mk/GMakerules.mak
