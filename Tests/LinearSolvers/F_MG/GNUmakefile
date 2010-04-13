NDEBUG :=
MPI    :=
OMP    :=
MPI    :=  
PROF   :=

NDEBUG :=

USE_MG_CPP :=  

COMP := g95
COMP := PathScale
COMP := Intel

MKVERBOSE := t

FPARALLEL := ../..

include $(FPARALLEL)/mk/GMakedefs.mak

programs   += main.f90

f90sources += wrapper.f90
f90sources += cc_single.f90
f90sources += cc_multi.f90
f90sources += nodal_single.f90
f90sources += nodal_multi.f90
f90sources += t_smoother.f90
#f90sources += bc_interp.f90

include $(FPARALLEL)/mg/test//GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/mg/test

include $(FPARALLEL)/boxlib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib

include $(FPARALLEL)/mg/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/mg

include $(FPARALLEL)/extern/random/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/random

all:	$(pnames)

include $(FPARALLEL)/mk/GMakerules.mak
