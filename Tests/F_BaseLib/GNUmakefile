BOXLIB_HOME := ../..

USE_RANDOM = TRUE

NDEBUG :=
MPI    := t
OMP    :=

COMP := Intel
COMP := gfortran

MKVERBOSE := t

PARTICLES := t

FPARALLEL := ../..

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

f90sources += t_main.f90
f90sources += t_knapsack.f90
f90sources += t_bx.f90
f90sources += t_cls.f90
f90sources += t_bxasc.f90
f90sources += t_bl_prof.f90
f90sources += t_particles.f90

programs   += main.f90
f90sources += mt19937ar.f90

all: $(pnames)

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS   += $(BOXLIB_HOME)/Src/F_BaseLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
