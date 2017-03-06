# AMREX_HOME defines the directory in which we will find all the BoxLib code
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../..

NDEBUG    := t
MPI       := t
OMP       :=
PROF      :=
COMP      := gfortran
MKVERBOSE := t

RANDOM := t  # this makes it include mt19937ar.f90
include $(AMREX_HOME)/Tools/F_mk/GMakedefs.mak

f90sources += main.f90
f90sources += constants.f90
f90sources += init_coeffs.f90
f90sources += init_rhs.f90
f90sources += cycle_tower.f90
f90sources += traverse.f90
f90sources += v_cycle.f90
f90sources += fmg_cycle.f90
f90sources += solve.f90

VPATH_LOCATIONS += .

include $(AMREX_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(AMREX_HOME)/Src/F_BaseLib

main.$(suf).exe: $(objects) 
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(AMREX_HOME)/Tools/F_mk/GMakerules.mak
