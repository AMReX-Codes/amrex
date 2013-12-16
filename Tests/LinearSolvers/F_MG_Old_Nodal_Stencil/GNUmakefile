BOXLIB_HOME := ../../..

NDEBUG := t
MPI    := t
OMP    :=
PROF   :=

USE_MG_CPP :=  

COMP := gfortran

MKVERBOSE := t

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

programs   += main.f90

f90sources += wrapper.f90
f90sources += cc_multi.f90
f90sources += cc_edge_coeffs.f90
f90sources += cc_rhs.f90
f90sources += nodal_multi.f90
f90sources += nodal_rhs.f90
f90sources += mt19937ar.f90
f90sources += init_cell_coeffs.f90

VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

include $(BOXLIB_HOME)/Src/LinearSolvers/F_MG/GPackage.mak
VPATH_LOCATIONS   += $(BOXLIB_HOME)/Src/LinearSolvers/F_MG
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/F_MG

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS   += $(BOXLIB_HOME)/Src/F_BaseLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

all:	$(pnames)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
