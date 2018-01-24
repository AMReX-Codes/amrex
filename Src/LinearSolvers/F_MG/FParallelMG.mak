# FParallelMG.mak
#
# Use this to access the F90 solvers from C++.
# Note: that you can do a parallel make if you are using this file if you
# are using the MODDEP dependency evaluater that is in Parallel/scripts 
# and is used in the default Parallel/mk/Make.{rules,defs}. Otherwise,
# you can't since the dependencies for the object/module files can not be
# inferred from the names or ordering of thes files.
#
# ------------- in your GNUmakefile --------------
# 
# FAMREX_HOME =
# FAMREX_HOME = ../../../fParallel 
# ifdef FAMREX_HOME
#   include FParallelMG.mak
#   Fdirs   := boxlib mg extern/LAPACK
#   Flocs   := $(foreach dir, $(Fdirs), $(FAMREX_HOME)/$(dir))
# endif
#
# VPATH_LOCATIONS   += . $(Blocs) $(Flocs)
# INCLUDE_LOCATIONS += . $(Blocs) $(Flocs)
#
# ------------------------------------------------

DEFINES += -DAMREX_USE_FBOXLIB_MG

f90EXE_sources += compute_defect.f90
f90EXE_sources += coarsen_coeffs.f90
f90EXE_sources += mg_prolongation.f90
f90EXE_sources += ml_prolongation.f90
f90EXE_sources += cc_mg_cpp.f90
f90EXE_sources += cc_applyop.f90
f90EXE_sources += cc_ml_resid.f90
f90EXE_sources += cc_smoothers.f90
f90EXE_sources += cc_stencil.f90
f90EXE_sources += cc_stencil_apply.f90
f90EXE_sources += cc_stencil_fill.f90
f90EXE_sources += cc_interface_stencil.f90
f90EXE_sources += cc_mg_tower_smoother.f90
f90EXE_sources += itsol.f90
f90EXE_sources += mg.f90
f90EXE_sources += mg_tower.f90
f90EXE_sources += ml_cc.f90
f90EXE_sources += ml_nd.f90
f90EXE_sources += ml_norm.f90

f90EXE_sources += tridiag.f90

f90EXE_sources += nodal_mg_cpp.f90
f90EXE_sources += nodal_mask.f90
f90EXE_sources += nodal_divu.f90
f90EXE_sources += nodal_interface_stencil.f90
f90EXE_sources += nodal_newu.f90
f90EXE_sources += nodal_stencil.f90
f90EXE_sources += nodal_stencil_fill.f90
f90EXE_sources += nodal_smoothers.f90
f90EXE_sources += nodal_stencil_apply.f90
f90EXE_sources += nodal_sum.f90
f90EXE_sources += nodal_mg_tower_smoother.f90

f90EXE_sources += stencil_types.f90

#f90EXE_sources += mg_hypre_solve.f90
f90EXE_sources += nodal_mg_cpp.f90
f90EXE_sources += nodal_sync_resid.f90
cEXE_headers   += mg_cpp_f.h

f90EXE_sources += stencil_util.f90

VPATH_LOCATIONS += $(AMREX_HOME)/Src/LinearSolvers/F_MG
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/LinearSolvers/F_MG
