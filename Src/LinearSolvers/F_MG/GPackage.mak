f90sources += stencil_types.f90
f90sources += compute_defect.f90
f90sources += coarsen_coeffs.f90
f90sources += itsol.f90

ifdef USE_MG_CPP
f90sources += cc_mg_cpp.f90
f90sources += nodal_mg_cpp.f90
f90sources += nodal_sync_resid.f90
endif

f90sources += cc_applyop.f90
f90sources += cc_smoothers.f90
f90sources += cc_interface_stencil.f90
f90sources += cc_ml_resid.f90
f90sources += cc_stencil.f90
f90sources += cc_stencil_apply.f90
f90sources += cc_stencil_fill.f90
f90sources += cc_mg_tower_smoother.f90

f90sources += mg.f90
f90sources += mg_tower.f90
f90sources += mg_prolongation.f90

f90sources += ml_prolongation.f90
f90sources += ml_nd.f90
f90sources += ml_cc.f90
f90sources += ml_norm.f90
f90sources += ml_solve.f90

f90sources += nodal_applyop.f90
f90sources += nodal_mask.f90
f90sources += nodal_mg_cpp.f90
f90sources += nodal_divu.f90
f90sources += nodal_enforce_dirichlet_rhs.f90
f90sources += nodal_interface_stencil.f90
f90sources += nodal_newu.f90
f90sources += nodal_stencil.f90
f90sources += nodal_stencil_fill.f90
f90sources += nodal_smoothers.f90
f90sources += nodal_stencil.f90
f90sources += nodal_stencil_apply.f90
f90sources += nodal_sum.f90
f90sources += nodal_mg_tower_smoother.f90

f90sources += tridiag.f90

f90sources += stencil_util.f90
 
#include $(FPARALLEL)/extern/SPARSKIT/GPackage.mak
#VPATH_LOCATIONS += $(FPARALLEL)/extern/SPARSKIT
 
#include $(FPARALLEL)/extern/BLAS/GPackage.mak
#VPATH_LOCATIONS += $(FPARALLEL)/extern/BLAS
