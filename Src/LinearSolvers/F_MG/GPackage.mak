f90sources += itsol.f90
f90sources += mg_cpp.f90
f90sources += mg.f90
f90sources += mg_smoother.f90
f90sources += mlmg.f90
f90sources += sparse_solve.f90
f90sources += stencil.f90
f90sources += st_coeffs.f90
f90sources += stencil_nodal.f90
f90sources += mg_defect.f90
f90sources += mg_prolongation.f90
f90sources += mg_restriction.f90
f90sources += ml_interface_stencil.f90
f90sources += ml_prolongation.f90
f90sources += ml_restriction.f90
f90sources += ml_util.f90

include $(FPARALLEL)/extern/SPARSKIT/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/SPARSKIT
