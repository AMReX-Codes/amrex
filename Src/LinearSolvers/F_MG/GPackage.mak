f90sources += itsol.f90
f90sources += mg_bc.f90
f90sources += mg_cpp.f90
f90sources += mg.f90
f90sources += mg_smoother.f90
f90sources += ml.f90
f90sources += mlmg.f90
f90sources += sparse_solve.f90
f90sources += stencil.f90
f90sources += st_coeffs.f90
f90sources += stencil_nodal.f90

include $(FPARALLEL)/extern/SPARSKIT/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/SPARSKIT
