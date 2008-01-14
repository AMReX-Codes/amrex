# FParallelMG.mak
#
# Use this to access fParallel's fortran solvers in Parallel/BoxLib code.
# See iamrlib/run2d/GNUmakefile for an example.
# Note: that you can do a parallel make if you are using this file if you
# are using the MODDEP dependency evaluater that is in Parallel/scripts 
# and is used in the default Parallel/mk/Make.{rules,defs}. Otherwise,
# you can't since the dependencies for the object/module files can not be
# inferred from the names or ordering of thes files.
#
# ------------- in your GNUmakefile --------------
# 
# FBOXLIB_HOME =
# FBOXLIB_HOME = ../../../fParallel 
# ifdef FBOXLIB_HOME
#   include FParallelMG.mak
#   DEFINES += -DMG_USE_FBOXLIB
#   Fdirs   := boxlib mg extern/SPARSKIT extern/LAPACK
#   Flocs   := $(foreach dir, $(Fdirs), $(FBOXLIB_HOME)/$(dir))
# endif
#
# VPATH_LOCATIONS   += . $(Blocs) $(Flocs)
# INCLUDE_LOCATIONS += . $(Blocs) $(Flocs)
#
# ------------------------------------------------

f90EXE_sources += bc.f90
f90EXE_sources += bl_IO.f90
f90EXE_sources += bl_types.f90
f90EXE_sources += vector_i.f90
f90EXE_sources += sort_d.f90
ifeq ($(USE_MPI),TRUE)
  f90EXE_sources += parallel.f90
else
  f90EXE_sources += parallel_stubs.f90
endif
f90EXE_sources += bl_constants.f90
f90EXE_sources += bl_error.f90
f90EXE_sources += knapsack.f90
#f90EXE_sources += interp.f90
f90EXE_sources += bl_string.f90
f90EXE_sources += bl_timer.f90
f90EXE_sources += bl_stream.f90
f90EXE_sources += bl_mem_stat.f90
f90EXE_sources += box.f90
f90EXE_sources += sort_box.f90
f90EXE_sources += list_box.f90
f90EXE_sources += fab.f90
f90EXE_sources += sort_i.f90
f90EXE_sources += boxarray.f90
f90EXE_sources += ml_boxarray.f90
f90EXE_sources += layout.f90
f90EXE_sources += multifab.f90
f90EXE_sources += ml_layout.f90
f90EXE_sources += bndry_reg.f90
f90EXE_sources += ml_multifab.f90
f90EXE_sources += fabio.f90
f90EXE_sources += bl_space.f90
f90EXE_sources += bl_prof_stubs.f90

f90EXE_sources += mg_prolongation.f90
f90EXE_sources += st_coeffs.f90
f90EXE_sources += ml_prolongation.f90
f90EXE_sources += stencil.f90
f90EXE_sources += stencil_nodal.f90
f90EXE_sources += sparse_solve.f90
f90EXE_sources += ml_interface_stencil.f90
f90EXE_sources += ml_util.f90
f90EXE_sources += mg_smoother.f90
f90EXE_sources += mg_restriction.f90
f90EXE_sources += mg_defect.f90
f90EXE_sources += itsol.f90
f90EXE_sources += ml_restriction.f90
f90EXE_sources += nodal_mask.f90
f90EXE_sources += nodal_divu.f90
f90EXE_sources += nodal_newu.f90
f90EXE_sources += mg.f90
f90EXE_sources += mlmg.f90
f90EXE_sources += ml_nd.f90
f90EXE_sources += ml_cc.f90
f90EXE_sources += mg_cpp.f90
f90EXE_sources += mg_nodal_cpp.f90
cEXE_headers   += mg_cpp_f.h

cEXE_sources   += fabio_c.c
cEXE_sources   += timer_c.c

fEXE_sources += daxpy.f
fEXE_sources += ddot.f
fEXE_sources += dnrm2.f
fEXE_sources += ilut.f
fEXE_sources += iters.f
fEXE_sources += sk_sup.f
