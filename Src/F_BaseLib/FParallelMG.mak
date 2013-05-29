# FParallelMG.mak
#
# Use this to access include the F_BaseLib files that you need
# when calling F90 code from C++.
#
# Note: you can do a parallel make if you are using this file if you
# are using the MODDEP dependency evaluater that is in Tools/C_scripts 
# and is used in the default Tools/C_mk/Make.{rules,defs}. Otherwise,
# you can't since the dependencies for the object/module files can not be
# inferred from the names or ordering of thes files.
#
# ------------------------------------------------

f90EXE_sources += bc.f90
f90EXE_sources += bl_constants.f90
f90EXE_sources += bl_error.f90
f90EXE_sources += bl_IO.f90
f90EXE_sources += bl_stream.f90
f90EXE_sources += bl_string.f90
f90EXE_sources += bl_timer.f90
f90EXE_sources += bl_types.f90
f90EXE_sources += bl_mem_stat.f90
f90EXE_sources += bl_space.f90
f90EXE_sources += bl_prof_stubs.f90
f90EXE_sources += bndry_reg.f90
f90EXE_sources += box_f.f90
f90EXE_sources += boxarray_f.f90
f90EXE_sources += fab.f90
f90EXE_sources += fabio.f90
f90EXE_sources += knapsack.f90
f90EXE_sources += layout.f90
f90EXE_sources += list_box.f90
f90EXE_sources += ml_boxarray.f90
f90EXE_sources += ml_layout.f90
f90EXE_sources += ml_multifab.f90
f90EXE_sources += multifab_f.f90
f90EXE_sources += sort_box.f90
f90EXE_sources += sort_d.f90
f90EXE_sources += sort_i.f90
f90EXE_sources += vector_i.f90

ifeq ($(USE_OMP),TRUE)
  f90EXE_sources += omp.f90
else
  f90EXE_sources += omp_stubs.f90
endif

ifeq ($(USE_MPI),TRUE)
  f90EXE_sources += parallel.f90
else
  f90EXE_sources += parallel_stubs.f90
endif

cEXE_sources   += fabio_c.c
cEXE_sources   += timer_c.c
