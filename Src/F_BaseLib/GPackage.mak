f90sources += boxlib_f.f90
f90sources += f2kcli$(f2kcli_suf).f90

f90sources += bl_constants.f90
f90sources += bl_IO.f90
ifdef BACKTRACE
  f90sources += bl_prof_backtrace.f90
  f90sources += bl_error_backtrace.f90
else
ifdef PROF
  f90sources += bl_prof.f90
else
  f90sources += bl_prof_stubs.f90
endif
f90sources += bl_error.f90
endif
f90sources += bl_mem_stat.f90
f90sources += bl_parmparse.f90
f90sources += bl_space.f90
f90sources += bl_stream.f90
f90sources += bl_string.f90
f90sources += bl_system.f90
f90sources += bl_timer.f90
f90sources += bl_types.f90

f90sources += box_f.f90
f90sources += knapsack.f90
f90sources += layout.f90
f90sources += boxarray_f.f90
f90sources += box_util.f90
f90sources += fab.f90
f90sources += multifab_f.f90
f90sources += fabio.f90
f90sources += plotfile.f90
f90sources += filler.f90
f90sources += cluster_f.f90
f90sources += interp.f90
f90sources += fourth_order_interp_coeffs.f90
f90sources += bc.f90
f90sources += bndry_reg.f90

f90sources += ml_boxarray.f90
f90sources += ml_layout.f90
f90sources += ml_multifab.f90

f90sources += create_umac_grown.f90
f90sources += define_bc_tower.f90
f90sources += fillpatch.f90
f90sources += multifab_fill_ghost_cells.f90
f90sources += multifab_physbc_edgevel.f90
f90sources += multifab_physbc.f90

f90sources += list_box.f90
f90sources += sort_box.f90
f90sources += vector_i.f90
f90sources += sort_d.f90
f90sources += sort_i.f90

f90sources += ppm_util.f90

f90sources += cutcells.f90

f90sources += make_new_grids.f90
f90sources += tag_boxes.f90

ifdef PARTICLES
  f90sources += particles_f.f90
endif

ifndef MPI
  f90sources += parallel_stubs.f90
else
  f90sources += parallel.f90
  ifeq ($(ARCH),Darwin)
#    include $(FPARALLEL)/extern/MacMPI/GPackage.mak
  endif
endif

ifdef OMP
  f90sources += omp.f90
else
  f90sources += omp_stubs.f90
endif

csources += fabio_c.c
csources += timer_c.c
csources += ppm_util_c.c
csources += system_util_c.c

ifeq ($(ARCH),AIX)
endif

ifeq ($(ARCH),OSF1)
endif
