set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")


set (CXXSRC backtrace_c.cpp)

set (CSRC  fabio_c.c  ppm_util_c.c  system_util_c.c  timer_c.c)

set (F90SRC
   f2kcli.f90 bc.f90 backtrace_f.f90 bc_functions.f90 bl_constants.f90
   bl_error.f90 bl_IO.f90 bl_mem_stat.f90 BLProfiler_f90.f90 bl_space.f90
   bl_stream.f90 bl_string.f90 bl_system.f90 bl_timer.f90 bl_types.f90
   bndry_reg.f90 boxarray_f.f90 box_f.f90 boxlib_f.f90 box_util.f90 cluster_f.f90
   create_umac_grown.f90 cutcells.f90 define_bc_tower.f90 fab.f90 amrex_fabio_c.f90
   amrex_timer_c.f90 fabio.f90 filler.f90 fillpatch.f90 fourth_order_interp_coeffs.f90
   interp.f90 knapsack.f90 layout.f90 list_box.f90 make_new_grids.f90 ml_boxarray.f90
   ml_layout.f90 ml_multifab.f90 ml_restrict_fill.f90 cc_restriction.f90 nodal_restriction.f90
   ml_cc_restriction.f90 ml_nd_restriction.f90 multifab_f.f90 multifab_fill_ghost_cells.f90
   multifab_physbc_edgevel.f90 multifab_physbc.f90 nodal_neumann_bcs.f90 nodal_stencil_bc.f90
   plotfile.f90 ppm_util.f90 regrid.f90 sort_box.f90 sort_d.f90 sort_i.f90 tag_boxes.f90
   vector_i.f90)

if (ENABLE_MPI)
  list (APPEND F90SRC parallel.f90)
else()
  list (APPEND F90SRC parallel_stubs.f90)
endif()

if (ENABLE_PROFILING)
  list (APPEND F90SRC bl_prof.f90)
else()
  list (APPEND F90SRC bl_prof_stubs.f90)
endif()

if (ENABLE_OMP)
  list (APPEND F90SRC omp.f90)
else()
  list (APPEND F90SRC omp_stubs.f90)
endif()

if (ENABLE_PARTICLES)
   list (APPEND F90SRC mt19937ar.f90 particles_f.f90)
endif ()


set ( ALLSRC  ${CXXSRC} ${F90SRC} ${CSRC} )

#
# Add absolute path (needed by target_sources)
# 
prepend ( ALLSRC ${CMAKE_CURRENT_LIST_DIR} )

target_sources ( amrex PRIVATE ${ALLSRC} )

