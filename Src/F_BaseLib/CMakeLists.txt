set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")


list ( APPEND ALLHEADERS )
list ( APPEND F90SRC     boxlib_f.f90 f2kcli.f90 BLProfiler_f90.f90 )
list ( APPEND F90SRC     bl_constants.f90 bl_IO.f90 )

#
# Profiling. In the GNUMake system, this is enabled if the legacy
# flag PROF is TRUE. Here, I enabled profiling if ENABLE_BASIC_PROFILE is ON.
# 
if ( ENABLE_BASIC_PROFILE )
   list ( APPEND F90SRC     bl_prof.f90 amrex_timer_c.f90 )
   list ( APPEND CXXSRC     timer_c.c )
else ()
   list ( APPEND F90SRC     bl_prof_stubs.f90 )
endif ()

list ( APPEND F90SRC     bl_error.f90 bl_mem_stat.f90 bl_parmparse.f90 bl_space.f90 )
list ( APPEND F90SRC     bl_stream.f90 bl_string.f90 bl_system.f90 bl_timer.f90 )
list ( APPEND F90SRC     amrex_timer_c.f90 bl_types.f90 box_f.f90 knapsack.f90 )
list ( APPEND F90SRC     layout.f90 boxarray_f.f90 box_util.f90 fab.f90 multifab_f.f90 )
list ( APPEND F90SRC     fabio.f90 amrex_fabio_c.f90 plotfile.f90 filler.f90 cluster_f.f90 )
list ( APPEND F90SRC     interp.f90 fourth_order_interp_coeffs.f90  bc.f90 bc_functions.f90 )
list ( APPEND F90SRC     bndry_reg.f90 ml_boxarray.f90 ml_layout.f90 ml_multifab.f90 )
list ( APPEND F90SRC     ml_restrict_fill.f90 cc_restriction.f90 nodal_restriction.f90 )
list ( APPEND F90SRC     ml_cc_restriction.f90 ml_nd_restriction.f90 nodal_neumann_bcs.f90 )
list ( APPEND F90SRC     nodal_stencil_bc.f90 create_umac_grown.f90 define_bc_tower.f90 )
list ( APPEND F90SRC     fillpatch.f90 multifab_fill_ghost_cells.f90 multifab_physbc_edgevel.f90 )
list ( APPEND F90SRC     multifab_physbc.f90 list_box.f90 sort_box.f90 vector_i.f90 )
list ( APPEND F90SRC     sort_d.f90 sort_i.f90 ppm_util.f90 cutcells.f90 regrid.f90 )
list ( APPEND F90SRC     make_new_grids.f90 tag_boxes.f90 )

if ( ENABLE_PARTICLES )
   list ( APPEND F90SRC     particles_f.f90 mt19937ar.f90 )
endif ()

if ( ENABLE_MPI )
   list ( APPEND F90SRC     parallel.f90 )
else ()
   list ( APPEND F90SRC     parallel_stubs.f90 )
endif ()

if (ENABLE_OMP)
  list ( APPEND  F90SRC     omp.f90 )
else()
  list ( APPEND  F90SRC     omp_stubs.f90 )
endif()


list ( APPEND CXXSRC     fabio_c.c timer_c.c ppm_util_c.c system_util_c.c backtrace_c.cpp )
list ( APPEND F90SRC     backtrace_f.f90 ) # DO NOT INCLUDE THIS ONE: CONFLICT amrex_fort_mod.f90 )

# IF CXX11 (always in AMReX)
list ( APPEND CXXSRC     bl_random_c.cpp )
list ( APPEND F90SRC     bl_random_f.f90 )


#
# Collect all source files
# 
set ( ALLSRC  ${CXXSRC} ${F90SRC} ${CSRC} )

#
# Add absolute path (needed by target_sources)
# 
prepend ( ALLSRC ${CMAKE_CURRENT_LIST_DIR} )

target_sources ( amrex PRIVATE ${ALLSRC} )

#
# Add to include dirs
#
accumulate_include_dirs (${CMAKE_CURRENT_LIST_DIR})

