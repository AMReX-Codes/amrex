# -*- mode: cmake -*-

include(TestManager)

#
# Define a project name
# After this command the following variables are defined
#   FBOXLIB_SOURCE_DIR
#   FBOXLIB_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(FBOXLIB)

include_directories(${FBOXLIB_INCLUDE_DIRS})

set(CXX_source_files)
set(C_source_files f2kgetcl.c  fabio_c.c  ppm_util_c.c  system_util_c.c  timer_c.c)
set(F77_source_files)
set(FPP_source_files)
set(F90_source_files bc.f90 bl_constants.f90 bl_error.f90 bl_IO.f90 bl_mem_stat.f90 bl_parmparse.f90 bl_prof.f90 bl_prof_stubs.f90 bl_space.f90 bl_stream.f90 bl_string.f90 bl_system.f90 bl_timer.f90 bl_types.f90 bndry_reg.f90 boxarray_f.f90 box_f.f90 boxlib_f.f90 box_util.f90 cluster_f.f90 create_umac_grown.f90 cutcells.f90 define_bc_tower.f90 f2kcli_crayx1.f90 f2kcli.f90 f2kcli_nag.f90 f2kcli_win32.f90 fab.f90 fabio.f90 filler.f90 fillpatch.f90 fourth_order_interp_coeffs.f90 interp.f90 knapsack.f90 layout.f90 list_box.f90 make_new_grids.f90 ml_boxarray.f90 ml_layout.f90 ml_multifab.f90 mt19937ar.f90 multifab_f.f90 multifab_fill_ghost_cells.f90 multifab_physbc_edgevel.f90 multifab_physbc.f90 omp.f90 omp_stubs.f90 particles_f.f90 pingpong.f90 plotfile.f90 ppm_util.f90 sort_box.f90 sort_d.f90 sort_i.f90 tag_boxes.f90 vector_i.f90)

if (ENABLE_MPI)
  list(APPEND F90_source_files parallel.f90)
else()
  list(APPEND F90_source_files parallel_stubs.f90)
endif()

set(CXX_header_files)
set(F77_header_files)
set(FPP_header_files)
set(F90_header_files)

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

set(local_source_files ${C_source_files} ${CXX_source_files} ${F77_source_files} ${FPP_out_files} ${F90_source_files})
set(local_header_files ${CXX_header_files} ${F77_header_files} ${FPP_header_files} ${F90_header_files})
add_library(box_f ${local_source_files})

add_install_include_file(${local_header_files})
add_install_library(box_f)

if (BUILD_TESTS AND BL_SPACEDIM EQUAL 3)
endif()
   
