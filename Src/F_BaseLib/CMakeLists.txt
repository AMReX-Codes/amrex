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
set(C_source_files fabio_c.c timer_c.c)
set(F77_source_files)
set(FPP_source_files)
set(F90_source_files bc.f90 bl_constants.f90 bl_error.f90 bl_IO.f90 bl_stream.f90 bl_string.f90 bl_timer.f90 bl_types.f90 bl_mem_stat.f90 bl_space.f90 bl_prof_stubs.f90 bndry_reg.f90 box_f.f90 boxarray_f.f90 fab.f90 fabio.f90 knapsack.f90 layout.f90 list_box.f90 ml_boxarray.f90 ml_layout.f90 ml_multifab.f90 multifab_f.f90 sort_box.f90 sort_d.f90 sort_i.f90 vector_i.f90)

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

install(TARGETS box_f
  LIBRARY DESTINATION ${CCSE_INSTALL_DIR}/lib${MPI_SUFFIX}${OMP_SUFFIX}
  ARCHIVE DESTINATION ${CCSE_INSTALL_DIR}/lib${MPI_SUFFIX}${OMP_SUFFIX}
  CONFIGURATIONS ${CMAKE_BUILD_TYPE}
  )

if (${CCSE_INSTALL_DIR} STREQUAL ${PROJECT_SOURCE_DIR})
else ()
    install(FILES ${local_header_files}
        DESTINATION ${CCSE_INSTALL_DIR}/include)
endif ()

if (BUILD_TESTS AND BL_SPACEDIM EQUAL 3)
endif()
   
