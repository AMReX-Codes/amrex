# -*- mode: cmake -*-

include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   CFMGLIB_SOURCE_DIR
#   CFMGLIB_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(CFMGLIB)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CBOXLIB_INCLUDE_DIRS})
include_directories(${FBOXLIB_INCLUDE_DIRS})

set(CXX_source_files MGT_Solver.cpp FMultiGrid.cpp)
set(FPP_source_files)
set(F77_source_files)
set(F90_source_files)

set(CXX_header_files MGT_Solver.H FMultiGrid.H stencil_types.H)
set(FPP_header_files)
set(F77_header_files)
set(F90_header_files)


preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

set(local_source_files ${FPP_out_files} ${F77_source_files} ${F90_source_files} ${CXX_source_files})
set(local_header_files ${FPP_header_files} ${F77_header_files} ${F90_header_files} ${CXX_header_files})
add_library(box_cfmg OBJECT ${local_source_files})

add_install_include_file(${local_header_files})

if (BUILD_TESTS)

endif()

