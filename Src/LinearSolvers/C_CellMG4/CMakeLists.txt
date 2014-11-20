# -*- mode: cmake -*-

include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   CMGLIB_SOURCE_DIR
#   CMGLIB_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(CMGLIB)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CBOXLIB_INCLUDE_DIRS})

set(CXX_source_files ABecLaplacian.cpp CGSolver.cpp Laplacian.cpp LinOp.cpp MultiGrid.cpp)
set(FPP_source_files ABec_${BL_SPACEDIM}D.F ABec_UTIL.F CG_${BL_SPACEDIM}D.F LO_${BL_SPACEDIM}D.F LP_${BL_SPACEDIM}D.F MG_${BL_SPACEDIM}D.F)
set(F77_source_files)
set(F90_source_files)

set(CXX_header_files ABecLaplacian.H CGSolver.H Laplacian.H LinOp.H MultiGrid.H)
set(FPP_header_files ABec_F.H CG_F.H LO_F.H LP_F.H MG_F.H)
set(F77_header_files lo_bctypes.fi)
set(F90_header_files)



preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

set(local_source_files ${FPP_out_files} ${F77_source_files} ${F90_source_files} ${CXX_source_files})
set(local_header_files ${FPP_header_files} ${F77_header_files} ${F90_header_files} ${CXX_header_files})
add_library(box_cmg OBJECT ${local_source_files})

add_install_include_file(${local_header_files})

if (BUILD_TESTS)

endif()

