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

set(CXX_source_files AMReX_ABecLaplacian.cpp AMReX_CGSolver.cpp AMReX_Laplacian.cpp AMReX_LinOp.cpp AMReX_MultiGrid.cpp)
set(FPP_source_files AMReX_ABec_${BL_SPACEDIM}D.F AMReX_ABec_UTIL.F AMReX_LO_${BL_SPACEDIM}D.F AMReX_LP_${BL_SPACEDIM}D.F AMReX_MG_${BL_SPACEDIM}D.F)
set(F77_source_files)
set(F90_source_files)

set(CXX_header_files AMReX_ABecLaplacian.H AMReX_CGSolver.H AMReX_Laplacian.H AMReX_LinOp.H AMReX_MultiGrid.H)
set(FPP_header_files AMReX_ABec_F.H AMReX_LO_F.H AMReX_LP_F.H AMReX_MG_F.H)
set(F77_header_files AMReX_lo_bctypes.fi)
set(F90_header_files)



preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

set(local_source_files ${FPP_out_files} ${F77_source_files} ${F90_source_files} ${CXX_source_files})
set(local_header_files ${FPP_header_files} ${F77_header_files} ${F90_header_files} ${CXX_header_files})
add_library(box_cmg OBJECT ${local_source_files})

install_include_files(${local_header_files})

if (BUILD_TESTS)

endif()

