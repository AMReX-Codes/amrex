# -*- mode: cmake -*-

include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   CTMGLIB_SOURCE_DIR
#   CTMGLIB_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(CTMGLIB)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CBOXLIB_INCLUDE_DIRS})

set(CXX_source_files AMReX_DivVis.cpp AMReX_MCCGSolver.cpp AMReX_MCInterpBndryData.cpp AMReX_MCLinOp.cpp AMReX_MCMultiGrid.cpp)
set(FPP_source_files AMReX_MCINTERPBNDRYDATA_${BL_SPACEDIM}D.F AMReX_MCLO_${BL_SPACEDIM}D.F)
if("${BL_SPACEDIM}" EQUAL 2)
  list(APPEND FPP_source_files "AMReX_DV_2D.F")
elseif("${BL_SPACEDIM}" EQUAL 3)
  list(APPEND FPP_source_files "AMReX_DV_3D1.F" "AMReX_DV_3D2.F" "AMReX_DV_3D3.F")
endif()
set(F77_source_files)
set(F90_source_files)

set(CXX_header_files AMReX_DivVis.H AMReX_MCCGSolver.H AMReX_MCInterpBndryData.H AMReX_MCLinOp.H AMReX_MCMultiGrid.H)
set(FPP_header_files AMReX_DivVis_F.H AMReX_MCINTERPBNDRYDATA_F.H AMReX_MCLO_F.H)
set(F77_header_files)
set(F90_header_files)

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

set(local_source_files ${FPP_out_files} ${F77_source_files} ${F90_source_files} ${CXX_source_files})
set(local_header_files ${FPP_header_files} ${F77_header_files} ${F90_header_files} ${CXX_header_files})
add_library(box_ctmg OBJECT ${local_source_files})

install_include_files(${local_header_files})

if (BUILD_TESTS)

endif()

