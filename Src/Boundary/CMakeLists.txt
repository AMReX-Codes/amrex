# -*- mode: cmake -*-

include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   CBNDRYLIB_SOURCE_DIR
#   CBNDRYLIB_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(CBNDRYLIB)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CBOXLIB_INCLUDE_DIRS})

set(CXX_source_files AMReX_BndryData.cpp AMReX_BndryRegister.cpp AMReX_FabSet.cpp AMReX_InterpBndryData.cpp AMReX_MacBndry.cpp AMReX_Mask.cpp AMReX_MultiMask.cpp)
set(FPP_source_files AMReX_INTERPBNDRYDATA_${BL_SPACEDIM}D.F AMReX_LO_UTIL.F )
set(F77_source_files)
set(F90_source_files)

set(CXX_header_files AMReX_BndryData.H AMReX_BndryRegister.H AMReX_BoundCond.H AMReX_FabSet.H AMReX_InterpBndryData.H AMReX_LO_BCTYPES.H AMReX_MacBndry.H AMReX_Mask.H AMReX_MultiMask.H)
set(FPP_header_files AMReX_INTERPBNDRYDATA_F.H)
set(F77_header_files)
set(F90_header_files)

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})
set(local_source_files ${CXX_source_files} ${F77_source_files} ${FPP_out_files} ${F90_source_files})
set(local_header_files ${CXX_header_files} ${F77_header_files} ${FPP_header_files} ${F90_header_files})
add_library(box_cbndry OBJECT ${local_source_files})

add_dependencies(box_cbndry box_c)


install_include_files(${local_header_files})

if (BUILD_TESTS)

endif()

