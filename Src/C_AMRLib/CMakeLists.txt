# -*- mode: cmake -*-

include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   CAMRLIB_SOURCE_DIR
#   CAMRLIB_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(CAMRLIB)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CBOXLIB_INCLUDE_DIRS})

set(CXX_source_files AMReX_Amr.cpp AMReX_AmrLevel.cpp AMReX_AuxBoundaryData.cpp AMReX_Derive.cpp AMReX_Extrapolater.cpp AMReX_SlabStat.cpp AMReX_StateData.cpp AMReX_StateDescriptor.cpp AMReX_StationData.cpp )
set(FPP_source_files AMReX_ARRAYLIM_${BL_SPACEDIM}D.F AMReX_SLABSTAT_${BL_SPACEDIM}D.F)
if(BL_SPACEDIM EQUAL 3)
  set(FPP_source_files ${FPP_source_files} AMReX_MAKESLICE_${BL_SPACEDIM}D.F)
endif()

set(CXX_header_files AMReX_Amr.H AMReX_AmrLevel.H AMReX_AuxBoundaryData.H AMReX_Derive.H AMReX_Extrapolater.H AMReX_LevelBld.H AMReX_PROB_AMR_F.H AMReX_SLABSTAT_F.H AMReX_SlabStat.H AMReX_StateData.H AMReX_StateDescriptor.H AMReX_StationData.H )
set(FPP_header_files AMReX_FLUSH_F.H )
if(BL_SPACEDIM EQUAL 3)
  list(APPEND FPP_header_files AMReX_MAKESLICE_F.H)
endif()
set(F77_header_files)
set(F90_header_files)

set(F77_source_files)
set(F90_source_files AMReX_extrapolater_${BL_SPACEDIM}d.f90)

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

set(local_source_files ${FPP_out_files} ${F77_source_files} ${F90_source_files} ${CXX_source_files})
set(local_header_files ${CXX_header_files} ${F77_header_files} ${FPP_header_files} ${F90_header_files})
add_library(box_camr OBJECT ${local_source_files})

add_install_include_file(${local_header_files})

if (BUILD_TESTS)

endif()
