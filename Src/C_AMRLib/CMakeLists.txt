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

set(CXX_source_files Amr.cpp AmrLevel.cpp AuxBoundaryData.cpp Cluster.cpp Derive.cpp ErrorList.cpp FillPatchUtil.cpp FluxRegister.cpp Interpolater.cpp Extrapolater.cpp SlabStat.cpp StateData.cpp StateDescriptor.cpp StationData.cpp TagBox.cpp)
set(FPP_source_files ARRAYLIM_${BL_SPACEDIM}D.F FILCC_${BL_SPACEDIM}D.F FLUXREG_${BL_SPACEDIM}D.F INTERP_${BL_SPACEDIM}D.F SLABSTAT_${BL_SPACEDIM}D.F)
if(BL_SPACEDIM EQUAL 3)
  set(FPP_source_files ${FPP_source_files} MAKESLICE_${BL_SPACEDIM}D.F)
endif()

set(CXX_header_files Amr.H AmrLevel.H AuxBoundaryData.H Cluster.H Derive.H ErrorList.H FillPatchUtil.H FluxRegister.H INTERP_F.H Interpolater.H Extrapolater.H LevelBld.H PROB_AMR_F.H SLABSTAT_F.H SlabStat.H StateData.H StateDescriptor.H StationData.H TagBox.H)
set(FPP_header_files FLUSH_F.H FLUXREG_F.H)
if(BL_SPACEDIM EQUAL 3)
  list(APPEND FPP_header_files MAKESLICE_F.H)
endif()
set(F77_header_files bc_types.fi)
set(F90_header_files)

set(F77_source_files)
set(F90_source_files extrapolater_${BL_SPACEDIM}d.f90)

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

set(local_source_files ${FPP_out_files} ${F77_source_files} ${F90_source_files} ${CXX_source_files})
set(local_header_files ${CXX_header_files} ${F77_header_files} ${FPP_header_files} ${F90_header_files})
add_library(box_camr OBJECT ${local_source_files})

add_install_include_file(${local_header_files})

if (BUILD_TESTS)

endif()
