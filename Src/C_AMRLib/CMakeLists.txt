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

set(CXX_source_files Amr.cpp AmrLevel.cpp AuxBoundaryData.cpp BCRec.cpp Cluster.cpp DatasetClient.cpp Derive.cpp ErrorList.cpp FluxRegister.cpp Interpolater.cpp SlabStat.cpp StateData.cpp StateDescriptor.cpp StationData.cpp TagBox.cpp)
set(FPP_source_files ARRAYLIM_${BL_SPACEDIM}D.F FILCC_${BL_SPACEDIM}D.F FLUXREG_${BL_SPACEDIM}D.F INTERP_${BL_SPACEDIM}D.F SLABSTAT_${BL_SPACEDIM}D.F)
if(BL_SPACEDIM EQUAL 3)
  set(FPP_source_files ${FPP_source_files} MAKESLICE_${BL_SPACEDIM}D.F)
endif()

set(CXX_header_files Amr.H AmrLevel.H AuxBoundaryData.H BCRec.H BC_TYPES.H Cluster.H DatasetClient.H Derive.H ErrorList.H FluxRegister.H INTERP_F.H Interpolater.H LevelBld.H PROB_AMR_F.H SLABSTAT_F.H SlabStat.H StateData.H StateDescriptor.H StationData.H TagBox.H)
set(FPP_header_files FLUSH_F.H FLUXREG_F.H)
if(BL_SPACEDIM EQUAL 3)
  list(APPEND FPP_header_files MAKESLICE_F.H)
endif()
set(F77_header_files bc_types.fi)
set(F90_header_files)

if(BL_USE_PARTICLES EQUAL 1)
  list(APPEND CXX_source_files Particles.cpp)
  list(APPEND CXX_header_files Particles.H)
  list(APPEND FPP_source_files Particles_${BL_SPACEDIM}D.F)
  list(APPEND FPP_header_files Particles_F.H)
endif()
set(F77_source_files)
set(F90_source_files)

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

set(local_source_files ${FPP_out_files} ${F77_source_files} ${F90_source_files} ${CXX_source_files})
set(local_header_files ${CXX_header_files} ${F77_header_files} ${FPP_header_files} ${F90_header_files})
add_library(box_camr OBJECT ${local_source_files})

add_install_include_file(${local_header_files})

if (BUILD_TESTS)

endif()
