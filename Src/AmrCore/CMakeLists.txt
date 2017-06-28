include (TestManager)

project (CAMRCORELIB)

include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories (${CBOXLIB_INCLUDE_DIRS})

set ( CXX_source_files
   AMReX_AmrCore.cpp AMReX_Cluster.cpp AMReX_ErrorList.cpp
   AMReX_FillPatchUtil.cpp AMReX_FluxRegister.cpp
   AMReX_Interpolater.cpp AMReX_TagBox.cpp AMReX_AmrMesh.cpp )

set ( FPP_source_files
   AMReX_FLUXREG_${BL_SPACEDIM}D.F AMReX_INTERP_${BL_SPACEDIM}D.F )

set ( F90PP_source_files AMReX_FillPatchUtil_${BL_SPACEDIM}d.F90 )

set ( CXX_header_files
   AMReX_AmrCore.H AMReX_Cluster.H AMReX_ErrorList.H
   AMReX_FillPatchUtil.H AMReX_FluxRegister.H AMReX_Interpolater.H
   AMReX_TagBox.H AMReX_AmrMesh.H AMReX_INTERP_F.H
   AMReX_FillPatchUtil_F.H)

set ( FPP_header_files AMReX_FLUXREG_F.H AMReX_FillPatchUtil_F.H )
set ( F77_header_files )
set ( F90_header_files )

set ( F77_source_files )

# These files are in Particle
if (ENABLE_PARTICLES EQUAL 1)
  list ( APPEND CXX_header_files AMReX_AmrParGDB.H AMReX_AmrParticles.H )
endif ()

set_F77_properties ( FPP_out_files ${FPP_source_files} )

set ( local_source_files
   ${FPP_out_files} ${F77_source_files} ${F90_source_files}
   ${CXX_source_files} )

set ( local_header_files
   ${CXX_header_files} ${F77_header_files} ${FPP_header_files}
   ${F90_header_files} )

add_library (box_camrcore OBJECT ${local_source_files})

install_include_files ( ${local_header_files} )

if (BUILD_TESTS)

endif()
