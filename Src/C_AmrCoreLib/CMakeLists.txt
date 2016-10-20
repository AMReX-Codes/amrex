# -*- mode: cmake -*-

include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   CAMRCORELIB_SOURCE_DIR
#   CAMRCORELIB_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(CAMRCORELIB)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CBOXLIB_INCLUDE_DIRS})

set(CXX_source_files AmrCore.cpp Cluster.cpp ErrorList.cpp FillPatchUtil.cpp FluxRegister.cpp Interpolater.cpp TagBox.cpp )
set(FPP_source_files FLUXREG_${BL_SPACEDIM}D.F INTERP_${BL_SPACEDIM}D.F )

set(CXX_header_files AmrCore.H Cluster.H ErrorList.H FillPatchUtil.H FluxRegister.H Interpolater.H TagBox.H INTERP_F.H)
set(FPP_header_files FLUXREG_F.H )
set(F77_header_files)
set(F90_header_files)

set(F77_source_files)

# These files are in C_ParticleLib
if(BL_USE_PARTICLES EQUAL 1)
  list(APPEND CXX_header_files AmrParGDB.H AmrParticles.H)
endif()

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

set(local_source_files ${FPP_out_files} ${F77_source_files} ${F90_source_files} ${CXX_source_files})
set(local_header_files ${CXX_header_files} ${F77_header_files} ${FPP_header_files} ${F90_header_files})
add_library(box_camrcore OBJECT ${local_source_files})

add_install_include_file(${local_header_files})

if (BUILD_TESTS)

endif()
