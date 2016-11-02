# -*- mode: cmake -*-

include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   CPARTLIB_SOURCE_DIR
#   CPARTLIB_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(CPARTLIB)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CBOXLIB_INCLUDE_DIRS})

set(CXX_source_files)
set(FPP_source_files)
set(F77_source_files)
set(F90_source_files)

set(CXX_header_files)
set(FPP_header_files)
set(F77_header_files)
set(F90_header_files)

if(BL_USE_PARTICLES EQUAL 1)
  list(APPEND CXX_source_files Particles.cpp TracerParticles.cpp)
  list(APPEND CXX_header_files Particles.H ParticleInit.H TracerParticles.H ParGDB.H)
  list(APPEND FPP_source_files Particles_${BL_SPACEDIM}D.F)
  list(APPEND FPP_header_files Particles_F.H)
endif()

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})
set(local_source_files ${CXX_source_files} ${F77_source_files} ${FPP_out_files} ${F90_source_files})
set(local_header_files ${CXX_header_files} ${F77_header_files} ${FPP_header_files} ${F90_header_files})
add_library(box_cpart OBJECT ${local_source_files})

add_install_include_file(${local_header_files})

if (BUILD_TESTS)

endif()

