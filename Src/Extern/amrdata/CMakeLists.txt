# -*- mode: cmake -*-

include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   AMRDATALIB_SOURCE_DIR
#   AMRDATALIB_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(AMRDATALIB)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CBOXLIB_INCLUDE_DIRS})

set(CXX_source_files AmrData.cpp DataServices.cpp)
set(FPP_source_files FABUTIL_${BL_SPACEDIM}D.F)
set(F77_source_files)
set(F90_source_files)

set(CXX_header_files AmrData.H AmrvisConstants.H DataServices.H)
set(FPP_header_files)
set(F77_header_files)
set(F90_header_files)

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})
set(local_source_files ${CXX_source_files} ${F77_source_files} ${FPP_out_files} ${F90_source_files})
set(local_header_files ${CXX_header_files} ${F77_header_files} ${FPP_header_files} ${F90_header_files})
add_library(box_camrdata ${local_source_files})

add_install_include_file(${local_header_files})
add_install_library(box_camrdata)

if (BUILD_TESTS)

endif()

