cmake_minimum_required(VERSION 3.18)

#
# Allow for MSVC Runtime library controls
#
if(POLICY CMP0091)
  cmake_policy(SET CMP0091 NEW)
endif()

#
# We use simple syntax in cmake_dependent_option, so we are compatible with the
# extended syntax in CMake 3.22+
# https://cmake.org/cmake/help/v3.22/policy/CMP0127.html
#
if(POLICY CMP0127)
    cmake_policy(SET CMP0127 NEW)
endif()

#
# CMake 3.18+: CMAKE_CUDA_ARCHITECTURES
# https://cmake.org/cmake/help/latest/policy/CMP0104.html
# We have to migrate there, but maybe the new "native" option (CMake 3.24+)
# is what we want to wait for:
# https://cmake.org/cmake/help/v3.24/prop_tgt/CUDA_ARCHITECTURES.html
if(POLICY CMP0104)
    cmake_policy(SET CMP0104 OLD)
endif()

#
# Prevent in-source builds
#
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
   message(FATAL_ERROR
      "\nin-source builds are not allowed: "
      "build directory cannot be in the source directory path!\n"
      "You MUST remove the file ${CMAKE_BINARY_DIR}/CMakeCache.txt and "
      " the directory ${CMAKE_BINARY_DIR}/CMakeFiles/ to be able to build again.")
endif ()


#
# Set search path for AMReX-specific CMake modules
#
set( AMREX_CMAKE_MODULES_PATH "${CMAKE_CURRENT_LIST_DIR}/Tools/CMake" CACHE INTERNAL "" )
set( CMAKE_MODULE_PATH ${AMREX_CMAKE_MODULES_PATH} )

#
# Retrieve amrex version
#
include( AMReXUtils )
get_amrex_version()


########################################################################
#
# AMReX project
#
########################################################################
project( AMReX
    DESCRIPTION "A software framework for massively parallel, block-structured adaptive mesh refinement (AMR) applications"
    VERSION  ${AMREX_PKG_VERSION}
    HOMEPAGE_URL "https://amrex-codes.github.io/amrex/"
    LANGUAGES C CXX
    )

message(STATUS "CMake version: ${CMAKE_VERSION}")


#
# Provide a default install directory
#
if ( CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR AND CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
    set ( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/installdir"
          CACHE PATH "AMReX installation directory" FORCE)
endif ()

message(STATUS "AMReX installation directory: ${CMAKE_INSTALL_PREFIX}")

#
# Check if CMAKE_BUILD_TYPE is given. If not, use default
#
if ( NOT CMAKE_BUILD_TYPE )
   set(CMAKE_CONFIGURATION_TYPES "Release;Debug;MinSizeRel;RelWithDebInfo")
   set(CMAKE_BUILD_TYPE Release
       CACHE STRING
       "Choose the build type, e.g. Release, Debug, or RelWithDebInfo." FORCE)
else ()
   message(STATUS "Build type set by user to '${CMAKE_BUILD_TYPE}'.")
endif()

#
# Include options, utilities and other stuff we need
#
include( AMReXOptions )
if (AMReX_GPU_BACKEND STREQUAL "CUDA")
   include(AMReXCUDAOptions)
endif ()

#
# Enable Fortran if requested
#
if(AMReX_FORTRAN)
   enable_language(Fortran)
endif ()

#
# Enable CUDA if requested
#
if (AMReX_CUDA)
    enable_language(CUDA)
    if(CMAKE_VERSION VERSION_LESS 3.20)
        include(AMReX_SetupCUDA)
    endif()
endif ()

#
# Check compiler version
#
set_mininum_compiler_version(CXX GNU 5.1)
set_mininum_compiler_version(CXX MSVC 19.23)

#
# Set CMAKE_<LANG>_FLAGS_<CONFIG> if not already defined
#
set_default_config_flags ()

#
# Source files for all binaries and libraries found under src
#
add_subdirectory(Src)
get_property(_amrex_targets
   DIRECTORY Src
   PROPERTY BUILDSYSTEM_TARGETS)

#
# Plotfile tools
#
option(AMReX_PLOTFILE_TOOLS "Enable Plotfile tools" NO)

if (AMReX_PLOTFILE_TOOLS)
   # If this get executed, it cannot be EXCLUDED_FROM_ALL
   # because it needs to get installed
   add_subdirectory(Tools/Plotfile)
   get_property(_plotfile_targets
      DIRECTORY Tools/Plotfile
      PROPERTY BUILDSYSTEM_TARGETS)
   list(APPEND _amrex_targets ${_plotfile_targets})
endif ()


#
# Install amrex  -- Export
#
include(AMReXInstallHelpers)
install_amrex_targets(${_amrex_targets})

if(AMReX_INSTALL)
    # Add a test_install target to smoke-test
    # the installation
    add_test_install_target(
       ${CMAKE_CURRENT_LIST_DIR}/Tests/CMakeTestInstall
       ${CMAKE_INSTALL_PREFIX})
endif()

#
# Enable CTests
#
option(AMReX_ENABLE_TESTS "Enable CTest suite for AMReX"  NO)
if (AMReX_ENABLE_TESTS)
   enable_testing()
   add_subdirectory(Tests)
endif ()
