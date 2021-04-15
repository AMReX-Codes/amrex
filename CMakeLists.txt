cmake_minimum_required(VERSION 3.14)

#
# Allow for MSVC Runtime library controls
#
if(POLICY CMP0091)
  cmake_policy(SET CMP0091 NEW)
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
include( AMReX_Utils )
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
    # CMake 3.18+: CMAKE_CUDA_ARCHITECTURES
    # https://cmake.org/cmake/help/latest/policy/CMP0104.html
    if(POLICY CMP0104)
        cmake_policy(SET CMP0104 OLD)
    endif()

    enable_language(CUDA)
    include(AMReX_SetupCUDA)
endif ()

#
# Check compiler version
#
set_mininum_cxx_compiler_version(GNU 4.8)
set_mininum_cxx_compiler_version(MSVC 19.23)

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
if(AMReX_INSTALL)
    include(AMReXInstallHelpers)
    install_amrex_targets(${_amrex_targets})

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
