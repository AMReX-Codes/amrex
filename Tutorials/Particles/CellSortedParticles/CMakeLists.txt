# Preamble ####################################################################
#
cmake_minimum_required(VERSION 3.14.0)
project(CellSortedParticles LANGUAGES CXX Fortran)


# Options and Variants ########################################################
#
# Change Default: Release (instead of Debug)
set(CMAKE_CONFIGURATION_TYPES "Release;Debug;MinSizeRel;RelWithDebInfo")
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the build type, e.g. Release or Debug." FORCE)
endif()


# Dependencies ################################################################
#
# AMReX
#   preset required, non-default AMReX options
#   https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#building-with-cmake
set(ENABLE_PARTICLES    ON CACHE INTERNAL "")
set(ENABLE_TINY_PROFILE ON CACHE INTERNAL "")

#   here we build AMReX on the fly instead of relying on a pre-build AMReX that
#   we could add via `find_package(AMReX REQUIRED CONFIG)`
add_subdirectory(
  ${PROJECT_SOURCE_DIR}/../../../
  ${PROJECT_BINARY_DIR}/_deps/amrex
)


# Targets #####################################################################
#
add_executable(CellSortedParticles
    CellSortedPC.cpp cell_sorted_3d.F90 main.cpp)

# Required C++ standard
target_compile_features(CellSortedParticles PUBLIC cxx_std_11)
set_target_properties(CellSortedParticles PROPERTIES
    CXX_EXTENSIONS OFF
    CXX_STANDARD_REQUIRED ON
)

# dependency: AMReX
target_link_libraries(CellSortedParticles PRIVATE AMReX::amrex)

# copy input file to build directory
file(COPY ${PROJECT_SOURCE_DIR}/inputs
     DESTINATION ${PROJECT_BINARY_DIR})


# Testing #####################################################################
#
enable_testing()
add_test(NAME CellSortedParticles
    COMMAND ${PROJECT_BINARY_DIR}/CellSortedParticles inputs
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
)


# Status Message for Build Options ############################################
#
message("")
message("CellSortedParticles build configuration:")
message("  C++ Compiler: ${CMAKE_CXX_COMPILER_ID} "
                        "${CMAKE_CXX_COMPILER_VERSION} "
                        "${CMAKE_CXX_COMPILER_WRAPPER}")
message("    ${CMAKE_CXX_COMPILER}")
message("  Fortran Compiler: ${CMAKE_Fortran_COMPILER_ID} "
                            "${CMAKE_Fortran_COMPILER_VERSION} "
                            "${CMAKE_Fortran_COMPILER_WRAPPER}")
message("    ${CMAKE_Fortran_COMPILER}")
message("")
message("  Build Type: ${CMAKE_BUILD_TYPE}")
message("")

