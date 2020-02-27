# Preamble ####################################################################
#
cmake_minimum_required(VERSION 3.14.0)
project(ParallelReduce LANGUAGES CXX CUDA)
# Note: do not enable CUDA here
#   https://amrex-codes.github.io/amrex/docs_html/GPU.html#building-with-cmake


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
set(ENABLE_CUDA         ON CACHE INTERNAL "")
set(ENABLE_TINY_PROFILE ON CACHE INTERNAL "")

if(ENABLE_CUDA)
    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../../../Tools/CMake")
    include(AMReX_SetupCUDA)
endif()

#   here we build AMReX on the fly instead of relying on a pre-build AMReX that
#   we could add via `find_package(AMReX REQUIRED CONFIG CUDA)`
add_subdirectory(
  ${PROJECT_SOURCE_DIR}/../../../
  ${PROJECT_BINARY_DIR}/_deps/amrex
)


# Targets #####################################################################
#
add_executable(parallel_reduce main.cpp)
if(ENABLE_CUDA)
    set_source_files_properties(main.cpp PROPERTIES LANGUAGE CUDA)
endif()

# Required C++ standard
target_compile_features(parallel_reduce PUBLIC cxx_std_11)
set_target_properties(parallel_reduce PROPERTIES
    CXX_EXTENSIONS OFF
    CXX_STANDARD_REQUIRED ON
)

# dependency: AMReX
target_link_libraries(parallel_reduce PRIVATE AMReX::amrex)


# Testing #####################################################################
#
enable_testing()
add_test(NAME parallel_reduce
    COMMAND ${PROJECT_BINARY_DIR}/parallel_reduce
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
)


# Status Message for Build Options ############################################
#
message("")
message("ParallelReduce build configuration:")
message("  C++ Compiler: ${CMAKE_CXX_COMPILER_ID} "
                        "${CMAKE_CXX_COMPILER_VERSION} "
                        "${CMAKE_CXX_COMPILER_WRAPPER}")
message("    ${CMAKE_CXX_COMPILER}")
message("  CUDA Compiler: ${CMAKE_CUDA_COMPILER_ID} "
                         "${CMAKE_CUDA_COMPILER_VERSION} "
                         "${CMAKE_CUDA_COMPILER_WRAPPER}")
message("    ${CMAKE_CUDA_COMPILER}")
message("")
message("  Build Type: ${CMAKE_BUILD_TYPE}")
message("")

