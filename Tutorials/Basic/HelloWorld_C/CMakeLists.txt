# Preamble ####################################################################
#
cmake_minimum_required(VERSION 3.14.0)
project(HelloWorld LANGUAGES CXX)


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
#   here we build AMReX on the fly instead of relying on a pre-build AMReX that
#   we could add via `find_package(AMReX REQUIRED CONFIG)`
add_subdirectory(
  ${PROJECT_SOURCE_DIR}/../../../
  ${PROJECT_BINARY_DIR}/_deps/amrex
)


# Targets #####################################################################
#
add_executable(hello_amrex main.cpp)

# Required C++ standard
target_compile_features(hello_amrex PUBLIC cxx_std_11)
set_target_properties(hello_amrex PROPERTIES
    CXX_EXTENSIONS OFF
    CXX_STANDARD_REQUIRED ON
)

# dependency: AMReX
target_link_libraries(hello_amrex PRIVATE AMReX::amrex)


# Testing #####################################################################
#
enable_testing()
add_test(NAME hello_amrex
    COMMAND ${PROJECT_BINARY_DIR}/hello_amrex
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
message("")
message("  Build Type: ${CMAKE_BUILD_TYPE}")
message("")

