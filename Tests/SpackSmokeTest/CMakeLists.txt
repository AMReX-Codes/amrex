# This is a stand-alone test developed for Spack. In Spack,
# this CMakeLists.txt file is used to build the test when the
# command, `spack test run amrex` is called.
#
# This CMake script builds the test, Tests/Amr/Advection_AmrCore,
# against a currently installed version of AMReX. The resulting
# executable can then be ran to test functionality.

cmake_minimum_required(VERSION 3.18)

project(amrex-test-install)

get_filename_component(_base_dir ${CMAKE_CURRENT_LIST_DIR} DIRECTORY)

set(_base_dir  ${_base_dir}/Amr/Advection_AmrCore)
set(_src_dir   ${_base_dir}/Source)
set(_input_dir ${_base_dir}/Exec)


if (AMReX_FORTRAN_INTERFACES)
    enable_language(Fortran)
endif()

if (AMReX_GPU_BACKEND STREQUAL "CUDA")
    enable_language(CUDA)
    find_package(AMReX REQUIRED CUDA)

elseif(AMReX_GPU_BACKEND STREQUAL "HIP")
    if(CMAKE_VERSION VERSION_LESS 3.20)
        message(FATAL_ERROR "HIP requires CMake version 3.20 or newer")
    endif()
    find_package(AMReX REQUIRED HIP)
    find_package(rocrand REQUIRED CONFIG)
    find_package(rocprim REQUIRED CONFIG)
    find_package(hiprand REQUIRED CONFIG)

else()
   # Use installed version of AMReX
   find_package(AMReX REQUIRED)
endif()

if (AMReX_3D_FOUND)
    set(_dim 3)
elseif (AMReX_2D_FOUND)
    set(_dim 2)
else ()
    message(FATAL_ERROR "Cannot find a 2D or 3D version of AMReX")
endif ()


add_executable(install_test)

target_include_directories(install_test
    PUBLIC
    ${_src_dir}
    ${_src_dir}/Src_K/
    ${_input_dir}
    )

target_sources(install_test
    PRIVATE
    ${_src_dir}/AdvancePhiAllLevels.cpp
    ${_src_dir}/AdvancePhiAtLevel.cpp
    ${_src_dir}/AmrCoreAdv.cpp
    ${_src_dir}/AmrCoreAdv.H
    ${_src_dir}/bc_fill.H
    ${_src_dir}/DefineVelocity.cpp
    ${_src_dir}/face_velocity.H
    ${_src_dir}/Kernels.H
    ${_src_dir}/main.cpp
    ${_src_dir}/Tagging.H
    ${_src_dir}/Src_K/Adv_K.H
    ${_src_dir}/Src_K/compute_flux_2D_K.H
    ${_src_dir}/Src_K/slope_K.H
    ${_input_dir}/Prob.H
    )

target_link_libraries(install_test PUBLIC AMReX::amrex)

# Additional CUDA configuration commands
if (AMReX_GPU_BACKEND STREQUAL "CUDA")
    set_target_properties(install_test PROPERTIES
        LANGUAGE CUDA
        CUDA_ARCHITECTURES "${AMREX_CUDA_ARCHS}"
        CUDA_SEPARABLE_COMPILATION  ON)

    get_target_property(Source_Files install_test SOURCES)
    list(FILTER Source_Files INCLUDE REGEX "\\.cpp$")
    set_source_files_properties(${Source_Files} PROPERTIES LANGUAGE CUDA)
endif()
