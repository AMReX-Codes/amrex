#
# Test the amrex installation by
# building and running the code
# in Tests/Amr/Advection_AmrCore/
#
cmake_minimum_required(VERSION 3.18)

project(amrex-test-install)

if(DEFINED CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
endif()

get_filename_component(_base_dir ${CMAKE_CURRENT_LIST_DIR} DIRECTORY)

set(_base_dir  ${_base_dir}/Amr/Advection_AmrCore)
set(_src_dir   ${_base_dir}/Source)
set(_input_dir ${_base_dir}/Exec)


find_package(AMReX REQUIRED)

foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
       continue()
    endif ()

    add_executable(install_test_${D}d)

    target_link_libraries(install_test_${D}d PUBLIC AMReX::amrex_${D}d)

    target_include_directories(install_test_${D}d
       PUBLIC
       ${_src_dir}
       ${_src_dir}/Src_K/
       ${_input_dir}
    )

    target_sources(install_test_${D}d
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

    file(TO_NATIVE_PATH "${_input_dir}/inputs" inputs_path)

    add_custom_command(
       TARGET install_test_${D}d
       POST_BUILD
       COMMAND ${CMAKE_COMMAND} -E echo "Running test project"
       COMMAND ${CMAKE_COMMAND} -E echo "Command: $<TARGET_FILE:install_test_${D}d> ${inputs_path} max_step=1 > out_${D}d.txt"
       COMMAND $<TARGET_FILE:install_test_${D}d> ${inputs_path} max_step=1 > out_${D}d.txt
       WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    )
endforeach()
