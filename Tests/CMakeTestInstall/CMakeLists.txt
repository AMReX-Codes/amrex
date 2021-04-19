#
# Test the amrex installation by
# building and running the code
# in Tests/Amr/Advection_AmrCore/
#
cmake_minimum_required(VERSION 3.14)

project(amrex-test-install)

get_filename_component(_base_dir ${CMAKE_CURRENT_LIST_DIR} DIRECTORY)

set(_base_dir  ${_base_dir}/Amr/Advection_AmrCore)
set(_src_dir   ${_base_dir}/Source)
set(_input_dir ${_base_dir}/Exec)


find_package(AMReX REQUIRED)

if (AMReX_3D_FOUND)
   set(_dim 3)
elseif (AMReX_2D_FOUND)
   set(_dim 2)
else ()
   message(FATAL_ERROR "Cannot find a 2D or 3D version of AMReX")
endif ()


add_executable(install_test)

target_link_libraries(install_test PUBLIC AMReX::amrex)

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

add_custom_command(
   TARGET install_test
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} -E echo "Running test project"
   COMMAND ${PROJECT_BINARY_DIR}/install_test ${_input_dir}/inputs max_step=1 > out.txt
   WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
   )
