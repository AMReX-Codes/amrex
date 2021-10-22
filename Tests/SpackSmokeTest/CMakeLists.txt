#
# Test the amrex installation by
# building and running the code
# in Tests/Amr/Advection_AmrCore/
#
cmake_minimum_required(VERSION 3.21)

project(SingleVortex_Test)

enable_language(Fortran)

if(DEFINED CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
endif()

get_filename_component(_base_dir ${CMAKE_CURRENT_LIST_DIR} DIRECTORY)

set(_base_dir  ${_base_dir}/Amr/Advection_AmrLevel)
set(_src_dir   ${_base_dir}/Source)
set(_input_dir ${_base_dir}/Exec/SingleVortex)

message("Base dir: ${_base_dir}")

find_package(AMReX REQUIRED)

if (AMReX_3D_FOUND)
   set(_dim 3)
elseif (AMReX_2D_FOUND)
   set(_dim 2)
else ()
   message(FATAL_ERROR "Cannot find a 2D or 3D version of AMReX")
endif ()


add_executable(SingleVortex)

target_link_libraries(SingleVortex PUBLIC AMReX::amrex)

target_include_directories(SingleVortex
   PUBLIC
   ${_src_dir}
   ${_src_dir}/Src_K/
   ${_input_dir}
   )

target_sources(SingleVortex
   PRIVATE
   ${_src_dir}/LevelBldAdv.cpp
   ${_src_dir}/Tagging_params.cpp
   ${_src_dir}/Kernels.H
   ${_src_dir}/bc_nullfill.cpp
   ${_src_dir}/AmrLevelAdv.cpp
   ${_src_dir}/AmrLevelAdv.H
   ${_src_dir}/Adv.cpp
   ${_src_dir}/main.cpp

   ${_src_dir}/Src_K/tagging_K.H
   ${_src_dir}/Src_K/slope_K.H
   ${_src_dir}/Src_K/flux_${AMReX_SPACEDIM}d_K.H
   ${_src_dir}/Src_K/Adv_K.H

   ${_input_dir}/Prob_Parm.H
   ${_input_dir}/Prob.cpp
   ${_input_dir}/Adv_prob.cpp
   ${_input_dir}/face_velocity_${AMReX_SPACEDIM}d_K.H
   )

#add_custom_command(
#   TARGET install_test
#   POST_BUILD
#   COMMAND ${CMAKE_COMMAND} -E echo "Running test project"
#   COMMAND ${PROJECT_BINARY_DIR}/install_test ${_input_dir}/inputs max_step=1 > out.txt
#   WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
#   )
