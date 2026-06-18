if ( (NOT AMReX_EB) OR (NOT AMReX_LINEAR_SOLVERS_INCFLO) OR NOT (2 IN_LIST AMReX_SPACEDIM))
   return()
endif ()

set(_sources main.cpp)
set(_input_files
    inputs.cart_eb_neumann_sphere
    inputs.cart_mms_rap
    inputs.rz_eb_inflow_sphere
    inputs.rz_eb_neumann_all_neumann
    inputs.rz_eb_neumann_all_neumann_offset
    inputs.rz_eb_neumann_mixed_bc
    inputs.rz_eb_neumann_plane
    inputs.rz_eb_neumann_sphere
    inputs.rz_eb_rap_zero
    inputs.rz_mms
    inputs.rz_mms_rap)

setup_test(2 _sources _input_files BASE_NAME LinearSolvers_NodeRZEBBaseline)

set(_test_exe_dir ${CMAKE_CURRENT_BINARY_DIR}/2d)
set(_test_exe ${_test_exe_dir}/Test_LinearSolvers_NodeRZEBBaseline_2d)
set(_base_test_name LinearSolvers_NodeRZEBBaseline_2d)

foreach(_input_file IN LISTS _input_files)
   if (_input_file STREQUAL "inputs.cart_eb_neumann_sphere")
      continue()
   endif ()

   get_filename_component(_input_name ${_input_file} NAME)
   string(REPLACE "inputs." "" _case_name ${_input_name})
   set(_test_name ${_base_test_name}_${_case_name})

   if (AMReX_OMP)
      add_test(
         NAME               ${_test_name}
         COMMAND            ${_test_exe} ${_input_name}
         WORKING_DIRECTORY  ${_test_exe_dir}
      )
      set_tests_properties(${_test_name} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=2)
   elseif (AMReX_MPI)
      add_test(
         NAME               ${_test_name}
         COMMAND            mpiexec -n 2 ${_test_exe} ${_input_name}
         WORKING_DIRECTORY  ${_test_exe_dir}
      )
   else ()
      add_test(
         NAME               ${_test_name}
         COMMAND            ${_test_exe} ${_input_name}
         WORKING_DIRECTORY  ${_test_exe_dir}
      )
   endif ()
endforeach()

unset(_test_exe_dir)
unset(_test_exe)
unset(_base_test_name)
unset(_input_name)
unset(_case_name)
unset(_test_name)
unset(_sources)
unset(_input_files)
