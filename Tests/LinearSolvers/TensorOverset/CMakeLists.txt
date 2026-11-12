# The driver is 3D and double precision only.
if (NOT AMReX_LINEAR_SOLVERS_INCFLO OR NOT 3 IN_LIST AMReX_SPACEDIM
    OR AMReX_PRECISION STREQUAL "SINGLE")
    return()
endif ()

set(_sources main.cpp MyTest.cpp MyTest.H MyTest_K.H)
set(_input_files)
set(_ctest_params n_cell=32 max_grid_size=16 max_error=4.e-2 do_plots=0)

setup_test(3 _sources _input_files
   CMDLINE_PARAMS ${_ctest_params})

# Single precision multigrid cycles, with and without the overset mask; the
# error must match the double precision run.
setup_test(3 _sources _input_files
   BASE_NAME LinearSolvers_TensorOverset_mixed
   RUNTIME_SUBDIR mixed
   CMDLINE_PARAMS ${_ctest_params} mixed_precision=1)

setup_test(3 _sources _input_files
   BASE_NAME LinearSolvers_TensorOverset_mixed_nomask
   RUNTIME_SUBDIR mixed_nomask
   CMDLINE_PARAMS ${_ctest_params} mixed_precision=1 do_overset=0)

unset(_sources)
unset(_input_files)
unset(_ctest_params)
