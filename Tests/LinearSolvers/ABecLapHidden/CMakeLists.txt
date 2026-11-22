if (NOT (3 IN_LIST AMReX_SPACEDIM))
    return()
endif ()

set(_sources main.cpp)
set(_input_files inputs)

setup_test(3 _sources _input_files)

# More cases on the same executable: each hidden direction, a singular
# operator, Jacobi smoothing, and semicoarsening along z and along x, whose
# line solve is host only.
set(_exe_dir ${CMAKE_CURRENT_BINARY_DIR}/3d)
set(_cases hidden0 hidden2 singular jacobi)
if (AMReX_GPU_BACKEND STREQUAL NONE)
    list(APPEND _cases semicoarsening semicoarsening_x)
endif ()
foreach (_case IN LISTS _cases)
    set(_name LinearSolvers_ABecLapHidden_${_case}_3d)
    if (_case STREQUAL "hidden0")
        set(_args hidden_direction=0 "n_cell=1 32 32" "max_grid_size=1 16 16")
    elseif (_case STREQUAL "hidden2")
        set(_args hidden_direction=2 "n_cell=32 32 1" "max_grid_size=16 16 1"
                  "bc_lo=dirichlet neumann neumann" "bc_hi=dirichlet dirichlet neumann")
    elseif (_case STREQUAL "singular")
        set(_args ascalar=0 "bc_lo=periodic neumann neumann" "bc_hi=periodic neumann neumann")
    elseif (_case STREQUAL "jacobi")
        set(_args use_gauss_seidel=0)
    elseif (_case STREQUAL "semicoarsening")
        set(_args semicoarsening=1 semicoarsening_direction=2 max_semicoarsening_level=3
                  "prob_hi=1.0 1.0 0.125")
    else ()
        set(_args semicoarsening=1 semicoarsening_direction=0 max_semicoarsening_level=3
                  "prob_hi=0.125 1.0 1.0")
    endif ()
    set(_cmd ${_exe_dir}/Test_LinearSolvers_ABecLapHidden_3d inputs ${_args})
    if (AMReX_OMP)
        add_test(NAME ${_name} COMMAND ${_cmd} WORKING_DIRECTORY ${_exe_dir})
        set_tests_properties(${_name} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=2)
    elseif (AMReX_MPI)
        add_test(NAME ${_name} COMMAND mpiexec -n 2 ${_cmd} WORKING_DIRECTORY ${_exe_dir})
    else ()
        add_test(NAME ${_name} COMMAND ${_cmd} WORKING_DIRECTORY ${_exe_dir})
    endif ()
endforeach ()
unset(_exe_dir)
unset(_cases)
unset(_name)
unset(_args)
unset(_cmd)

unset(_sources)
unset(_input_files)
