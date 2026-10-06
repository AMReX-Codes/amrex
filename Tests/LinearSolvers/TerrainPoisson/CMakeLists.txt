if (NOT (3 IN_LIST AMReX_SPACEDIM))
    return()
endif ()

set(_sources main.cpp)
set(_input_files inputs)

setup_test(3 _sources _input_files)

# Boxes split in z, with each line solver, on the same executable.
set(_exe_dir ${CMAKE_CURRENT_BINARY_DIR}/3d)
foreach (_solver IN ITEMS spike column)
    set(_name LinearSolvers_TerrainPoisson_zsplit_${_solver}_3d)
    set(_cmd ${_exe_dir}/Test_LinearSolvers_TerrainPoisson_3d inputs
        "max_grid_size_v=16 16 8" zsplit_solver=${_solver})
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
unset(_name)
unset(_cmd)

unset(_sources)
unset(_input_files)
