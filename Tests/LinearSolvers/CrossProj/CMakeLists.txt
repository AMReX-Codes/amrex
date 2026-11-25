foreach(D IN LISTS AMReX_SPACEDIM)
    if (NOT D EQUAL 2)
        continue()
    endif ()

    set(_sources main.cpp)
    set(_input_files)

    setup_test(${D} _sources _input_files)

    # Box lengths 2 mod 4: the reduced level has boxes at odd indices
    setup_test(${D} _sources _input_files
       BASE_NAME LinearSolvers_CrossProj_OddBoxes
       CMDLINE_PARAMS n_cell_x=72 n_cell_y=72 max_grid_size=2 max_grid_size_2=6)

    unset(_sources)
    unset(_input_files)
endforeach()
