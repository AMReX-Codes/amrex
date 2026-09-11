# The hypre comparison is skipped at run time when hypre is not available.
if ( (NOT AMReX_EB) OR (NOT AMReX_LINEAR_SOLVERS_EM) )
   return()
endif ()

foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
       continue()
    endif ()

    set(_sources main.cpp)
    set(_input_files inputs_${D}d)

    setup_test(${D} _sources _input_files)

    unset(_sources)
    unset(_input_files)
endforeach()
