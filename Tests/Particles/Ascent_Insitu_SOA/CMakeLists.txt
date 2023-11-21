if ( NOT AMReX_ASCENT )
   return ()
endif ()

foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources     main.cpp)
    set(_input_files inputs.rt  )

    setup_test(${D} _sources _input_files NTASKS 2)

    unset(_sources)
    unset(_input_files)
endforeach()
