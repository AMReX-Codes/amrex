# This test requires particle support
if (NOT AMReX_PARTICLES)
   return()
endif ()

foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources main.cpp)

    setup_test(${D} _sources _input_files)

    unset(_sources)
endforeach()
