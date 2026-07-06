# This test requires particle support
if (NOT AMReX_PARTICLES)
   return()
endif ()

# We can enable it again, once the test is fixed for GPU.
if (NOT AMReX_GPU_BACKEND STREQUAL NONE)
   return()
endif ()

foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources main.cpp)

    setup_test(${D} _sources _input_files)

    unset(_sources)
endforeach()
