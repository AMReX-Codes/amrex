if (WIN32 AND (AMReX_GPU_BACKEND STREQUAL "CUDA"))
    return()
endif()

# This test solves on float MultiFabs set up from Real data, i.e., it mixes
# precisions, which is not supported with AMReX_PRECISION=FP32MP2.
if (AMReX_PRECISION STREQUAL "FP32MP2")
    return()
endif()

foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
       continue()
    endif ()

    set(_sources
       main.cpp
       MyTest.cpp
       initProb.cpp
       MyTestPlotfile.cpp
       MyTest.H
       initProb_K.H)

    set(_input_files inputs)

    setup_test(${D} _sources _input_files)

    unset(_sources)
    unset(_input_files)
endforeach()
