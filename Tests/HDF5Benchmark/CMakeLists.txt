string(REGEX REPLACE "Tests.*" "" _base_name ${CMAKE_CURRENT_LIST_DIR})

# This tests requires particle support
if (NOT AMReX_PARTICLES)
   return()
endif ()

set(_sources     main.cpp)
set(_input_files inputs  )

setup_test(_sources _input_files)

target_include_directories(Test_HDF5Benchmark PRIVATE ${_base_name}/Src/Extern/SWFFT)

unset(_sources)
unset(_input_files)
unset(_base_name)
