# This turorials requires some EB functionalities which are active
# only if Fortran is enabled
if (NOT DIM EQUAL 3 OR NOT CMAKE_Fortran_COMPILER_LOADED)
   return()
endif ()

set(_sources Src/main.cpp)
set(_input_files Exec/inputs)

setup_tutorial(_sources _input_files)

unset(_sources)
unset(_input_files)
