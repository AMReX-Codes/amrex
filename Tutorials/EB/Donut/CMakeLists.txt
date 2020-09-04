if (NOT DIM EQUAL 3)
   return()
endif ()

set(_sources Src/main.cpp)
set(_input_files Exec/inputs)

setup_tutorial(_sources _input_files)

unset(_sources)
unset(_input_files)
