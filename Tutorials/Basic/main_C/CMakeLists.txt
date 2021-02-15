if (NOT AMReX_MPI)
   return()
endif ()

set(_sources main.cpp)
set(_input_files)

setup_tutorial(_sources  _input_files)

unset(_sources)
unset(_input_files)
