if (AMReX_MPI)
   return()
endif()

set(_sources     main.cpp)
set(_input_files inputs  )

setup_test(_sources _input_files)

unset(_sources)
unset(_input_files)
