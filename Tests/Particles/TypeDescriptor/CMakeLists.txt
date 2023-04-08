if (AMReX_MPI)
   return ()
endif()

set(_sources main.cpp)

setup_test(_sources FALSE)

unset(_sources)
