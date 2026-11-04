if (AMReX_MPI)
   return ()
endif()

# The test does not depend on the dimension, so build it once.
list(GET AMReX_SPACEDIM 0 D)
set(_sources main.cpp)

setup_test(${D} _sources FALSE)

unset(_sources)
unset(D)
