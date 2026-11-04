# The test does not depend on the dimension, so build it once.
list(GET AMReX_SPACEDIM 0 D)
set(_sources     main.cpp)
set(_input_files inputs  )

setup_test(${D} _sources _input_files)

unset(_sources)
unset(_input_files)
unset(D)
