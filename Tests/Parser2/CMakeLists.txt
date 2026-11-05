# The test does not depend on the dimension, so build it once.
list(GET AMReX_SPACEDIM 0 D)
set(_sources     main.cpp fn.cpp)
set(_input_files fn.cpp)

setup_test(${D} _sources _input_files)

unset(_sources)
unset(_input_files)
unset(D)
