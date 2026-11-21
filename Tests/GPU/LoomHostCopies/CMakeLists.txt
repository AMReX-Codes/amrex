# LOOM backend only: asynchronous copies to and from pageable host memory
# behave as with hipMemcpyAsync, and pinned copies stay asynchronous.
if (NOT AMReX_LOOM)
    return()
endif ()

# The test does not depend on the dimension, so build it once.
list(GET AMReX_SPACEDIM 0 D)
set(_sources main.cpp)

setup_test(${D} _sources _input_files
    BASE_NAME GPU_LoomHostCopies)

unset(_sources)
unset(_input_files)
unset(D)
