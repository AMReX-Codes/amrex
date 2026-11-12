# LOOM backend only: code objects load at the first launch of one of their
# kernels, and kernel arguments of every width reach the device.
if (NOT AMReX_LOOM)
    return()
endif ()

# The test does not depend on the dimension, so build it once.
list(GET AMReX_SPACEDIM 0 D)
set(_sources main.cpp)

setup_test(${D} _sources _input_files
    BASE_NAME GPU_LoomCodeObjects)

unset(_sources)
unset(_input_files)
unset(D)
