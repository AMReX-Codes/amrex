# LOOM backend only: pinned host memory that kernels read and write directly
# stays coherent across synchronizations in every submission mode.
if (NOT AMReX_LOOM)
    return()
endif ()

# The test does not depend on the dimension, so build it once.
list(GET AMReX_SPACEDIM 0 D)
set(_sources main.cpp)

setup_test(${D} _sources _input_files
    BASE_NAME GPU_LoomPinnedVisibility)

unset(_sources)
unset(_input_files)
unset(D)
