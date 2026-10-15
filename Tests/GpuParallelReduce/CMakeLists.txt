foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources     main.cpp)
    set(_input_files)

    setup_test(${D} _sources _input_files)

    unset(_sources)
    unset(_input_files)
endforeach()

if (AMReX_NCCL_DEVICE_API AND AMReX_MPI AND AMReX_GPU_BACKEND STREQUAL CUDA)
    find_path(AMReX_NCCL_INCLUDE_DIR
        NAMES nccl.h
        HINTS ${AMReX_NCCL_ROOT}
        PATH_SUFFIXES include)
    find_path(AMReX_NCCL_DEVICE_INCLUDE_DIR
        NAMES nccl_device.h
        HINTS ${AMReX_NCCL_ROOT}
        PATH_SUFFIXES include)
    find_library(AMReX_NCCL_LIBRARY
        NAMES nccl
        HINTS ${AMReX_NCCL_ROOT}
        PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu)

    if (NOT AMReX_NCCL_INCLUDE_DIR OR
        NOT AMReX_NCCL_DEVICE_INCLUDE_DIR OR
        NOT AMReX_NCCL_LIBRARY)
        message(FATAL_ERROR
            "AMReX_NCCL_DEVICE_API requires nccl.h, nccl_device.h, and libnccl")
    endif ()

    foreach(D IN LISTS AMReX_SPACEDIM)
        set(_target Test_GpuParallelReduce_NCCLDevice_${D}d)
        add_executable(${_target} nccl_device_api.cu)
        target_include_directories(${_target} PRIVATE
            ${AMReX_NCCL_INCLUDE_DIR}
            ${AMReX_NCCL_DEVICE_INCLUDE_DIR})
        target_link_libraries(${_target} PRIVATE
            AMReX::amrex_${D}d
            ${AMReX_NCCL_LIBRARY})
        set_target_properties(${_target} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${D}d
            CUDA_SEPARABLE_COMPILATION ON)
        setup_target_for_cuda_compilation(${_target})

        unset(_target)
    endforeach()
endif ()
