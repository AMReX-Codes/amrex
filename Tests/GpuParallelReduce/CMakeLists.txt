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

        set(_gin_target Test_GpuParallelReduce_NCCLGIN_${D}d)
        add_executable(${_gin_target} nccl_gin.cpp nccl_gin.cu)
        target_include_directories(${_gin_target} PRIVATE
            ${AMReX_NCCL_INCLUDE_DIR} ${AMReX_NCCL_DEVICE_INCLUDE_DIR})
        target_link_libraries(${_gin_target} PRIVATE
            AMReX::amrex_${D}d ${AMReX_NCCL_LIBRARY})
        setup_target_for_cuda_compilation(${_gin_target})
        # Keep the numerical reference on the host compiler, separate from CUDA.
        set_source_files_properties(nccl_gin.cpp PROPERTIES LANGUAGE CXX)
        set_target_properties(${_gin_target} PROPERTIES
            CXX_STANDARD 23 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF
            CUDA_STANDARD 20 CUDA_STANDARD_REQUIRED ON
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${D}d)
        target_compile_options(${_gin_target} PRIVATE
            "$<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>:-fno-fast-math>")
        add_test(NAME GpuParallelReduce_NCCLGIN_${D}d
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${_gin_target}>
                ${MPIEXEC_POSTFLAGS})
        set_tests_properties(GpuParallelReduce_NCCLGIN_${D}d PROPERTIES
            TIMEOUT 120 PROCESSORS 2 LABELS "GPU;NCCL;multinode")
        unset(_gin_target)
    endforeach()
endif ()
