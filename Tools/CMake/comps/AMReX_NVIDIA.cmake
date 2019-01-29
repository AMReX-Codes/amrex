#
# Setup for NVIDIA compile, AKA NVCC
# This file should be included by 'setup_amrex_compilers()' only,
# AFTER the setup for ALL the other supported compilers have been loaded.
#
# target_compile_option() for CUDA does not work yet as expected since it
# does not propagate the NVCC flags to the device linking phase.
# Consequently we have to use CMAKE_CUDA_FLAGS to make sure that all the
# CUDA compilation/linking steps use the selected flags.
# Hopefully this behavior will be fixed in future CMake versions.
# All the NVCC flags are to be considered REQUIRED for the application code.
#
string( REPLACE "." ";" VERSION_LIST ${CMAKE_CUDA_COMPILER_VERSION})
list( GET VERSION_LIST 0 NVCC_VERSION_MAJOR )
list( GET VERSION_LIST 1 NVCC_VERSION_MINOR )

target_compile_definitions( amrex PUBLIC
   AMREX_NVCC_VERSION=${CMAKE_CUDA_COMPILER_VERSION}
   AMREX_NVCC_MAJOR_VERSION=${NVCC_VERSION_MAJOR}
   AMREX_NVCC_MINOR_VERSION=${NVCC_VERSION_MINOR} )

if (CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "8.0")
   message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION}."
      "This is unsupported. Please use CUDA toolkit version 8.0 or newer.")
endif ()

# string(REPLACE ";" " " NVCC_ARCH_FLAGS "${NVCC_ARCH_FLAGS}")
# set(cuda_flags "--expt-relaxed-constexpr --expt-extended-lambda --std=c++11 -dc" )
# set(cuda_flags "${cuda_flags} -Wno-deprecated-gpu-targets -m64 ${NVCC_ARCH_FLAGS} -maxrregcount=${CUDA_MAXREGCOUNT}")
# set(cuda_flags "${cuda_flags} -Xcompiler=--std=c++11")

# if (ENABLE_CUDA_FASTMATH)
#    set(cuda_flags "${cuda_flags} --use_fast_math")
# endif ()

target_compile_options(amrex PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr --expt-extended-lambda>)
target_compile_options(amrex PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:${cuda_flags} -Wno-deprecated-gpu-targets -m64 ${NVCC_ARCH_FLAGS} -maxrregcount=${CUDA_MAXREGCOUNT} -Xcompiler=--std=c++11 --use_fast_math>)

# set( CMAKE_CUDA_FLAGS "${cuda_flags}" CACHE
#    STRING "Flags used by the CUDA compiler during all build types.")

# set( CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS_RELEASE} -lineinfo --ptxas-options=-O3,-v" CACHE
#    STRING "Flags used by the CUDA compiler during RELEASE builds.")

# set( CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -G" CACHE
#    STRING "Flags used by the CUDA compiler during DEBUG builds.")

print(CMAKE_CUDA_SEPARABLE_COMPILATION)
# set_property(TARGET amrex PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
set_target_properties(amrex
   PROPERTIES
   CUDA_SEPARABLE_COMPILATION ON  # This adds -dc
   CUDA_STANDARD 11               # Adds -std=c++11
   CUDA_STANDARD_REQUIRED ON
   CUDA_RESOLVE_DEVICE_SYMBOLS OFF
   )


