###############################################

# Here we define the default config options   #
# that can be overwritten by the user         #

###############################################
include_guard(GLOBAL)

include(CMakeDependentOption)

#
# Define a macro to print active options
#
macro (print_option _var)
   if (${_var})
      message( STATUS "   ${_var}")
   endif ()
endmacro ()


#
# Dimensionality of the build  ===============================================
#
set(AMReX_SPACEDIM_VALUES 1 2 3)
set(AMReX_SPACEDIM 3 CACHE STRING "Dimension of AMReX build: <1,2,3>")
set_property(CACHE AMReX_SPACEDIM PROPERTY STRINGS ${AMReX_SPACEDIM_VALUES})
if(NOT AMReX_SPACEDIM IN_LIST AMReX_SPACEDIM_VALUES)
   message(FATAL_ERROR "AMReX_SPACEDIM=${AMReX_SPACEDIM} is not allowed."
      " Must be one of ${AMReX_SPACEDIM_VALUES}")
endif()
message( STATUS "Building AMReX with AMReX_SPACEDIM = ${AMReX_SPACEDIM}")

#
# Populate the cache and check the value of the user-definable options
#
message(STATUS "Configuring AMReX with the following options enabled: ")

#
# This is the option to enable/disable xSDK mode
#
# To handle both XDSK options and correponding plain AMReX options,
# we make use of policy CMP0077 introduced as a default in CMake 3.13
# Under policy CMP0077, normal variables prevent option()
# to set internal variables of the same name.
# Example: if XSDK mode is ON and XSDK_ENABLE_Fortran=ON, a normal
# variable AMReX_FORTRAN will be created and set to ON.
# This will stop the subsequent option( AMReX_FORTRAN "Enable Fortran language" ON )
# from being executed and no entry AMReX_FORTRAN will be created in the cache
#
option( USE_XSDK_DEFAULTS "Enable xSDK mode"  OFF )
print_option( USE_XSDK_DEFAULTS )


#
# Option to control the type of library: static vs shared  ===================
#
if ( USE_XSDK_DEFAULTS )
   option( BUILD_SHARED_LIBS "Build AMReX shared library" ON )
else ()
   option( BUILD_SHARED_LIBS "Build AMReX shared library" OFF )
endif ()
print_option( BUILD_SHARED_LIBS )

#
# Option to control generation of install targets
#
option(AMReX_INSTALL "Generate Install Targets" ON)


#
# Option to control if Fortran must be enabled  ==============================
#
if ( USE_XSDK_DEFAULTS )
   option( XSDK_ENABLE_Fortran "Enable Fortran language" OFF )
   set( AMReX_FORTRAN ${XSDK_ENABLE_Fortran} )
   print_option(XSDK_ENABLE_Fortran)
endif ()

option( AMReX_FORTRAN "Enable Fortran language" OFF )
print_option( AMReX_FORTRAN )

#
# Option to control precision of the build  ==================================
#
if ( USE_XSDK_DEFAULTS )
   set( XSDK_PRECISION "DOUBLE" CACHE STRING "Precision:<SINGLE,DOUBLE>" )
   set( AMReX_PRECISION ${XSDK_PRECISION})
   set( AMReX_PARTICLES_PRECISION ${XSDK_PRECISION})
   message( STATUS "   XSDK_PRECISION = ${XSDK_PRECISION}")
endif ()


set(AMReX_PRECISION_VALUES SINGLE DOUBLE)
set(AMReX_PRECISION DOUBLE CACHE STRING "Precision of AMReX build: <SINGLE,DOUBLE>")
set_property(CACHE AMReX_PRECISION PROPERTY STRINGS ${AMReX_PRECISION_VALUES})
if(NOT AMReX_PRECISION IN_LIST AMReX_PRECISION_VALUES)
   message(FATAL_ERROR "AMReX_PRECISION=${AMReX_PRECISION} not supported."
      " Must be one of ${AMReX_PRECISION_VALUES}")
endif()
message( STATUS "   AMReX_PRECISION = ${AMReX_PRECISION}")


#
# GPU backends    =============================================================
#
set(AMReX_GPU_BACKEND_VALUES NONE SYCL CUDA HIP)
set(AMReX_GPU_BACKEND NONE CACHE STRING "On-node, accelerated GPU backend: <NONE,SYCL,CUDA,HIP>")
set_property(CACHE AMReX_GPU_BACKEND PROPERTY STRINGS ${AMReX_GPU_BACKEND_VALUES})
if (NOT AMReX_GPU_BACKEND IN_LIST AMReX_GPU_BACKEND_VALUES)
   message(FATAL_ERROR "AMReX_GPU_BACKEND=${AMReX_GPU_BACKEND} is not allowed."
      " Must be one of ${AMReX_GPU_BACKEND_VALUES}")
endif ()

if (NOT AMReX_GPU_BACKEND STREQUAL NONE)
   message( STATUS "   AMReX_GPU_BACKEND = ${AMReX_GPU_BACKEND}")
endif ()

# Legacy variables for internal use only
if (AMReX_GPU_BACKEND STREQUAL SYCL)
   set(AMReX_DPCPP ON )
   set(AMReX_CUDA  OFF)
   set(AMReX_HIP   OFF)
elseif (AMReX_GPU_BACKEND STREQUAL CUDA)
   set(AMReX_DPCPP OFF)
   set(AMReX_CUDA  ON )
   set(AMReX_HIP   OFF)
elseif (AMReX_GPU_BACKEND STREQUAL HIP)
   set(AMReX_DPCPP OFF)
   set(AMReX_CUDA  OFF)
   set(AMReX_HIP   ON )
else ()
   set(AMReX_DPCPP OFF)
   set(AMReX_CUDA  OFF)
   set(AMReX_HIP   OFF)
endif ()

# --- SYCL ---
if (AMReX_DPCPP)
   set(_valid_dpcpp_compiler_ids Clang IntelClang IntelDPCPP)
   if (NOT (CMAKE_CXX_COMPILER_ID IN_LIST _valid_dpcpp_compiler_ids) )
      message(WARNING "\nAMReX_GPU_BACKEND=${AMReX_GPU_BACKEND} is tested with "
         "DPCPP. Verify '${CMAKE_CXX_COMPILER_ID}' is correct and potentially "
         "set CMAKE_CXX_COMPILER=dpcpp.")
   endif ()
   unset(_valid_dpcpp_compiler_ids)
endif ()

cmake_dependent_option( AMReX_DPCPP_AOT  "Enable DPCPP ahead-of-time compilation (WIP)"  OFF
   "AMReX_GPU_BACKEND STREQUAL SYCL" OFF)
print_option( AMReX_DPCPP_AOT )

cmake_dependent_option( AMReX_DPCPP_SPLIT_KERNEL "Enable DPCPP kernel splitting"  ON
   "AMReX_GPU_BACKEND STREQUAL SYCL" OFF)
print_option(  AMReX_DPCPP_SPLIT_KERNEL )

# --- HIP ----
if (AMReX_HIP)
   set(AMReX_AMD_ARCH "IGNORE" CACHE STRING
      "AMD GPU architecture (Must be provided if AMReX_HIP=ON)")
   if (NOT AMReX_AMD_ARCH)
      message(FATAL_ERROR "\nMust specify AMReX_AMD_ARCH if AMReX_HIP=ON\n")
   endif ()
endif ()

if (AMReX_CUDA OR AMReX_HIP)
   set(GPUS_PER_SOCKET "IGNORE" CACHE STRING "Number of GPUs per socket" )
   print_option(GPUS_PER_SOCKET)

   set(GPUS_PER_NODE "IGNORE" CACHE STRING "Number of GPUs per node" )
   print_option(GPUS_PER_NODE)
endif ()


#
# Parallel backends    ========================================================
#
option( AMReX_MPI  "Enable MPI"  ON )
print_option( AMReX_MPI )

cmake_dependent_option( AMReX_MPI_THREAD_MULTIPLE
   "whether to initialize MPI so that multiple threads can make MPI calls at the same time"  OFF
   "AMReX_MPI" OFF)
print_option( AMReX_MPI_THREAD_MULTIPLE )

option( AMReX_OMP  "Enable OpenMP" OFF)
print_option( AMReX_OMP )


#
# AMReX components selection  ================================================
#
option( AMReX_AMRLEVEL  "Build AmrLevel class" ON )
print_option( AMReX_AMRLEVEL )

cmake_dependent_option( AMReX_EB "Build with Embedded Boundary support" OFF
   "NOT AMReX_SPACEDIM EQUAL 1" OFF )
print_option(AMReX_EB)

cmake_dependent_option( AMReX_FORTRAN_INTERFACES "Build Fortran API" OFF
   "AMReX_FORTRAN" OFF )
print_option(AMReX_FORTRAN_INTERFACES)

option( AMReX_LINEAR_SOLVERS  "Build AMReX Linear solvers" ON )
print_option( AMReX_LINEAR_SOLVERS )

cmake_dependent_option( AMReX_AMRDATA "Build data services" OFF
   "AMReX_FORTRAN" OFF )
print_option( AMReX_AMRDATA )

option( AMReX_PARTICLES "Build particle classes" OFF)
print_option( AMReX_PARTICLES )

if (AMReX_PARTICLES)
   set(AMReX_PARTICLES_PRECISION_VALUES SINGLE DOUBLE)
   set(AMReX_PARTICLES_PRECISION ${AMReX_PRECISION}
      CACHE STRING "Precision of reals in particle classes: <SINGLE,DOUBLE>")
   set_property(CACHE AMReX_PARTICLES_PRECISION PROPERTY STRINGS ${AMReX_PARTICLES_PRECISION_VALUES})
   if(NOT AMReX_PARTICLES_PRECISION IN_LIST AMReX_PARTICLES_PRECISION_VALUES)
      message(FATAL_ERROR "AMReX_PARTICLES_PRECISION=${AMReX_PRECISION} not supported."
         " Must be one of ${AMReX_PARTICLES_PRECISION_VALUES}")
   endif()
   message( STATUS "   AMReX_PARTICLES_PRECISION = ${AMReX_PARTICLES_PRECISION}")
endif ()


#
# External packages  =========================================================
#

# sensei
option( AMReX_SENSEI "Enable SENSEI in situ infrastructure" OFF )
print_option( AMReX_SENSEI )

# Conduit (requires CONDUIT_DIR)
option( AMReX_CONDUIT "Enable Conduit support" OFF )
print_option( AMReX_CONDUIT )

# Ascent
cmake_dependent_option( AMReX_ASCENT "Enable Ascent support" OFF
   "AMReX_CONDUIT" OFF )
print_option( AMReX_ASCENT )

# Hypre
cmake_dependent_option(AMReX_HYPRE "Enable Hypre interfaces" OFF
   "AMReX_LINEAR_SOLVERS" OFF)
print_option(AMReX_HYPRE)

# PETSc
cmake_dependent_option(AMReX_PETSC "Enable PETSc interfaces" OFF
   "AMReX_LINEAR_SOLVERS" OFF )
print_option(AMReX_PETSC)

# HDF5
option(AMReX_HDF5 "Enable HDF5-based I/O" OFF)
print_option(AMReX_HDF5)

cmake_dependent_option(AMReX_HDF5_ASYNC "Enable asynchronous writes in the HDF5-based IO" OFF
   "AMReX_HDF5" OFF )
print_option(AMReX_HDF5_ASYNC)

if (AMReX_HDF5_ASYNC)
   message(FATAL_ERROR "\nAMReX_HDF5_ASYNC not yet supported\n")
endif ()


#
# Miscellanoues options  =====================================================
#
option( AMReX_PIC "Build position-independent code" OFF)
print_option( AMReX_PIC )

option( AMReX_IPO "Enable interprocedural optimization (IPO/LTO)" OFF)
print_option( AMReX_IPO )

option(AMReX_FPE "Enable Floating Point Exceptions checks" OFF)
print_option( AMReX_FPE )

if ( "${CMAKE_BUILD_TYPE}" MATCHES "Debug" )
   option( AMReX_ASSERTIONS "Enable assertions" ON)
else ()
   option( AMReX_ASSERTIONS "Enable assertions" OFF)
endif ()

print_option( AMReX_ASSERTIONS )

option(AMReX_BOUND_CHECK  "Enable bound checking in Array4 class" OFF)
print_option( AMReX_BOUND_CHECK )

#
# Profiling options  =========================================================
#
option( AMReX_BASE_PROFILE "Enable basic profiling" OFF )
print_option( AMReX_BASE_PROFILE )

cmake_dependent_option( AMReX_TINY_PROFILE "Enable tiny profiling" OFF
   "NOT AMReX_BASE_PROFILE" OFF)
print_option( AMReX_TINY_PROFILE )

cmake_dependent_option( AMReX_TRACE_PROFILE "Enable trace-profiling" OFF
   "AMReX_BASE_PROFILE" OFF)
print_option( AMReX_TRACE_PROFILE )

option( AMReX_MEM_PROFILE   "Enable memory profiling" OFF )
print_option( AMReX_MEM_PROFILE )

cmake_dependent_option( AMReX_COMM_PROFILE  "Enable communicator-profiling" OFF
   "AMReX_BASE_PROFILE" OFF)
print_option( AMReX_COMM_PROFILE )

cmake_dependent_option(AMReX_PROFPARSER "Enable profile parser" OFF
   "AMReX_BASE_PROFILE;AMReX_TRACE_PROFILE;AMReX_AMRDATA" OFF)
print_option( AMReX_PROFPARSER )

set(AMReX_TP_PROFILE_VALUES IGNORE CRAYPAT FORGE VTUNE)
set(AMReX_TP_PROFILE IGNORE CACHE STRING "Third-party profiling options: <CRAYPAT,FORGE,VTUNE>")
set_property(CACHE AMReX_TP_PROFILE PROPERTY STRINGS ${AMReX_TP_PROFILE_VALUES})
if (NOT AMReX_TP_PROFILE IN_LIST AMReX_TP_PROFILE_VALUES)
   message(FATAL_ERROR "AMReX_TP_PROFILE (${AMReX_TP_PROFILE}) must be one of ${AMReX_TP_PROFILE_VALUES}")
endif()
if (AMReX_TP_PROFILE)
   message( STATUS "   AMReX_TP_PROFILE = ${AMReX_TP_PROFILE}")
endif ()

# Check profile options
if ( NOT ( CMAKE_C_COMPILER_ID STREQUAL "Intel" ) AND
      ( AMReX_TP_PROFILE STREQUAL "VTUNE") )
   message( FATAL_ERROR "VTUNE cannot be used with ${CMAKE_C_COMPILER_ID} compiler" )
endif ()

if (  ( ( AMReX_TP_PROFILE STREQUAL "CRAYPAT" ) OR
        ( AMReX_TP_PROFILE STREQUAL "FORGE"   ) OR
        ( AMReX_TP_PROFILE STREQUAL "VTUNE"   )   ) AND
     (AMReX_BASE_PROFILE OR AMReX_TINY_PROFILE) )
   message(WARNING "This configuration should only be used to profile BL_PROFILE!")
endif()


#
# Extra options  =========================================================
#
option(AMReX_DIFFERENT_COMPILER
   "Allow an application to use a different compiler than the one used to build AMReX" OFF)
print_option(AMReX_DIFFERENT_COMPILER)

if (BUILD_SHARED_LIBS AND NOT (CMAKE_SYSTEM_NAME STREQUAL "Linux") )
   option(AMReX_PROBINIT "Enable support for probin file" OFF)
else ()
   option(AMReX_PROBINIT "Enable support for probin file" ON)
endif ()
print_option(AMReX_PROBINIT)
