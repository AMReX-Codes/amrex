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
set (DIM 3 CACHE STRING "Dimension of AMReX build")
if ( (DIM GREATER 3) OR (DIM LESS 1) )
   message( FATAL_ERROR "DIM must be either 1, 2 or 3.")
endif ()
message( STATUS "Building AMReX with DIM = ${DIM}")


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
# variable ENABLE_FORTRAN will be created and set to ON.
# This will stop the subsequent option( ENABLE_FORTRAN "Enable Fortran language" ON )
# from being executed and no entry ENABLE_FORTRAN will be created in the cache
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
# Option to control if Fortran must be enabled  ==============================
#
if ( USE_XSDK_DEFAULTS )
   option( XSDK_ENABLE_Fortran "Enable Fortran language" OFF )
   set( ENABLE_FORTRAN ${XSDK_ENABLE_Fortran} )
   print_option(XSDK_ENABLE_Fortran)
endif ()

option( ENABLE_FORTRAN "Enable Fortran language" ON )
print_option( ENABLE_FORTRAN )

#
# Option to control precision of the build  ==================================
#
if ( USE_XSDK_DEFAULTS )
   set( XSDK_PRECISION "DOUBLE" CACHE STRING "Precision:<SINGLE,DOUBLE>" )
   if ( XSDK_PRECISION STREQUAL "DOUBLE" )
      set( ENABLE_DP ON )
      set( ENABLE_DP_PARTICLES ON )
   elseif ( XSDK_PRECISION STREQUAL "SINGLE" )
      set( ENABLE_DP OFF )
      set( ENABLE_DP_PARTICLES OFF )
   else ()
      message( FATAL_ERROR "\nUnsupported precision ${XSDK_PRECISION}\n" )
   endif ()
   print_option( XSDK_PRECISION )
endif ()

option( ENABLE_DP "Enable double precision" ON )
print_option( ENABLE_DP )



#
# Parallel backends    ========================================================
#

# For the time being ENABLE_DPCPP is defined before project() is called
# Check whether the C++ compiler is dpcpp
print_option(ENABLE_DPCPP)
if (ENABLE_DPCPP AND (NOT (CMAKE_CXX_COMPILER MATCHES "dpcpp") ) )
   message(FATAL_ERROR "\nENABLE_DPCPP=${ENABLE_DPCPP} but CXX compiler is not dpcpp\n")
endif ()

cmake_dependent_option( ENABLE_DPCPP_AOT  "Enable DPCPP ahead-of-time compilation (WIP)"  OFF
   "ENABLE_DPCPP" OFF)
print_option( ENABLE_DPCPP_AOT )

cmake_dependent_option( ENABLE_DPCPP_SPLIT_KERNEL "Enable DPCPP kernel splitting"  ON
   "ENABLE_DPCPP" OFF)
print_option(  ENABLE_DPCPP_SPLIT_KERNEL )

cmake_dependent_option( ENABLE_MPI  "Enable MPI"  ON
   "NOT ENABLE_DPCPP" OFF)
print_option( ENABLE_MPI )

cmake_dependent_option( ENABLE_MPI_THREAD_MULTIPLE
   "whether to initialize MPI so that multiple threads can make MPI calls at the same time"  OFF
   "ENABLE_MPI" OFF)
print_option( ENABLE_MPI_THREAD_MULTIPLE )

option( ENABLE_OMP  "Enable OpenMP" OFF)
print_option( ENABLE_OMP )

cmake_dependent_option( ENABLE_CUDA "Enable GPU support via CUDA" OFF
   "NOT ENABLE_DPCPP" OFF)
print_option( ENABLE_CUDA )

option( ENABLE_ACC  "Enable GPU support via OpenACC" OFF )
print_option( ENABLE_ACC )

if (ENABLE_CUDA OR ENABLE_ACC)
   set(GPUS_PER_SOCKET "IGNORE" CACHE STRING "Number of GPUs per socket" )
   print_option(GPUS_PER_SOCKET)

   set(GPUS_PER_NODE "IGNORE" CACHE STRING "Number of GPUs per node" )
   print_option(GPUS_PER_NODE)
endif ()


#
# AMReX components selection  ================================================
#
option( ENABLE_EB "Build EB Code" OFF )
print_option(ENABLE_EB)

cmake_dependent_option( ENABLE_FORTRAN_INTERFACES "Build Fortran API" OFF
   "ENABLE_FORTRAN" OFF )
print_option(ENABLE_FORTRAN_INTERFACES)

option( ENABLE_LINEAR_SOLVERS  "Build AMReX Linear solvers" ON )
print_option( ENABLE_LINEAR_SOLVERS )

option( ENABLE_AMRDATA "Build data services" OFF)
print_option( ENABLE_AMRDATA )

option( ENABLE_PARTICLES "Build particle classes" OFF)
print_option( ENABLE_PARTICLES )

cmake_dependent_option( ENABLE_DP_PARTICLES "Enable double-precision particle data" ON
   "ENABLE_PARTICLES" OFF )
print_option( ENABLE_DP_PARTICLES )


#
# External packages  =========================================================
#

# sensei
option( ENABLE_SENSEI_INSITU "Enable SENSEI in situ infrastructure" OFF )
print_option( ENABLE_SENSEI_INSITU )

# Conduit (requires CONDUIT_DIR)
option( ENABLE_CONDUIT "Enable Conduit support" OFF )
print_option( ENABLE_CONDUIT )

# Ascent
cmake_dependent_option( ENABLE_ASCENT "Enable Ascent support" OFF
   "ENABLE_CONDUIT" OFF )
print_option( ENABLE_ASCENT )

# SUNDIALS
cmake_dependent_option(ENABLE_SUNDIALS "Enable SUNDIALS4 interfaces" OFF
   "ENABLE_FORTRAN_INTERFACES" OFF)
print_option(ENABLE_SUNDIALS)

# Hypre
cmake_dependent_option(ENABLE_HYPRE "Enable Hypre interfaces" OFF
   "ENABLE_LINEAR_SOLVERS" OFF)
print_option(ENABLE_HYPRE)

# PETSc
cmake_dependent_option(ENABLE_PETSC "Enable PETSc interfaces" OFF
   "ENABLE_LINEAR_SOLVERS" OFF )
print_option(ENABLE_PETSC)

# HDF5
option(ENABLE_HDF5 "Enable HDF5-based I/O" OFF)
print_option(ENABLE_HDF5)

cmake_dependent_option(ENABLE_HDF5_ASYNC "Enable asynchronous writes in the HDF5-based IO" OFF
   "ENABLE_HDF5" OFF )
print_option(ENABLE_HDF5_ASYNC)

if (ENABLE_HDF5_ASYNC)
   message(FATAL_ERROR "\nENABLE_HDF5_ASYNC not yet supported\n")
endif ()


#
# Miscellanoues options  =====================================================
#
option( ENABLE_PIC "Build position-independent code" OFF)
print_option( ENABLE_PIC )

option(ENABLE_FPE "Enable Floating Point Exceptions checks" OFF)
print_option( ENABLE_FPE )

if ( "${CMAKE_BUILD_TYPE}" MATCHES "Debug" )
   option( ENABLE_ASSERTIONS "Enable assertions" ON)
else ()
   option( ENABLE_ASSERTIONS "Enable assertions" OFF)
endif ()

print_option( ENABLE_ASSERTIONS )


#
# Profiling options  =========================================================
#
option( ENABLE_BASE_PROFILE "Enable basic profiling" OFF )
print_option( ENABLE_BASE_PROFILE )

cmake_dependent_option( ENABLE_TINY_PROFILE "Enable tiny profiling" OFF
   "NOT ENABLE_BASE_PROFILE" OFF)
print_option( ENABLE_TINY_PROFILE )

cmake_dependent_option( ENABLE_TRACE_PROFILE "Enable trace-profiling" OFF
   "ENABLE_BASE_PROFILE" OFF)
print_option( ENABLE_TRACE_PROFILE )

option( ENABLE_MEM_PROFILE   "Enable memory profiling" OFF )
print_option( ENABLE_MEM_PROFILE )

cmake_dependent_option( ENABLE_COMM_PROFILE  "Enable communicator-profiling" OFF
   "ENABLE_BASE_PROFILE" OFF)
print_option( ENABLE_COMM_PROFILE )

cmake_dependent_option(ENABLE_PROFPARSER "Enable profile parser" OFF
   "ENABLE_BASE_PROFILE;ENABLE_TRACE_PROFILE;ENABLE_AMRDATA" OFF)
print_option( ENABLE_PROFPARSER )

set(TP_PROFILE_VALUES IGNORE CRAYPAT FORGE VTUNE)
set(TP_PROFILE IGNORE CACHE STRING "Third-party profiling options: <CRAYPAT,FORGE,VTUNE>")
set_property(CACHE TP_PROFILE PROPERTY STRINGS ${TP_PROFILE_VALUES})
if(NOT TP_PROFILE IN_LIST TP_PROFILE_VALUES)
    message(FATAL_ERROR "TP_PROFILE (${TP_PROFILE}) must be one of ${TP_PROFILE_VALUES}")
endif()
print_option( TP_PROFILE )

# Check profile options
if ( NOT ( CMAKE_C_COMPILER_ID STREQUAL "Intel" ) AND
      ( TP_PROFILE STREQUAL "VTUNE") )
   message( FATAL_ERROR "VTUNE cannot be used with ${CMAKE_C_COMPILER_ID} compiler" )
endif ()

if (  ( ( TP_PROFILE STREQUAL "CRAYPAT" ) OR
        ( TP_PROFILE STREQUAL "FORGE"   ) OR
        ( TP_PROFILE STREQUAL "VTUNE"   )   ) AND
     (ENABLE_BASE_PROFILE OR ENABLE_TINY_PROFILE) )
   message(WARNING "This configuration should only be used to profile BL_PROFILE!")
endif()


#
# Extra options  =========================================================
#
option(ALLOW_DIFFERENT_COMPILER
    "Allow an application to use a different compiler than the one used to build AMReX" OFF)
