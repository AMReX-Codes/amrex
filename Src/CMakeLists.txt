#
# Fist, define the library we want to add
# Think of this as a constructor for object "amrex"
# In the following, we will use setters to change
# the properties of this object, like setting the sources,
# setting the compile definitions and so on
#
add_library( amrex )

# Where to store Fortran modules
set_target_properties( amrex
   PROPERTIES
   Fortran_MODULE_DIRECTORY
   ${PROJECT_BINARY_DIR}/mod_files
   INTERFACE_INCLUDE_DIRECTORIES
   $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/mod_files>
   )

# Load Flags targets and use them if no user-defined flags is given
include(AMReXFlagsTargets)

if (NOT CMAKE_Fortran_FLAGS)
   target_link_libraries(amrex
      PUBLIC
      $<BUILD_INTERFACE:Flags_Fortran>
      )
endif ()

if (NOT CMAKE_CXX_FLAGS)
   target_link_libraries(amrex
      PUBLIC
      $<BUILD_INTERFACE:Flags_CXX>
      )
endif ()

if (ENABLE_FPE)
   target_link_libraries(amrex
      PUBLIC
      $<BUILD_INTERFACE:Flags_FPE>
      )
endif ()

# General configuration
include ( AMReX_Config )
configure_amrex ()

#
# Core components
#
add_subdirectory(Base)
add_subdirectory(Boundary)
add_subdirectory(AmrCore)
add_subdirectory(Amr)

#
# Optional components
#
if (ENABLE_EB)
   if (ENABLE_CUDA AND ( "${NVCC_VERSION_MAJOR}.${NVCC_VERSION_MINOR}" VERSION_EQUAL "9.2") )
      message(FATAL_ERROR "EB component of AMReX is not compatible with CUDA 9.2")
   endif ()
   add_subdirectory(EB)
endif ()

if (ENABLE_LINEAR_SOLVERS)
   add_subdirectory(LinearSolvers)
   target_link_libraries(amrex
      PUBLIC
      $<BUILD_INTERFACE:Flags_Fortran_REQUIRED>
      )
endif ()

if (ENABLE_FORTRAN_INTERFACES)
   add_subdirectory(F_Interfaces)
endif ()

if (ENABLE_PARTICLES)
   add_subdirectory(Particle)
endif ()

#
# Optional external components
#
if (ENABLE_AMRDATA)
   add_subdirectory(Extern/amrdata)
   target_link_libraries(amrex
      PUBLIC
      $<BUILD_INTERFACE:Flags_Fortran_REQUIRED>
      )
endif()

if (ENABLE_PROFPARSER)
   add_subdirectory(Extern/ProfParser)
endif ()

if (ENABLE_SENSEI_INSITU)
   add_subdirectory(Extern/SENSEI)
endif ()

if (ENABLE_SUNDIALS)
   add_subdirectory(Extern/SUNDIALS4)
endif ()

if (ENABLE_CONDUIT)
   add_subdirectory(Extern/Conduit)
endif ()

if (ENABLE_HYPRE)
   add_subdirectory(Extern/HYPRE)
endif ()


#
# Print out summary -- do it here so we already linked all
# libs at this point
#  
print_amrex_configuration_summary ()

#
# Make all headers as PUBLIC HEADERS so that they get installed
#
get_target_property(AMREX_SOURCES amrex SOURCES)

set(AMREX_PUBLIC_HEADERS ${AMREX_SOURCES})
list(FILTER AMREX_PUBLIC_HEADERS INCLUDE REGEX "\\.H")
set_target_properties( amrex PROPERTIES PUBLIC_HEADER "${AMREX_PUBLIC_HEADERS}")

#
# If ENABLE_CUDA, make C++ files be compiled as CUDA sources
#
if (ENABLE_CUDA)
   set(AMREX_CUDA_SOURCES ${AMREX_SOURCES})
   list(FILTER AMREX_CUDA_SOURCES INCLUDE REGEX "\\.cpp")
   set_source_files_properties(${AMREX_CUDA_SOURCES} PROPERTIES LANGUAGE CUDA )
endif ()


#
# Install amrex  -- Export
#
include(AMReXInstallHelpers)
install_amrex(Flags_CXX Flags_Fortran Flags_FPE)
