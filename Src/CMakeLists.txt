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
   add_subdirectory(EB)
endif ()

if (ENABLE_LINEAR_SOLVERS)
   add_subdirectory(LinearSolvers)
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
endif()

if (ENABLE_PROFPARSER)
   add_subdirectory(Extern/ProfParser)
endif ()

if (ENABLE_3D_NODAL_MLMG)
   add_subdirectory(Extern/Algoim)
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
# Here we generate AMReX_BuildInfo.cpp
# If Python < 2.7, script won't work so
# we do not include this in the library
#
#
# Check for Python >= 2.7: this is needed to build AMReX_BuildInfo.cpp
# This is not a required package. This provides PYTHONINTERP_FOUND and
# PYTHON_VERSION_STRING, PYTHON_VERSION_MAJOR and PYTHON_VERSION_MINOR
#
find_package (PythonInterp QUIET)

if ( (NOT (PYTHON_VERSION_STRING VERSION_LESS "2.7") ) AND PYTHONINTERP_FOUND )
   add_custom_command(
      COMMAND  ${PROJECT_SOURCE_DIR}/Tools/C_scripts/makebuildinfo_C.py
      --amrex_home "${PROJECT_SOURCE_DIR}"
      --COMP ${CMAKE_C_COMPILER_ID} --COMP_VERSION ${CMAKE_C_COMPILER_VERSION}
      --FCOMP ${CMAKE_Fortran_COMPILER_ID} --FCOMP_VERSION ${CMAKE_C_COMPILER_VERSION}
      #--GIT ". ${AMREX_GIT}"
      OUTPUT AMReX_buildInfo.cpp
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating AMReX_buildInfo.cpp" )
  
   target_sources( amrex PRIVATE  ${CMAKE_CURRENT_BINARY_DIR}/AMReX_buildInfo.cpp )
endif ()

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
# Install amrex
# 
include(AMReXInstallHelpers)
install_amrex()


