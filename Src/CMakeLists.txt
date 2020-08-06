#
# Fist, define the library we want to add
# Think of this as a constructor for object "amrex"
# In the following, we will use setters to change
# the properties of this object, like setting the sources,
# setting the compile definitions and so on
#
add_library( amrex )
add_library( AMReX::amrex ALIAS amrex )

#
# Setup parallel backends (any type: inter/intra-node)
#
include(AMReXParallelBackends)

#
# Find and link third party libraries if needed
#
include(AMReXThirdPartyLibraries)

# Where to store Fortran modules
if (CMAKE_Fortran_COMPILER_LOADED)
   set_target_properties( amrex
      PROPERTIES
      Fortran_MODULE_DIRECTORY
      ${PROJECT_BINARY_DIR}/mod_files
      INTERFACE_INCLUDE_DIRECTORIES
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/mod_files>
      )
   # adding those is a work-around that avoids that downstream C++ projects have
   # to `enable_language(Fortran)`
   target_link_libraries(amrex PUBLIC ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
endif ()

# Load Flags targets and use them if no user-defined flags is given
include(AMReXFlagsTargets)

if (CMAKE_Fortran_COMPILER_LOADED AND (NOT CMAKE_Fortran_FLAGS) )
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
    #
    # Generate BLProfParser.tab.H, BLProfParser.tab.cpp and
    # BLProfParser.lex.yy.cpp
    # Must be done here where target "amrex" is defined because
    # add_custom_command create rules which are seen only within the scope
    # from which add_custom_command is called
    #
    set(_dir ${CMAKE_CURRENT_LIST_DIR}/Extern/ProfParser)
    set(_odir ${CMAKE_CURRENT_BINARY_DIR})
    add_custom_command(
        OUTPUT  ${_dir}/BLProfParser.tab.cpp ${_dir}/BLProfParser.tab.H
        COMMAND cat   BLProfParser.y $(SED0) $(SED1) > BLProfParserNC.y
        COMMAND bison --defines=BLProfParser.tab.H
        --output=BLProfParser.tab.cpp BLProfParserNC.y
        COMMAND rm    BLProfParserNC.y
        DEPENDS ${_dir}/BLProfParser.y
        WORKING_DIRECTORY ${_dir}
        COMMENT "Generating BLProfParser.tab.H and BLProfParser.tab.cpp" )

    add_custom_command(
        OUTPUT ${_dir}/BLProfParser.lex.yy.cpp
        COMMAND flex --outfile=${_dir}/BLProfParser.lex.yy.cpp BLProfParser.l
        DEPENDS ${_dir}/BLProfParser.l
        WORKING_DIRECTORY ${_dir}
        COMMENT "Generating BLProfParser.lex.yy.cpp" )

    target_sources( amrex
        PRIVATE
        ${_dir}/BLProfParser.tab.cpp
        ${_dir}/BLProfParser.lex.yy.cpp
        ${_dir}/BLProfParser.tab.H
        )
    unset(_dir)

    add_subdirectory(Extern/ProfParser)
endif ()

if (ENABLE_SENSEI_INSITU)
   add_subdirectory(Extern/SENSEI)
endif ()

if (ENABLE_SUNDIALS)
   add_subdirectory(Extern/SUNDIALS)
endif ()

if (ENABLE_CONDUIT OR ENABLE_ASCENT)
   add_subdirectory(Extern/Conduit)
endif ()

if (ENABLE_HYPRE OR ENABLE_PETSC)
   add_subdirectory(Extern/HYPRE)
endif ()

if (ENABLE_PETSC)
   add_subdirectory(Extern/PETSc)
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
# If ENABLE_CUDA, C++ files will be compiled as CUDA sources
#
if (ENABLE_CUDA)
   setup_target_for_cuda_compilation( amrex )
endif ()


#
# Install amrex  -- Export
#
include(AMReXInstallHelpers)
set(_amrex_targets Flags_CXX Flags_Fortran Flags_FPE )
if (ENABLE_DPCPP)
   list(APPEND _amrex_targets SYCL)
endif()
install_amrex(${_amrex_targets})
