#
#
#
function ( generate_amrex_config_header )

   get_target_property(_defines_list amrex COMPILE_DEFINITIONS)
   eval_genex(_defines_list NONE NONE INTERFACE BUILD)

   # set variables from list of defines
   # configure_file() will use these variables to
   # define or undefine items in AMReX_Config.H
   foreach(_define IN LISTS _defines_list)
      string(FIND "${_define}" "=" idx)
      if (idx EQUAL -1) # define without value
         set(${_define} " ")
      else ()   # define with value
         string(SUBSTRING "${_define}" 0      ${idx}  _def)
         string(SUBSTRING "${_define}" ${idx} -1      _val)
         string(REPLACE "=" " " _val "${_val}")
         string(STRIP "${_val}" _val)
         set(${_def} ${_val})
      endif ()

   endforeach ()

   set(COMP_DECLS)
   if (NOT ALLOW_DIFFERENT_COMPILER)
       if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
           set(COMPILER_ID_MACRO __GNUC__)
       elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
           set(COMPILER_ID_MACRO __INTEL_COMPILER)
       elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Cray" )
           set(COMPILER_ID_MACRO  __CRAYC)
       elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI" )
           set(COMPILER_ID_MACRO  __PGI)
       elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" )
           set(COMPILER_ID_MACRO  __llvm__)
       elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" )
           set(COMPILER_ID_MACRO  __llvm__)
       elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC" )
           set(COMPILER_ID_MACRO  _MSC_VER)
       else ()
           message(FATAL_ERROR "Compiler '${CMAKE_CXX_COMPILER_ID}' not supported by AMReX developers! "
               "Try to configure with -DALLOW_DIFFERENT_COMPILER=ON")
       endif ()
       set(msg "libamrex was built with ${CMAKE_CXX_COMPILER_ID}. To avoid this error, reconfigure with -DALLOW_DIFFERENT_COMPILER=ON")
       set(COMP_DECLS "\n#ifndef ${COMPILER_ID_MACRO}\nstatic_assert(false,\"${msg}\");\n#endif")
   endif()

   if (ENABLE_OMP)
      set(OMP_DECLS "#ifndef _OPENMP\nstatic_assert(false,\"libamrex was built with OpenMP\");\n#endif")
   else ()
      set(OMP_DECLS "#ifdef _OPENMP\nstatic_assert(false,\"libamrex was built without OpenMP\");\n#endif")
   endif ()

   configure_file(${AMREX_CMAKE_MODULES_PATH}/AMReX_Config.H.in
      "${CMAKE_CURRENT_BINARY_DIR}/AMReX_Config.H")

   install(FILES "${CMAKE_CURRENT_BINARY_DIR}/AMReX_Config.H" DESTINATION include)

endfunction ()


#
# Manage AMReX installation process
# Extra arguments are the targets to install
#
function (install_amrex)

   # Check if the given arguments are valid targets
   set(_targets amrex)
   foreach (_arg IN LISTS ARGV)
      if (TARGET ${_arg})
         set(_targets ${_targets} ${_arg})
      else ()
         message(WARNING "Target ${_arg} does not exist: skipping installation")
      endif ()
   endforeach ()

   # Write and install configure file
   include(CMakePackageConfigHelpers)

   #
   # Setup for install and export of build tree
   #
   # Relative path to top level installation/build-tree
   if(WIN32)
       set(CMAKE_FILES_DIR   cmake)   
   else()
       set(CMAKE_FILES_DIR   lib/cmake/AMReX)
   endif()
   set(MODULE_PATH       Tools/CMake)       # Relative path to top level installation/build-tree

   # Write Config file -- this is designed to work for both install and build trees
   configure_package_config_file(${AMREX_CMAKE_MODULES_PATH}/AMReXConfig.cmake.in
      ${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIR}/AMReXConfig.cmake
      INSTALL_DESTINATION ${CMAKE_FILES_DIR}
      PATH_VARS MODULE_PATH)

   write_basic_package_version_file(
       ${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIR}/AMReXConfigVersion.cmake
       COMPATIBILITY AnyNewerVersion )

   install( FILES
      ${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIR}/AMReXConfig.cmake
      ${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIR}/AMReXConfigVersion.cmake
      DESTINATION ${CMAKE_FILES_DIR} )

   #
   # Export install-tree
   #
   install(
      TARGETS       ${_targets}
      EXPORT        AMReXTargets
      ARCHIVE       DESTINATION lib
      LIBRARY       DESTINATION lib
      INCLUDES      DESTINATION include # Adds proper directory to INTERFACE_INCLUDE_DIRECTORIES
      PUBLIC_HEADER DESTINATION include
      )

   install( EXPORT AMReXTargets
      NAMESPACE AMReX::
      DESTINATION lib/cmake/AMReX )

   # Install fortran modules if Fortran is enabled
   get_property(_lang GLOBAL PROPERTY ENABLED_LANGUAGES)
   if ("Fortran" IN_LIST _lang )
      get_target_property(_mod_dir amrex Fortran_MODULE_DIRECTORY )
      install( DIRECTORY ${_mod_dir}/ DESTINATION include ) # Trailing backslash is crucial here!
   endif ()

   # Generate config header
   generate_amrex_config_header()

   # Install Tools directory
   install(
      DIRECTORY
      ${PROJECT_SOURCE_DIR}/Tools/CMake
      ${PROJECT_SOURCE_DIR}/Tools/C_scripts
      ${PROJECT_SOURCE_DIR}/Tools/typechecker
      DESTINATION
      Tools
      USE_SOURCE_PERMISSIONS
      )

   # Modify installed headers by calling external script: add #include<AMReX_Config.H>
   install(SCRIPT "${AMREX_CMAKE_MODULES_PATH}/modify_installed_headers.cmake" )

   #
   # Export build-tree
   #
   export( EXPORT AMReXTargets NAMESPACE AMReX::
      FILE ${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIR}/AMReXTargets.cmake )

   # Copy Tools directory to build tree
   file(
      COPY
      ${PROJECT_SOURCE_DIR}/Tools/CMake
      ${PROJECT_SOURCE_DIR}/Tools/C_scripts
      ${PROJECT_SOURCE_DIR}/Tools/typechecker
      DESTINATION
      ${PROJECT_BINARY_DIR}/Tools
      )


endfunction ()
