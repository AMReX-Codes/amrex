#
#
#
function ( generate_amrex_config_header )

   get_target_property(_defines amrex COMPILE_DEFINITIONS)
   evaluate_genex(_defines _defines_list )

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
   endif ()

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

   # Check if the given arguments are vvalid targets
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

   configure_package_config_file(${AMREX_CMAKE_MODULES_PATH}/AMReXConfig.cmake.in
      ${PROJECT_BINARY_DIR}/export/AMReXConfig.cmake
      INSTALL_DESTINATION bin/cmake/AMReX )

   write_basic_package_version_file(
       ${PROJECT_BINARY_DIR}/export/AMReXConfigVersion.cmake
       COMPATIBILITY AnyNewerVersion )

   install( FILES
      ${PROJECT_BINARY_DIR}/export/AMReXConfig.cmake
      ${PROJECT_BINARY_DIR}/export/AMReXConfigVersion.cmake
      DESTINATION lib/cmake/AMReX )

   # Setup for target amrex installation
   install(
      TARGETS       ${_targets}
      EXPORT        AMReXTargets
      ARCHIVE       DESTINATION lib
      LIBRARY       DESTINATION lib
      INCLUDES      DESTINATION include # Adds proper directory to INTERFACE_INCLUDE_DIRECTORIES
      PUBLIC_HEADER DESTINATION include
      )

   # Setup for export AMReXtargets installation
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
   install(DIRECTORY ${PROJECT_SOURCE_DIR}/Tools/CMake
      ${PROJECT_SOURCE_DIR}/Tools/C_scripts
      DESTINATION Tools
      USE_SOURCE_PERMISSIONS
      )

   # Modify installed headers by calling external script: add #include<AMReX_Config.H>
   install(SCRIPT "${AMREX_CMAKE_MODULES_PATH}/modify_installed_headers.cmake" )

endfunction ()
