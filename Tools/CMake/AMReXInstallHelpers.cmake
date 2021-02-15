#
# Manage AMReX installation process
# Extra arguments are the targets to install
#
function (install_amrex_targets)

   # Check if the given arguments are valid targets
   set(_targets)
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
      PATH_VARS MODULE_PATH
      NO_CHECK_REQUIRED_COMPONENTS_MACRO)  # We have our own check_required_components

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
      RUNTIME       DESTINATION bin
      )

   install( EXPORT AMReXTargets
      NAMESPACE AMReX::
      DESTINATION lib/cmake/AMReX )

   # Install fortran modules if Fortran is enabled
   get_property(_lang GLOBAL PROPERTY ENABLED_LANGUAGES)
   if ("Fortran" IN_LIST _lang AND "amrex" IN_LIST _targets)
      get_target_property(_mod_dir amrex Fortran_MODULE_DIRECTORY )
      install( DIRECTORY ${_mod_dir}/ DESTINATION include ) # Trailing backslash is crucial here!
   endif ()

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
