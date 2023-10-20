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
   set(MODULE_PATH       ${CMAKE_FILES_DIR}/AMReXCMakeModules)       # Relative path to top level installation/build-tree

   # Write Config file -- this is designed to work for both install and build trees
   configure_package_config_file(${AMREX_CMAKE_MODULES_PATH}/AMReXConfig.cmake.in
      ${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIR}/AMReXConfig.cmake
      INSTALL_DESTINATION ${CMAKE_FILES_DIR}
      PATH_VARS MODULE_PATH
      NO_CHECK_REQUIRED_COMPONENTS_MACRO)  # We have our own check_required_components

   write_basic_package_version_file(
       ${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIR}/AMReXConfigVersion.cmake
       COMPATIBILITY AnyNewerVersion )

   if(AMReX_INSTALL)
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
          DESTINATION ${CMAKE_FILES_DIR} )

       #
       # alias: last element will be legacy target
       #
       list(LENGTH AMReX_SPACEDIM list_len)
       math(EXPR list_last "${list_len} - 1")
       list(GET AMReX_SPACEDIM ${list_last} AMReX_SPACEDIM_LAST)

       # legacy symlink for: libamrex.[so|a] / amrex.[dll.lib]
       #   escape spaces for generated cmake_install.cmake file
       file(TO_CMAKE_PATH "${CMAKE_INSTALL_PREFIX}/lib" ABS_INSTALL_LIB_DIR)
       install(CODE "file(CREATE_LINK
           $<TARGET_FILE_NAME:amrex_${AMReX_SPACEDIM_LAST}d>
           \"${ABS_INSTALL_LIB_DIR}/$<TARGET_FILE_PREFIX:amrex_${AMReX_SPACEDIM_LAST}d>amrex$<TARGET_FILE_SUFFIX:amrex_${AMReX_SPACEDIM_LAST}d>\"
           COPY_ON_ERROR SYMBOLIC)"
       )

       # Install fortran modules if Fortran is enabled
       get_property(_lang GLOBAL PROPERTY ENABLED_LANGUAGES)
       if ("Fortran" IN_LIST _lang AND "amrex_${AMReX_SPACEDIM_LAST}d" IN_LIST _targets)
          get_target_property(_mod_dir "amrex_${AMReX_SPACEDIM_LAST}d" Fortran_MODULE_DIRECTORY )
          install( DIRECTORY ${_mod_dir}/ DESTINATION include ) # Trailing backslash is crucial here!
       endif ()

       # Install Tools directory
       install(
          DIRECTORY
            ${PROJECT_SOURCE_DIR}/Tools/C_scripts
            ${PROJECT_SOURCE_DIR}/Tools/typechecker
          DESTINATION
            share/amrex
          USE_SOURCE_PERMISSIONS
          )
       install(
          DIRECTORY
            ${PROJECT_SOURCE_DIR}/Tools/CMake/
          DESTINATION
            ${MODULE_PATH}
          USE_SOURCE_PERMISSIONS
          )
   endif()

   #
   # Export build-tree
   #
   export(TARGETS ${_targets} NAMESPACE AMReX::
          FILE ${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIR}/AMReXTargets.cmake)

   # Copy Tools directory to build tree
   file(
      COPY
        ${PROJECT_SOURCE_DIR}/Tools/C_scripts
        ${PROJECT_SOURCE_DIR}/Tools/typechecker
      DESTINATION
        ${PROJECT_BINARY_DIR}/share/amrex
      )
   file(
      COPY
        ${PROJECT_SOURCE_DIR}/Tools/CMake/
      DESTINATION
        ${PROJECT_BINARY_DIR}/${MODULE_PATH}
      )


endfunction ()


#
# Create a test_install target
#
# _dir        is the project directory
# _amrex_root is the amrex installation dir
#
macro( add_test_install_target _dir _amrex_root )

   get_filename_component( _dirname ${_dir} NAME )
   set(_builddir  ${CMAKE_CURRENT_BINARY_DIR}/${_dirname})
   set(_enable_fortran)
   if(AMReX_FORTRAN)
      set(_enable_fortran -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER})
   endif()
   add_custom_target(test_install
      COMMAND ${CMAKE_COMMAND} -E echo ""
      COMMAND ${CMAKE_COMMAND} -E echo "------------------------------------"
      COMMAND ${CMAKE_COMMAND} -E echo "     Testing AMReX installation     "
      COMMAND ${CMAKE_COMMAND} -E echo "------------------------------------"
      COMMAND ${CMAKE_COMMAND} -E echo ""
      COMMAND ${CMAKE_COMMAND} -E make_directory ${_builddir}
      COMMAND ${CMAKE_COMMAND} -E echo "Configuring test project"
      COMMAND ${CMAKE_COMMAND} -S ${_dir} -B ${_builddir} -DAMReX_ROOT=${_amrex_root} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} ${_enable_fortran}
      COMMAND ${CMAKE_COMMAND} -E echo "Building test project"
      COMMAND ${CMAKE_COMMAND} --build ${_builddir}
      COMMAND ${CMAKE_COMMAND} -E echo ""
      COMMAND ${CMAKE_COMMAND} -E echo "------------------------------------"
      COMMAND ${CMAKE_COMMAND} -E echo "   AMReX is installed correctly"
      COMMAND ${CMAKE_COMMAND} -E echo "              Enjoy!           "
      COMMAND ${CMAKE_COMMAND} -E echo "------------------------------------"
      COMMAND ${CMAKE_COMMAND} -E echo ""
      COMMAND ${CMAKE_COMMAND} -E remove_directory ${_builddir} # So we can run it again from scratch
      )

endmacro()
