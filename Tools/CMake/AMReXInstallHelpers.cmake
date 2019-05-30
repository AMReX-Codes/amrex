#
#
# 
function ( generate_amrex_config_header )

   get_target_property(_defines amrex COMPILE_DEFINITIONS)
   evaluate_genex(_defines _defines_list )

   set(_config_fname "${CMAKE_CURRENT_BINARY_DIR}/AMReX_Config.H")

   file(WRITE  "${_config_fname}" "#ifndef AMREX_CONFIG_H_\n")
   file(APPEND "${_config_fname}" "#define AMREX_CONFIG_H_\n" )

   foreach(_define IN LISTS _defines_list)
      string(REPLACE "=" " " _define "${_define}")
      file(APPEND "${_config_fname}" "#define ${_define}\n" )      
   endforeach ()

   file(APPEND "${_config_fname}" "#ifdef __cplusplus\n")

   if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
      file(APPEND "${_config_fname}" "#ifndef __GNUC__\n")
      file(APPEND "${_config_fname}" "static_assert(false,\"libamrex was built with GNU\");\n")
      file(APPEND "${_config_fname}" "#endif\n" )
   elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
      file(APPEND "${_config_fname}" "#ifndef __INTEL_COMPILER\n")
      file(APPEND "${_config_fname}" "static_assert(false,\"libamrex was built with Intel\");\n")
      file(APPEND "${_config_fname}" "#endif\n" )
   elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Cray" )
      file(APPEND "${_config_fname}" "#ifndef __CRAYC\n")
      file(APPEND "${_config_fname}" "static_assert(false,\"libamrex was built with Cray\");\n")
      file(APPEND "${_config_fname}" "#endif\n" )
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI" )
      file(APPEND "${_config_fname}" "#ifndef __PGI\n")
      file(APPEND "${_config_fname}" "static_assert(false,\"libamrex was built with PGI\");\n")
      file(APPEND "${_config_fname}" "#endif\n" )
   elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" )
      file(APPEND "${_config_fname}" "#ifndef __llvm__\n")
      file(APPEND "${_config_fname}" "static_assert(false,\"libamrex was built with Clang/LLVM\");\n")
      file(APPEND "${_config_fname}" "#endif\n" )   
   endif ()

   if (ENABLE_OMP)
      file(APPEND "${_config_fname}" "#ifndef _OPENMP\n")
      file(APPEND "${_config_fname}" "static_assert(false,\"libamrex was built with OpenMP\");\n")
      file(APPEND "${_config_fname}" "#endif\n" )      
   else ()
      file(APPEND "${_config_fname}" "#ifdef _OPENMP\n")
      file(APPEND "${_config_fname}" "static_assert(false,\"libamrex was built without OpenMP\");\n")
      file(APPEND "${_config_fname}" "#endif\n" )        
   endif ()

   file(APPEND "${_config_fname}" "#endif\n" )   
   file(APPEND "${_config_fname}" "#endif\n" )

   install(FILES ${_config_fname} DESTINATION include)
   
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


   # Install Version file
   # Package version is a modified form of AMREX_GIT_VERSION
   if (AMREX_GIT_VERSION)
      string(FIND "${AMREX_GIT_VERSION}" "-" _idx REVERSE)
      string(SUBSTRING "${AMREX_GIT_VERSION}" 0 "${_idx}" _pkg_version )
      string(FIND "${_pkg_version}" "-" _idx REVERSE)
      string(SUBSTRING "${_pkg_version}" 0 "${_idx}" _pkg_version )
      string(REPLACE "-" "." _pkg_version "${_pkg_version}")
   endif ()

   write_basic_package_version_file( ${PROJECT_BINARY_DIR}/export/AMReXConfigVersion.cmake
      VERSION ${_pkg_version}
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

   # Install fortran modules
   get_target_property(_mod_dir amrex Fortran_MODULE_DIRECTORY )
   install( DIRECTORY ${_mod_dir}/ DESTINATION include ) # Trailing backslash is crucial here!

   # # This header in a weird path has to be copied to install includes
   # install( FILES ${PROJECT_SOURCE_DIR}/Tools/C_scripts/AMReX_buildInfo.H
   #    DESTINATION include )

   # Generate config header
   generate_amrex_config_header()
   
   # Install Tools directory
   install(DIRECTORY ${PROJECT_SOURCE_DIR}/Tools/  DESTINATION Tools 
      USE_SOURCE_PERMISSIONS )

   # Modify installed headers by calling external script: add #include<AMReX_Config.H>
   install(SCRIPT "${AMREX_CMAKE_MODULES_PATH}/modify_installed_headers.cmake" )

endfunction ()
