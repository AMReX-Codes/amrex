#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "AMReX::amrex_3d" for configuration "Release"
set_property(TARGET AMReX::amrex_3d APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(AMReX::amrex_3d PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libamrex_3d.a"
  )

list(APPEND _cmake_import_check_targets AMReX::amrex_3d )
list(APPEND _cmake_import_check_files_for_AMReX::amrex_3d "${_IMPORT_PREFIX}/lib/libamrex_3d.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
