#[=======================================================================[:
FindHYPRE
-------

Finds the HYPRE library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target, if found:

``HYPRE``
  The HYPRE library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``HYPRE_FOUND``
  True if the hypre library has been found.
``HYPRE_VERSION``
  The version of the HYPRE library which was found.
``HYPRE_INCLUDE_DIRS``
  Include directories needed to use HYPRE.
``HYPRE_LIBRARIES``
  Libraries needed to link to HYPRE.
#]=======================================================================]

# Find include directories
find_path(HYPRE_INCLUDE_DIRS NAMES HYPRE.h)

# Find libraries
find_library(HYPRE_LIBRARIES NAMES HYPRE)

# Get version for config file
find_path(_config_h_path  NAMES HYPRE_config.h)

if (_config_h_path)
   file( STRINGS ${_config_h_path}/HYPRE_config.h
      _version_string REGEX "HYPRE_RELEASE_VERSION")
   string(REGEX MATCHALL "[0-9]+"  HYPRE_VERSION "${_version_string}")
   string(REPLACE ";" "." HYPRE_VERSION "${HYPRE_VERSION}")
endif ()
unset(_config_h_path CACHE)
unset(_version_string CACHE)


include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(HYPRE
   REQUIRED_VARS
   HYPRE_LIBRARIES
   HYPRE_INCLUDE_DIRS
   VERSION_VAR
   HYPRE_VERSION
   )

mark_as_advanced(HYPRE_VERSION HYPRE_LIBRARIES HYPRE_INCLUDE_DIRS)

# Create imported target
if (HYPRE_FOUND AND NOT TARGET HYPRE)
   add_library(HYPRE UNKNOWN IMPORTED)
   set_target_properties(HYPRE
      PROPERTIES
      IMPORTED_LOCATION "${HYPRE_LIBRARIES}"
      INTERFACE_INCLUDE_DIRECTORIES "${HYPRE_INCLUDE_DIRS}"
      )
endif ()  



