#[=======================================================================[:
FindSUNDIALS
-------

Finds the SUNDIALS libraries.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``SUNDIALS::<comp>``
  The <comp> library

where <comp> is the selected component.
If no component is specified, this module will search for ``nvecserial`` only.

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``SUNDIALS_FOUND``
  True if all the selected components have been found.
``SUNDIALS_VERSION``
  The version of the SUNDIALS library which was found.
``SUNDIALS_<COMP>_INCLUDE_DIRS``
  Include directories needed to use component <COMP>.
``SUNDIALS_<COMP>_LIBRARIES``
  Libraries needed to link to component <COMP>.
#]=======================================================================]
include(FindPackageHandleStandardArgs)

# Get version for config file
find_path(_config_h_path
   NAMES          sundials_config.h
   PATH_SUFFIXES  sundials )

if (_config_h_path)
   file(STRINGS ${_config_h_path}/sundials_config.h
      _version_string REGEX "SUNDIALS_VERSION ") # space at the end of SUNDIAL_VERSION is needed!
   string(REGEX MATCHALL "[0-9]+" SUNDIALS_VERSION "${_version_string}")
   string(REPLACE ";" "." SUNDIALS_VERSION "${SUNDIALS_VERSION}")
endif ()
unset(_version_file_path CACHE)
unset(_version_string CACHE)


# Include directory is only the top level one, i.e. "include"
# Find path by using directory "sundials" as reference
find_path(_sundials_include_path NAMES sundials PATH_SUFFIXES include)


# Valid components
set(_valid_components
   nvecserial
   nvecparallel
   nvecopenmp
   nvecopenmpdev
   nveccuda
   cvode
   arkode
   )


# Search for a default library (nvecserial) or for all the
# required components
if (NOT SUNDIALS_FIND_COMPONENTS)
   set(_sundials_findlist nvecserial )
else ()
   set(_sundials_findlist ${SUNDIALS_FIND_COMPONENTS})
endif ()

list(REMOVE_DUPLICATES _sundials_findlist)

#
# Setup one imported target per component
#
set(_SUNDIALS_REQUIRED_VARS)
foreach(_comp IN LISTS _sundials_findlist)
   
   string( TOLOWER "${_comp}" _comp_lower  )
   string( TOUPPER "${_comp}" _comp_upper  )

   # Check whether component is in valid list
   if (NOT (${_comp_lower} IN_LIST _valid_components))
      message(FATAL_ERROR "Invalid component ${_comp_lower}")
   endif ()

   # Include path is always the path to the top-level "include" directory in the install tree
   set(SUNDIALS_${_comp_upper}_INCLUDE_DIRS ${_sundials_include_path})
   find_library(SUNDIALS_${_comp_upper}_LIBRARIES NAMES sundials_${_comp_lower} PATH_SUFFIXES lib lib64)

   find_package_handle_standard_args(SUNDIALS_${_comp_upper}
      REQUIRED_VARS
      SUNDIALS_${_comp_upper}_LIBRARIES
      SUNDIALS_${_comp_upper}_INCLUDE_DIRS
      VERSION_VAR SUNDIALS_VERSION
      )

   mark_as_advanced(SUNDIALS_${_comp_upper}_LIBRARIES SUNDIALS_${_comp_upper}_INCLUDE_DIRS)

   list(APPEND _SUNDIALS_REQUIRED_VARS "SUNDIALS_${_comp_upper}_FOUND")
   
   # Create imported target
   set(_target SUNDIALS::${_comp_lower})
   if (SUNDIALS_${_comp_upper}_FOUND AND NOT TARGET ${_target})
      add_library(${_target} UNKNOWN IMPORTED)
      set_target_properties(${_target} PROPERTIES
         IMPORTED_LOCATION "${SUNDIALS_${_comp_upper}_LIBRARIES}"
         INTERFACE_INCLUDE_DIRECTORIES "${SUNDIALS_${_comp_upper}_INCLUDE_DIRS}"
         )
   endif ()  

endforeach()

unset(_sundials_include_path  CACHE)
unset(_valid_components       CACHE)
unset(_sundials_findlist      CACHE)


# Declare full package found
find_package_handle_standard_args(SUNDIALS
   REQUIRED_VARS ${_SUNDIALS_REQUIRED_VARS}
   VERSION_VAR   SUNDIALS_VERSION
   )

unset(_SUNDIALS_REQUIRED_VARS CACHE)
