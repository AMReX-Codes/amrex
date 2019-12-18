#[=======================================================================[:
FindPETSc
-------

Finds the PETSc library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target, if found:

``PETSc``
  The PETSc library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``PETSC_FOUND``
  True if the hypre library has been found.
``PETSC_VERSION``
  The version of the PETSc library which was found.
``PETSC_INCLUDE_DIRS``
  Include directories needed to use PETSc.
``PETSC_LIBRARIES``
  Libraries needed to link to PETSc.
#]=======================================================================]

# Find include directories
message(STATUS "kae")
find_path(PETSC_INCLUDE_DIRS_BASE PATHS ${PETSC_DIR}/include NAMES petsc.h )
message(STATUS "base ${PETSC_INCLUDE_DIRS_BASE}")
find_path(PETSC_INCLUDE_DIRS_ARCH PATHS ${PETSC_DIR}/${PETSC_ARCH}/include NAMES petscconf.h)
message(STATUS "arch ${PETSC_INCLUDE_DIRS_ARCH}")
set(PETSC_INCLUDE_DIRS ${PETSC_INCLUDE_DIRS_BASE} ${PETSC_INCLUDE_DIRS_ARCH})
message(STATUS "joe ${PETSC_INCLUDE_DIRS}")

#Set Fortran include directories TODO: currently overwrites any given values
set(CMAKE_Fortran_FLAGS ${PETSC_INCLUDE_DIRS})

# Find libraries
find_library(PETSC_LIBRARIES PATHS ${PETSC_DIR}/${PETSC_ARCH}/lib NAMES petsc)

# Get version for config file
find_path(_config_h_path  PATHS ${PETSC_DIR}/include NAMES petscversion.h)

if (_config_h_path)
   file( STRINGS ${_config_h_path}/petscversion.h _version_string REGEX "PETSC_VERSION_MAJOR")
   string(REGEX MATCHALL "[0-9]+"  PETSC_VERSION_MAJOR "${_version_string}")
   string(REPLACE ";" "." PETSC_VERSION_MAJOR "${PETSC_VERSION_MAJOR}")
   file( STRINGS ${_config_h_path}/petscversion.h _version_string REGEX "PETSC_VERSION_MINOR")
   string(REGEX MATCHALL "[0-9]+"  PETSC_VERSION_MINOR "${_version_string}")
   string(REPLACE ";" "." PETSC_VERSION_MINOR "${PETSC_VERSION_MINOR}")
   file( STRINGS ${_config_h_path}/petscversion.h _version_string REGEX "PETSC_VERSION_SUBMINOR")
   string(REGEX MATCHALL "[0-9]+"  PETSC_VERSION_SUBMINOR "${_version_string}")
   string(REPLACE ";" "." PETSC_VERSION_SUBMINOR "${PETSC_VERSION_SUBMINOR}")
   string(JOIN . PETSC_VERSION ${PETSC_VERSION_MAJOR} ${PETSC_VERSION_MINOR} ${PETSC_VERSION_SUBMINOR})
endif ()
unset(_config_h_path CACHE)
unset(_version_string CACHE)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(PETSc
   REQUIRED_VARS
   PETSC_LIBRARIES
   PETSC_INCLUDE_DIRS
   VERSION_VAR
   PETSC_VERSION
   )

mark_as_advanced(PETSC_VERSION PETSC_LIBRARIES PETSC_INCLUDE_DIRS)

# Create imported target
if (PETSC_FOUND AND NOT TARGET PETSC)
   add_library(PETSC UNKNOWN IMPORTED GLOBAL)
   set_target_properties(PETSC
      PROPERTIES
      IMPORTED_LOCATION "${PETSC_LIBRARIES}"
      INTERFACE_INCLUDE_DIRECTORIES "${PETSC_INCLUDE_DIRS}"
      )
endif ()  



