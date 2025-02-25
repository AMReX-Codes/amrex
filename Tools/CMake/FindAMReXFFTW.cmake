#[=======================================================================[:
FindAMReXFFTW
-------

Finds the FFTW library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target:

``AMReX::FFTW``
  The FFTW library

Options/Control Variables
^^^^^^^^^^^^^^^^^^^^^^^^^

``AMReX_FFTW_SEARCH``
  FFTW search method (PKGCONFIG/CMAKE).
  Defaults to CMake config packages on Windows and to PkgConfig pc files on Linux/macOS.

``AMReX_FFTW_IGNORE_OMP``
  Ignore FFTW3 OpenMP support, even if found.

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``AMReXFFTW_FOUND``
  True if the FFTW library has been found.

This will also create an imported target, AMReX::FFTW.

#]=======================================================================]

# Helper Functions ############################################################
#
option(AMReX_FFTW_IGNORE_OMP "Ignore FFTW3's OpenMP support, even if found" OFF)
mark_as_advanced(AMReX_FFTW_IGNORE_OMP)

# Set the AMReX_FFTW_OMP=1 define on AMReX::FFTW if TRUE and print
# a message
#
function(fftw_add_define HAS_FFTW_OMP_LIB)
    if(HAS_FFTW_OMP_LIB)
        message(STATUS "FFTW: Found OpenMP support")
        target_compile_definitions(AMReX::FFTW INTERFACE AMReX_FFTW_OMP=1)
    else()
        message(STATUS "FFTW: Could NOT find OpenMP support")
    endif()
endfunction()

# Check if the found FFTW install location has an _omp library, e.g.,
# libfftw3(f)_omp.(a|so) shipped and if yes, set the AMReX_FFTW_OMP=1 define.
#
function(fftw_check_omp library_paths fftw_precision_suffix)
    find_library(HAS_FFTW_OMP_LIB${fftw_precision_suffix} fftw3${fftw_precision_suffix}_omp
        PATHS ${library_paths}
        # this is intentional, so we don't mix different FFTW installs
        # and only check what is in the location hinted by the
        # "library_paths" variable
        NO_DEFAULT_PATH
        NO_PACKAGE_ROOT_PATH
        NO_CMAKE_PATH
        NO_CMAKE_ENVIRONMENT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        NO_CMAKE_FIND_ROOT_PATH
    )
    if(HAS_FFTW_OMP_LIB${fftw_precision_suffix})
        # the .pc files here forget to link the _omp.a/so files
        # explicitly - we add those manually to avoid any trouble,
        # e.g., in static builds.
        target_link_libraries(AMReX::FFTW INTERFACE ${HAS_FFTW_OMP_LIB${fftw_precision_suffix}})
    endif()

    fftw_add_define("${HAS_FFTW_OMP_LIB${fftw_precision_suffix}}")
endfunction()


# Central FFTW3 Search ###############################################
#
# On Windows, try searching for FFTW3(f)Config.cmake files first
#   Installed .pc files wrongly and unconditionally add -lm
#   https://github.com/FFTW/fftw3/issues/236

# On Linux & macOS, note Autotools install bug:
#   https://github.com/FFTW/fftw3/issues/235
# Thus, rely on .pc files

set(AMReX_FFTW_SEARCH_VALUES PKGCONFIG CMAKE)
set(AMReX_FFTW_SEARCH_DEFAULT PKGCONFIG)
if(WIN32)
    set(AMReX_FFTW_SEARCH_DEFAULT CMAKE)
endif()
set(AMReX_FFTW_SEARCH ${AMReX_FFTW_SEARCH_DEFAULT}
        CACHE STRING "FFTW search method (PKGCONFIG/CMAKE)")
set_property(CACHE AMReX_FFTW_SEARCH PROPERTY STRINGS ${AMReX_FFTW_SEARCH_VALUES})
if(NOT AMReX_FFTW_SEARCH IN_LIST AMReX_FFTW_SEARCH_VALUES)
    message(FATAL_ERROR "AMReX_FFTW_SEARCH (${AMReX_FFTW_SEARCH}) must be one of ${AMReX_FFTW_SEARCH_VALUES}")
endif()
mark_as_advanced(AMReX_FFTW_SEARCH)

function(fftw_find_precision HFFTWp)
    if(AMReX_FFTW_SEARCH STREQUAL CMAKE)
        find_package(FFTW3${HFFTWp} CONFIG REQUIRED)
        set(AMReX_FFTW_LIBRARY_DIRS "${FFTW3${HFFTWp}_LIBRARY_DIRS}")
        message(STATUS "Found FFTW: ${FFTW3${HFFTWp}_DIR} (found version \"${FFTW3${HFFTWp}_VERSION}\")")
    else()
        find_package(PkgConfig REQUIRED QUIET)
        pkg_check_modules(fftw3${HFFTWp} REQUIRED IMPORTED_TARGET fftw3${HFFTWp})
        message(STATUS "Found FFTW: ${fftw3${HFFTWp}_PREFIX}")
        if(fftw3${HFFTWp}_LIBRARY_DIRS)
            set(AMReX_FFTW_LIBRARY_DIRS "${fftw3${HFFTWp}_LIBRARY_DIRS}")
        else()
            set(AMReX_FFTW_LIBRARY_DIRS "${fftw3${HFFTWp}_LIBDIR}")
        endif()
    endif()

    if(AMReX_FFTW_SEARCH STREQUAL CMAKE)
        target_link_libraries(AMReX::FFTW INTERFACE FFTW3::fftw3${HFFTWp})
    else()
        target_link_libraries(AMReX::FFTW INTERFACE PkgConfig::fftw3${HFFTWp})
    endif()

    if(AMReX_OMP)
        if(AMReX_FFTW_IGNORE_OMP)
            message(STATUS "FFTW: Requested to IGNORE OpenMP support")
        else()
            fftw_check_omp("${AMReX_FFTW_LIBRARY_DIRS}" "${HFFTWp}")
        endif()
    else()
        message(STATUS "FFTW: Did NOT search for OpenMP support (AMReX_OMP is not set)")
    endif()
endfunction()

if(NOT TARGET AMReX::FFTW)
    # Create imported target
    add_library(AMReX::FFTW INTERFACE IMPORTED GLOBAL)

    # floating point precision suffixes: we request float and double precision
    fftw_find_precision("")
    fftw_find_precision("f")
endif()

# Vars for CMake config
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AMReXFFTW
    HANDLE_COMPONENTS)
