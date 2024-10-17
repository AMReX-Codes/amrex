#[=======================================================================[:
FindAMReXFFTW
-------

Finds the FFTW library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target, if found:

``FFTW``
  The FFTW library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``AMReXFFTW_FOUND``
  True if the hypre library has been found.
``FFTW_INCLUDES``
  Include directories needed to use FFTW.
``FFTW_LIBRARIES``
  Libraries needed to link to FFTW.

This will also create an imported target, AMReX::FFTW.
#]=======================================================================]

if (NOT FFTW_INCLUDES)
    find_path(FFTW_INCLUDES NAMES "fftw3.h" HINTS ${FFTW_ROOT}/include)
endif()

if (NOT FFTW_LIBRARIES)
    find_library(FFTW_LIBRARY NAMES "fftw3" HINTS ${FFTW_ROOT}/lib)
    find_library(FFTWF_LIBRARY NAMES "fftw3f" HINTS ${FFTW_ROOT}/lib)
    set(FFTW_LIBRARIES ${FFTW_LIBRARY} ${FFTWF_LIBRARY})
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(AMReXFFTW
    REQUIRED_VARS FFTW_LIBRARIES FFTW_INCLUDES)

mark_as_advanced(FFTW_LIBRARIES FFTW_INCLUDES)

# Create imported target
add_library(AMReX::FFTW INTERFACE IMPORTED GLOBAL)
target_link_libraries(AMReX::FFTW INTERFACE ${FFTW_LIBRARIES})
set_target_properties(AMReX::FFTW PROPERTIES
	INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDES}")
