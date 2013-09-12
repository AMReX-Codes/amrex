# -*- mode: cmake -*-
#
# Amanzi CCSE Find Module
#
# Usage:
#    Control the search through CCSE_DIR or setting environment variable
#    CCSE_ROOT to the ccse installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    CCSE_FOUND            (BOOL)       Flag indicating if CCSE was found
#    CCSE_INCLUDE_DIR      (PATH)       Path to the CCSE include file
#    CCSE_INCLUDE_DIRS     (LIST)       List of all required include files
#    CCSE_LIBRARY_DIR      (PATH)       Path to the CCSE library
#    CCSE_LIBRARIES        (LIST)       List of all required CCSE libraries
#    CCSE_PERL_DIR         (PATH)       Path to CCSE Perl scripts
#    CCSE_EXT_LIBRARIES    (LIST)
#    CCSE_EXT_LIBRARY_DIRS (LIST)
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

if (ENABLE_PETSC)
    find_package(X11 QUIET REQUIRED)
    if (NOT X11_FOUND)
        message(WARNING "If ENABLE_PETSC, need X11 but cmake failed to find"
                        " Build will likely fail.")
    endif()
endif()

include(CCSEOptions)


if ( CCSE_LIBRARIES AND CCSE_INCLUDE_DIRS AND CCSE_PERL_DIR)

    # Do nothing. Variables are set. No need to search again

else(CCSE_LIBRARIES AND CCSE_INCLUDE_DIRS AND CCSE_PERL_DIR)

    # Cache variables
    if(CCSE_DIR)
        set(CCSE_DIR "${CCSE_DIR}" CACHE PATH "Path to search for CCSE include and library files")
    endif()

    if(CCSE_INCLUDE_DIR)
        set(CCSE_INCLUDE_DIR "${CCSE_INCLUDE_DIR}" CACHE PATH "Path to search for CCSE include files")
    endif()

    if(CCSE_LIBRARY_DIR)
        set(CCSE_LIBRARY_DIR "${CCSE_LIBRARY_DIR}" CACHE PATH "Path to search for CCSE library files")
    endif()
    
    if(CCSE_PERL_DIR)
        set(CCSE_PERL_DIR "${CCSE_PERL_DIR}" CACHE PATH "Path to search for CCSE perl scripts")
    endif()

    
    # Search for include files
    # Search order preference:
    #  (1) CCSE_INCLUDE_DIR - check existence of path AND if the include files exist
    #  (2) CCSE_DIR/include
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(ccse_inc_names "BoxLib.H")

    if (CCSE_INCLUDE_DIR)

        if (EXISTS "${CCSE_INCLUDE_DIR}")

            find_path(ccse_test_include_path
                      NAMES ${ccse_inc_names}
                      HINTS ${CCSE_INCLUDE_DIR}
                      NO_DEFAULT_PATH)
            if(NOT ccse_test_include_path)
                message(SEND_ERROR "Cannot locate ${ccse_inc_names} in ${CCSE_INCLUDE_DIR}")
            endif()
            set(CCSE_INCLUDE_DIR ${ccse_test_include_path} )
        else()
            message(SEND_ERROR "CCSE_INCLUDE_DIR=${CCSE_INCLUDE_DIR} does not exist")
            set(CCSE_INCLUDE_DIR "CCSE_INCLUDE_DIR-NOTFOUND")
        endif()

    else() 

        set(ccse_inc_suffixes "include")
        if(CCSE_DIR)

            if (EXISTS "${CCSE_DIR}" )

                find_path(CCSE_INCLUDE_DIR
                          NAMES ${ccse_inc_names}
                          HINTS ${CCSE_DIR}
                          PATH_SUFFIXES ${ccse_inc_suffixes}
                          NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "CCSE_DIR=${CCSE_DIR} does not exist")
                 set(CCSE_INCLUDE_DIR "CCSE_INCLUDE_DIR-NOTFOUND")
            endif()    


        else()

            find_path(CCSE_INCLUDE_DIR
                      NAMES ${ccse_inc_names}
                      PATH_SUFFIXES ${ccse_inc_suffixes})

        endif()

    endif()

    if ( NOT CCSE_INCLUDE_DIR )
        message(SEND_ERROR "Cannot locate CCSE include directory")
    else()
        set(CCSE_INCLUDE_DIRS ${CCSE_INCLUDE_DIR})
    endif()

    # Search for libraries 
    # Search order preference:
    #  (1) CCSE_LIBRARY_DIR - check existence of path AND if the include files exist
    #  (2) CCSE_DIR/<lib,Lib>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(ccse_lib_name cboxlib)
    if (CCSE_LIBRARY_DIR)

        if (EXISTS "${CCSE_LIBRARY_DIR}")

            find_library(CCSE_LIBRARY
                         NAMES ${ccse_lib_name}
                         HINTS ${CCSE_LIBRARY_DIR}
                         NO_DEFAULT_PATH)
        else()
            message(SEND_ERROR "CCSE_LIBRARY_DIR=${CCSE_LIBRARY_DIR} does not exist")
            set(CCSE_LIBRARY "CCSE_LIBRARY-NOTFOUND")
        endif()

    else() 

        if(CCSE_DIR)

            if (EXISTS "${CCSE_DIR}" )

                find_library(CCSE_LIBRARY
                             NAMES ${ccse_lib_name}
                             HINTS ${CCSE_DIR}
                             PATH_SUFFIXES "lib"
                             NO_DEFAULT_PATH)
                get_filename_component(CCSE_LIBRARY_DIR "${CCSE_LIBRARY}" PATH)

            else()
                 message(SEND_ERROR "CCSE libs not found in CCSE_DIR/lib${CCSE_LIBDIR_MPI_SUFFIX}${CCSE_LIBDIR_OMP_SUFFIX}, (CCSE_DIR=${CCSE_DIR})")
                 set(CCSE_LIBRARY "CCSE_LIBRARY-NOTFOUND")
            endif()    

        else()

            find_library(CCSE_LIBRARY
                         NAMES ${ccse_lib_name})

            get_filename_component(CCSE_LIBRARY_DIR "${CCSE_LIBRARY}" PATH)
            
        endif()

    endif()

    # Now, make sure the rest are in the same place
    set(CCSE_LIBRARIES cboxlib;fboxlib;cfboxlib;box_camrdata;gslib)

    foreach (L ${CCSE_LIBRARIES})

            find_library(CCSE_LIBRARY
                         NAMES ${L}
                         HINTS ${CCSE_LIBRARY_DIR}
                         NO_DEFAULT_PATH)

            if ( NOT CCSE_LIBRARY )
                message(SEND_ERROR "Cannot locate CCSE library: ${L}")
            endif()

    endforeach()
    set(CCSE_LIBRARY_DIRS ${CCSE_LIBRARY_DIR})

    # Search for perl scripts
    # Search order preference:
    #  (1) CCSE_PERL_DIR - check existence of path AND if the perl script exist
    #  (2) CCSE_DIR/perl
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(ccse_perl_name strip72)

    if (CCSE_PERL_DIR)

        if (EXISTS "${CCSE_PERL_DIR}")

            find_path(CCSE_PERL
                      NAMES ${ccse_perl_name}
                      HINTS ${CCSE_PERL_DIR}
                      NO_DEFAULT_PATH)
        else()
            message(SEND_ERROR "CCSE_PERL_DIR=${CCSE_PERL_DIR} does not exist")
            set(CCSE_PERL "CCSE_PERL-NOTFOUND")
        endif()

    else() 

        if(CCSE_DIR)

            if (EXISTS "${CCSE_DIR}" )

                find_path(CCSE_PERL
                          NAMES ${ccse_perl_name}
                          HINTS "${CCSE_DIR}/perl"
                          NO_DEFAULT_PATH)

                if ( NOT CCSE_PERL )
                    message(SEND_ERROR "CCSE Perl scripts not in CCSE_DIR/perl=${CCSE_DIR}/perl")
                else()
                    set(CCSE_PERL_DIR ${CCSE_PERL})
                endif()    

            else()

                message(SEND_ERROR "CCSE_DIR/perl=${CCSE_DIR}/perl does not exist")

            endif()

        endif()

    endif()

    if (ENABLE_PETSC)

        set(PETSC_DIR $ENV{PETSC_DIR})
        if ("${PETSC_DIR}" STREQUAL "")
            message(SEND_ERROR "Must define env variable PETSC_DIR if ENABLE_PETSC=ON")
        endif()

        message(STATUS "CCSE requires PETSc and X11 since ENABLE_PETSC=ON")
        message(STATUS "     using PETSC_DIR=${PETSC_DIR}")

        # NOTE: Not sure why we have to explicitly include X11 stuff, since FindX11 was supposed to do it...
        set(CCSE_EXT_LIBRARIES petsc ${X11_LIBRARIES})
        list(APPEND CCSE_INCLUDE_DIRS ${PETSC_DIR}/include ${X11_INCLUDE_DIR})

    else()

        set(CCSE_EXT_LIBRARIES "")

    endif()

endif(CCSE_LIBRARIES AND CCSE_INCLUDE_DIRS AND CCSE_PERL_DIR)    

# Send useful message if everything is found
find_package_handle_standard_args(CCSE DEFAULT_MSG
                                  CCSE_LIBRARIES
                                  CCSE_LIBRARY_DIRS
                                  CCSE_INCLUDE_DIRS
                                  CCSE_PERL_DIR)

mark_as_advanced(
  CCSE_INCLUDE_DIR
  CCSE_INCLUDE_DIRS
  CCSE_LIBRARIES
  CCSE_LIBRARY_DIR
  CCSE_LIBRARY_DIRS
  CCSE_PERL_DIR
)
