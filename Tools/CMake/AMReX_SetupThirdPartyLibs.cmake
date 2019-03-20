
if (NOT TARGET amrex)
   message(FATAL_ERROR "Target amrex must be defined before including AMReX_SetupThirdpartyLibs.cmake")
endif ()


#
# Check for Conduit Support
#
if ( ENABLE_CONDUIT )
    if (NOT CONDUIT_DIR)
        message(FATAL_ERROR "Conduit support requires CONDUIT_DIR")
    endif ()

    if(NOT EXISTS ${CONDUIT_DIR}/lib/cmake/ConduitConfig.cmake)
        MESSAGE(FATAL_ERROR "Could not find Conduit CMake include file (${CONDUIT_DIR}/lib/cmake/ConduitConfig.cmake)")
    endif()

    #
    # Use CMake's find_package to import conduit's targets
    #
    find_package(Conduit REQUIRED QUIET
                 NO_DEFAULT_PATH
                 PATHS ${CONDUIT_DIR}/lib/cmake)

    target_link_libraries( amrex PUBLIC conduit::conduit)

    if( ENABLE_MPI )
        target_link_libraries( amrex PUBLIC conduit::conduit_mpi )
    endif ()

    MESSAGE(STATUS "Found Conduit at ${CONDUIT_DIR}" )

endif ()


#
# Check for Ascent Support
#
if ( ENABLE_ASCENT )
    if (NOT ASCENT_DIR)
        message(FATAL_ERROR "Ascent support requires ASCENT_DIR")
    endif ()

    if(NOT EXISTS ${ASCENT_DIR}/lib/cmake/AscentConfig.cmake)
        MESSAGE(FATAL_ERROR "Could not find Ascent CMake include file (${ASCENT_DIR}/lib/cmake/AscentConfig.cmake)")
    endif()

    #
    # Use CMake's find_package to import ascent's targets
    #
    find_package(Ascent REQUIRED QUIET
                 NO_DEFAULT_PATH
                 PATHS ${ASCENT_DIR}/lib/cmake)


    if( NOT ENABLE_MPI )
        target_link_libraries( amrex PUBLIC ascent::ascent )
    else()
        target_link_libraries( amrex PUBLIC ascent::ascent_mpi )
    endif ()

    MESSAGE(STATUS "Found Ascent at ${ASCENT_DIR}" )

endif ()
