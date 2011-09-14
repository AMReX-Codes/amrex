# -*- mode: cmake -*-
#
# #############################################################################
#    
# CCSE Configuration Report
#
# #############################################################################
include(FeatureSummary)


# Grab global information about enabled languages and packages
# that have returned successfully from a find_package call
get_property(lang_enabled GLOBAL PROPERTY ENABLED_LANGUAGES)
get_property(pack_found GLOBAL PROPERTY PACKAGES_FOUND)
get_property(pack_not_found GLOBAL PROPERTY PACKAGES_NOT_FOUND)

# Define the build timestamp
set(build_timestamp "Not available on this platform")
if (UNIX)
    execute_process(COMMAND "date"
                    RESULT_VARIABLE _ret_code
                    OUTPUT_VARIABLE _stdout
                    ERROR_VARIABLE  _stderr
                    )
    string(REGEX REPLACE "[\n\r]" "" build_timestamp ${_stdout})
endif()    

# Useful macros
macro(_write_to_log line)
    file(APPEND ${CCSE_CONFIG_LOG} "${line}\n")
endmacro(_write_to_log)

macro(_write_blank_line)
    _write_to_log("")
endmacro(_write_blank_line)    

macro(_create_log_file timestamp)
    file(WRITE ${CCSE_CONFIG_LOG} "${PROJECT_NAME} Configuration\n")
    _write_blank_line()
    _write_to_log("Version ${CCSE_VERSION}")
    _write_blank_line()
    _write_to_log("Build timestamp: ${timestamp}")
    _write_blank_line()
endmacro(_create_log_file)    

# Create the log file
_create_log_file(${build_timestamp})

# Write System Information
_write_blank_line()
_write_to_log( "System Information")
_write_to_log( "\tSystem:         ${CMAKE_SYSTEM}")
_write_to_log( "\tSystem Name:    ${CMAKE_SYSTEM_NAME}")
_write_to_log( "\tSystem Version: ${CMAKE_SYSTEM_VERSION}")
if ( APPLE )
    _write_to_log( "\tSystem Type:    Mac OSX")
endif(APPLE)    
if ( WIN32 )
    _write_to_log( "\tSystem Type:    Windows")
endif(WIN32)    
if ( UNIX )
    _write_to_log( "\tSystem Type:    UNIX based platform")
endif(UNIX)   

# Write Compiler Information
_write_blank_line()
_write_to_log( "Compilers")
_write_to_log( "\tEnabled Languages: ${lang_enabled}")
_write_to_log( "\tC COMPILER      ${CMAKE_C_COMPILER}")
_write_to_log( "\tC COMPILER ID   ${CMAKE_C_COMPILER_ID}")
_write_to_log( "\tCXX COMPILER    ${CMAKE_CXX_COMPILER}")
_write_to_log( "\tCXX COMPILER ID ${CMAKE_CXX_COMPILER_ID}")
if (CMAKE_FORTRAN_COMPILER_LOADED)
    _write_to_log( "\tFortran Compiler ${CMAKE_FORTRAN_COMPILER}")
    _write_to_log( "\tFortran Compiler ID ${CMAKE_FORTRAN_COMPILER_ID}")
endif()    
_write_blank_line()
_write_to_log( "Build type ${CMAKE_BUILD_TYPE}")
_write_to_log( "")
_write_to_log( "Compile Flags")
_write_to_log("\tNo information available at this time")
#_write_to_log( "\tFlags:       ${CMAKE_REQUIRED_FLAGS}")
#_write_to_log( "\tDefinitions: ${CMAKE_REQUIRED_DEFINITIONS}")
#_write_to_log( "\tIncludes:    ${CMAKE_REQUIRED_INCLUDES}")
#_write_to_log( "\tLibraries:   ${CMAKE_REQUIRED_LIBRARIES}")
_write_blank_line()

# Write Package Information
_write_to_log("Third Party Libraries")
if(pack_found)
    list(LENGTH pack_found num_pack_found)
    _write_to_log("\tNumber Of Packages Found: ${num_pack_found}")
    _write_to_log("\tPackages: ${pack_found}")
    _write_blank_line()

    foreach(pack ${pack_found})
        _write_to_log( "\t${pack}")
        _write_to_log( "\t\t${pack}_INCLUDE_DIRS=${${pack}_INCLUDE_DIRS}")
        _write_to_log( "\t\t${pack}_INCLUDE_DIR=${${pack}_INCLUDE_DIR}")
        _write_to_log( "\t\t${pack}_LIBRARY=${${pack}_LIBRARY}")
        _write_to_log( "\t\t${pack}_LIBRARIES=${${pack}_LIBRARIES}")
        _write_to_log( "\t\t${pack}_LIBRARY_DIR=${${pack}_LIBRARY_DIR}")
        _write_blank_line()
    endforeach()    
endif(pack_found)
if(pack_not_found)
    _write_blank_line()
    list(LENGTH pack_not_found num_pack_not_found)
    _write_to_log("\tNumber Of Packages Not Found: ${num_pack_not_found}") 
    if ( num_pack_not_found GREATER 0 )
        _write_to_log("\tDid not find packages: ${pack_not_found}")
    else()
        _write_to_log("\tLocated all required and requested packages")
    endif()    
else()    
    _write_blank_line()
    _write_to_log("\tLocated all required and requested packages")
endif()

