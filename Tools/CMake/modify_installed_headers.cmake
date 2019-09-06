#
# Add line #include<AMReX_Config.H> to all installed headers
#
message(STATUS "Modifying installed headers")
file(GLOB _installed_headers "${CMAKE_INSTALL_PREFIX}/include/*.H")
list(REMOVE_ITEM _installed_headers "${CMAKE_INSTALL_PREFIX}/include/AMReX_Config.H")

# Do not modify  headers ending in "I.H": these headers
# only contain templates and are only and always included by other headers
list(FILTER _installed_headers EXCLUDE REGEX "(I+\\.+H)$")


foreach( _header IN LISTS _installed_headers)
  file(READ ${_header} _contents)
  # Regex looks for:
  # - Anything before include guard
  # - #ifndef CPP_IDENTIFIER
  # - #define CPP_IDENTIFIER
  # - Rest of file
  if (_contents MATCHES "^(.*)#ifndef +([A-Za-z0-9_]+)\n#define +([A-Za-z0-9_]+) *\n(.+)")
    # We have found an include guard
    if (${CMAKE_MATCH_2} STREQUAL ${CMAKE_MATCH_3})
      # The guard defines match and are hence OK
    else()
      message(FATAL_ERROR "Include guard is wrong up for header file ${_header}")
    endif()

    # Insert the #include <AMReX_Config.H> between include guard and the rest of the file
    string(CONCAT "_contents" "${CMAKE_MATCH_1}" "#ifndef " "${CMAKE_MATCH_2}" "\n#define " "${CMAKE_MATCH_3}" "\n#include <AMReX_Config.H>\n" "${CMAKE_MATCH_4}")

  else()
    message(WARNING "No include guard in header ${_header}")
    string(PREPEND _contents "#include <AMReX_Config.H>\n")
  endif()
  # Rewrite
  file(WRITE ${_header} "${_contents}")
endforeach()
