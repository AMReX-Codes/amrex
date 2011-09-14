# -*- mode: cmake -*-
#
# Functions for building a link-line for linking to CCSE.
#
# See the modified CMakeLists.txt files for usage. 
#
# This functionality should be merged with InstallManager.

include(ParseLibraryList)

# From CCSEConfigReport.cmake:
set(build_timestamp "Not available on this platform")
if (UNIX)
  execute_process(COMMAND "date"
    RESULT_VARIABLE _ret_code
    OUTPUT_VARIABLE _stdout
    ERROR_VARIABLE  _stderr
    )
  string(REGEX REPLACE "[\n\r]" "" build_timestamp ${_stdout})
endif()

function(_parse_add_libraries library_list libraries_to_add)

if (library_list)
  parse_library_list(
    ${library_list}
    FOUND   libraries_split
    DEBUG   debug_libraries
    OPT     opt_libraries
    GENERAL general_libraries)

  if (libraries_split)
    message("Libraries for ${package} are present in multiple debug and/or opt versions")
    if(${CMAKE_BUILD_TYPE} MATCHES "debug")
      message("Adding debug libraries")
      set(${libraries_to_add} "${debug_libraries}" PARENT_SCOPE)
    else()
      message("Adding optimized libraries")
      set(${libraries_to_add} "${opt_libraries}" PARENT_SCOPE)
    endif()
  else()
    set(${libraries_to_add} "${library_list}" PARENT_SCOPE)
  endif()

endif()
endfunction()





macro(_add_to_link_line)
  set_property(GLOBAL APPEND PROPERTY CCSE_LINK_LINE ${ARGV})
endmacro()

macro(_add_to_target_list)
  set_property(GLOBAL APPEND PROPERTY CCSE_LIBRARY_TARGETS ${ARGV})
endmacro()



macro(_add_directories_to_link_line)
  foreach(directory ${ARGV})
    _add_to_link_line("-L${directory}")
  endforeach()
endmacro()

macro(_add_libraries_to_link_line)
  foreach(library ${ARGV})
    if(EXISTS ${library})  
      _add_to_link_line("${library}")   # If it's a filename, add it as given.
    else()
      _add_to_link_line("-l${library}") # Else, add it as a library to be looked up.
    endif()
  endforeach()
endmacro()



macro(add_package_libraries)

  # Grab the project name to find the dependent libraries
  SET(package ${PROJECT_NAME})

  # Add the directory locations of libraries it depends on.
  _add_directories_to_link_line("${${package}_LIBRARY_DIR}")
  _add_directories_to_link_line("${${package}_LIBRARY_DIRS}")

  # ${package}_LIBRARIES may contain debug and opt keywords, so parse the list into to_add:
  _parse_add_libraries("${${package}_LIBRARIES}" to_add)

  _add_libraries_to_link_line("${to_add}")

  # Add the accumulated contents of value_list to the link line property
  set_property(GLOBAL APPEND PROPERTY CCSE_LINK_LINE ${value_list})

endmacro()



macro(add_ccse_libraries libraries)
  _add_libraries_to_link_line(${libraries})
  _add_to_target_list(${libraries})
endmacro()

