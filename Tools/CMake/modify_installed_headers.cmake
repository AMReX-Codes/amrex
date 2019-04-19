#
# Add line #include<AMReX_Config.H> to all installed headers
#
message(STATUS "Adding #include <AMReX_Config.H> to installed headers")
file(GLOB _installed_headers "${CMAKE_INSTALL_PREFIX}/include_orig/*.H")

list(REMOVE_ITEM _installed_headers "${CMAKE_INSTALL_PREFIX}/include/AMReX_Config.H")

foreach( _header IN LISTS _installed_headers )
  get_filename_component(_header_fname "${_header}" NAME)
  set(_header_dest "${CMAKE_INSTALL_PREFIX}/include/${_header_fname}")
  file(READ "${_header}" HEADER_CONTENTS)
  file(WRITE "${_header_dest}" "#include <AMReX_Config.H>\n")
  file(APPEND "${_header_dest}" "${HEADER_CONTENTS}")
endforeach()
