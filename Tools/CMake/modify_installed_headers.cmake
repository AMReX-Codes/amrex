#
# Prepend line #include<AMReX_Config.H> to all installed headers
#
message(STATUS "Adding #include <AMReX_Config.H> to installed headers")
file(GLOB _installed_headers "${CMAKE_INSTALL_PREFIX}/include/*.H")
list(REMOVE_ITEM _installed_headers "${CMAKE_INSTALL_PREFIX}/include/AMReX_Config.H")

foreach( _header IN LISTS _installed_headers )
  file(READ "${_header}" HEADER_CONTENTS)
  file(WRITE "${_header}" "#include <AMReX_Config.H>\n")
  file(APPEND "${_header}" "${HEADER_CONTENTS}")
endforeach()
