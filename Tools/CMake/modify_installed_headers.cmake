r#
# Add line #include<AMReX_Config.H> to all installed headers
# 
message(STATUS "Modifying installed headers")
file(GLOB _installed_headers "${CMAKE_INSTALL_PREFIX}/include/*.H")
list(REMOVE_ITEM _installed_headers "${CMAKE_INSTALL_PREFIX}/include/AMReX_Config.H")

foreach( _header IN LISTS _installed_headers)
   file(READ ${_header} _contents)
   string(PREPEND _contents "#include <AMReX_Config.H>\n")
   file(WRITE ${_header} "${_contents}")
endforeach()
