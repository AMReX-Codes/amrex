#
# Add line #include<AMReX_Config.H> to all installed headers
# 
message(STATUS "Modifying installed headers")
file(GLOB _installed_headers "${CMAKE_INSTALL_PREFIX}/include/*.H")
list(REMOVE_ITEM _installed_headers "${CMAKE_INSTALL_PREFIX}/include/AMReX_Config.H")

foreach( _header IN LISTS _installed_headers)
   # This command only affect the first line of each file. Therefore, it can be used countless
   # times on the same file and it always gives the same result (i.e. it does not keep prepending)
   execute_process(COMMAND bash "-c" "sed -i '1s/^/#include <AMReX_Config.H>\\n\\n/' ${_header}" )
endforeach()
