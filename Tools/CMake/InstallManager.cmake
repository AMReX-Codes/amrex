# # -*- mode: cmake -*-

#
# Functions for managing the install targets
#


include(CMakeParseArguments)
include(CCSELinkLine)

export(PACKAGE CCSE)

#
# Usage: ADD_INSTALL_INCLUDE_FILE( file1 file2 file3 ... )
#
# Arguments:
#  A list of files that will be installed in the CCSE_INSTALL_INCLUDE_DIR
#
#
function ( ADD_INSTALL_INCLUDE_FILE )

  foreach(_inc_file ${ARGV})
    install(
      FILES ${_inc_file}
      DESTINATION include
      )
  endforeach()

endfunction( ADD_INSTALL_INCLUDE_FILE )

#
# Usage: ADD_INSTALL_LIBRARY( lib1 lib2 lib3 ... )
#
# Arguments:
#  A list of libraries that will be installed in the CCSE_INSTALL_LIB_DIR
#
#
function ( ADD_INSTALL_LIBRARY )

  install(
    TARGETS ${ARGV}
    EXPORT CCSETargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

  # Add the libraries to our global list
  add_ccse_libraries(${ARGV})

  # Add dependency libaries as determined by the pacakge definition.
  add_package_libraries()

endfunction( ADD_INSTALL_LIBRARY )


#
# Usage: ADD_INSTALL_SHELLSCRIPT( script1 script2 ... )
#
# Arguments:
#  A list of shell scripts that will be installed in the CCSE_INSTALL_BIN_DIR
#
#
function ( ADD_INSTALL_SHELLSCRIPT )

foreach(_shellscript_file ${ARGV})
  install(
    FILES ${_shellscript_file}
    DESTINATION bin
    )
endforeach()

endfunction( ADD_INSTALL_SHELLSCRIPT )


#
# Usage: ADD_INSTALL_PERLSCRIPT( script1 script2 ... )
#
# Arguments:
#  A list of perl scripts that will be installed in the CCSE_INSTALL_PERL_DIR
#
#
function ( ADD_INSTALL_PERLSCRIPT )

foreach(_perlscript_file ${ARGV})
  install(
    FILES ${_perlscript_file}
    DESTINATION perl
    )
endforeach()

endfunction( ADD_INSTALL_PERLSCRIPT )


#
# Usage: ADD_INSTALL_CMAKE_FILES( file1 file2 ... )
#
# Arguments:
#  A list of cmake module files that will be installed in CCSE_INSTALL_PREFIX/cmake
#
#
function ( ADD_INSTALL_CMAKE_FILES )

foreach(_cmake_file ${ARGV})
  install(
    FILES ${_cmake_file}
    DESTINATION cmake
    )
endforeach()

endfunction( ADD_INSTALL_CMAKE_FILES )



#
# Usage: ADD_INSTALL_BINARY( exe1 exe2 ... )
#
# Arguments:
#  A list of executables that will be installed in the CCSE_INSTALL_BIN_DIR
#
#
function ( ADD_INSTALL_BINARY )

foreach(_bin_file ${ARGV})
  install(
    TARGETS ${_bin_file}
    EXPORT CCSETargets
    DESTINATION bin
    )
endforeach()

endfunction( ADD_INSTALL_BINARY )


#
# Usage: makefile_include_dirs(CMAKE_INCLUDE_LIST in_list
#                               MAKE_INCLUDE_LIST out_list)
#
# Arguments:
#          CMAKE_INCLUDE_LIST List of include directories
#          MAKE_INCLUDE_LIST  List include directories for make
#
function(makefile_include_dirs)

    cmake_parse_arguments(PARSE_ARGS "" "MAKE_INCLUDE_LIST" "CMAKE_INCLUDE_LIST" ${ARGN})
    #print_variable(PARSE_ARGS_CMAKE_INCLUDE_LIST)
    #print_variable(PARSE_ARGS_MAKE_INCLUDE_LIST)

    set(tmp_inc_list)
    set(loop_list ${PARSE_ARGS_CMAKE_INCLUDE_LIST})
    list(REMOVE_DUPLICATES loop_list)
    foreach( dir  ${loop_list})
      set(i_path "-I${dir} ")
      list(APPEND tmp_inc_list ${i_path})
    endforeach()

    set(tmp_make_list)
    string(REGEX REPLACE ";" "" tmp_make_list ${tmp_inc_list})
    set(${PARSE_ARGS_MAKE_INCLUDE_LIST} "${tmp_make_list}" PARENT_SCOPE)

endfunction(makefile_include_dirs)

#
# Usage: makefile_library_dirs(CMAKE_LIB_LIST in_list
#                              MAKE_LIB_LIST out_list)
#
# Arguments:
#          CMAKE_LIB_LIST List of library directories
#          MAKE_LIB_LIST  List library directories for make
#
function(makefile_library_dirs)

    cmake_parse_arguments(PARSE_ARGS "" "MAKE_LIB_LIST" "CMAKE_LIB_LIST" ${ARGN})
    set(tmp_lib_list)
    set(loop_list ${PARSE_ARGS_CMAKE_LIB_LIST})
    list(REVERSE loop_list)
    list(REMOVE_DUPLICATES loop_list)
    list(REVERSE loop_list)
    foreach( dir  ${loop_list})
      set(l_path "-L${dir} ")
      list(APPEND tmp_lib_list ${l_path})
    endforeach()

    set(tmp_make_list)
    string(REGEX REPLACE ";" "" tmp_make_list ${tmp_lib_list})
    set(${PARSE_ARGS_MAKE_LIB_LIST} "${tmp_make_list}" PARENT_SCOPE)

endfunction(makefile_library_dirs)

#
# Usage: create_exports
#
# Arguments: None
#
#
function (CREATE_EXPORTS)

# Template file located in the CMake module directory

# Find the packages found for CCSE
get_property(LINK_LINE GLOBAL PROPERTY CCSE_LINK_LINE)

# Define CCSE_INCLUDE_DIRS and CCSE_LIBRARY_DIRS
set(CCSE_INCLUDE_DIRS "${CMAKE_INSTALL_PREFIX}/include")
set(CCSE_LIBRARY_DIRS "${CMAKE_INSTALL_PREFIX}/lib")

list(REMOVE_DUPLICATES CCSE_INCLUDE_DIRS)
list(REMOVE_DUPLICATES CCSE_LIBRARY_DIRS)

# Write the export Makefile and add to the include install list
makefile_include_dirs(CMAKE_INCLUDE_LIST ${CCSE_INCLUDE_DIRS}
                      MAKE_INCLUDE_LIST CCSE_MAKE_INCLUDE_DIRS)
makefile_library_dirs(CMAKE_LIB_LIST ${CCSE_LIBRARY_DIRS}
                      MAKE_LIB_LIST CCSE_MAKE_LIBRARY_DIRS)
set(in_makefile  "${CCSE_MODULE_PATH}/MakefileConfig.export.in")
set(out_makefile "${CCSE_BINARY_DIR}/Makefile.export")
configure_file("${in_makefile}" "${out_makefile}")
install(FILES "${out_makefile}" DESTINATION lib)

set(in_config   "${CCSE_MODULE_PATH}/CCSEConfig-install.cmake.in")
set(out_config   "${CCSE_BINARY_DIR}/CCSEConfig.cmake")
configure_file(${in_config} ${out_config})
install(FILES ${out_config} DESTINATION lib)

# Write the CCSEConfigVersion.cmake file
set(in_config   "${CCSE_MODULE_PATH}/CCSEConfigVersion-install.cmake.in")
set(out_config   "${CCSE_BINARY_DIR}/CCSEConfigVersion.cmake")
configure_file(${in_config} ${out_config} @ONLY)
install(FILES ${out_config} DESTINATION lib)

# Write the CMake configuration target file
message(STATUS "Writing target file")
install(EXPORT CCSETargets
        DESTINATION lib
	NAMESPACE ccse_
	FILE CCSETargets.cmake)

endfunction()
