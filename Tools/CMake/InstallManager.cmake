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
# Usage: create_exports
#
# Arguments: None
#
#
function (CREATE_EXPORTS)

# Template file located in the CMake module directory

# Find the packages found for Amanzi
get_property(LINK_LINE GLOBAL PROPERTY CCSE_LINK_LINE)

# Convert the link line to a space deliminated string
foreach (arg ${LINK_LINE})
  set(LINK_LINE_STRING "${LINK_LINE_STRING} ${arg}")
endforeach()

# Write and install the link-line file
file(WRITE ${CCSE_LINK_LINE_FILE} ${LINK_LINE_STRING})
install(FILES ${CCSE_LINK_LINE_FILE} DESTINATION lib)

# Write the export Makefile and add to the include install list
set(in_makefile  "${CCSE_MODULE_PATH}/MakefileConfig.export.in")
set(out_makefile "${CCSE_BINARY_DIR}/Makefile.export")
configure_file("${in_makefile}" "${out_makefile}")
install(FILES "${out_makefile}" DESTINATION lib)

# Write the CCSEConfig.cmake file
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
