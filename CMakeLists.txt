# -*- mode: cmake -*-
# Require cmake 2.8.8 or higher
cmake_minimum_required(VERSION 2.8.8)

#set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
#set(CMAKE_DISABLE_SOURCE_CHANGES ON)


if ("${CMAKE_PATCH_VERSION}" GREATER 3)
  set(ADJUST_POLICY "True")
endif()

if (${ADJUST_POLICY})
  cmake_policy(SET CMP0017 NEW)
endif()

# Enable testing, ctest needs this
# all add_test commands are ignored unless this is called!
enable_testing()

# Useful variables pointing to directories in the source tree
#set(CCSE_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(CCSE_TOOLS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Tools")
set(CCSE_PERL_DIR "${CCSE_TOOLS_DIR}/C_scripts")
#set(CCSE_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/install")
set(CCSE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}")

# Set the module search path so find_package and include commands
# can locate files in <root source tree>/tools/cmake
set(CCSE_MODULE_PATH "${CCSE_TOOLS_DIR}/CMake")
set(CMAKE_MODULE_PATH  ${CCSE_MODULE_PATH})

# Addtional build options
include(CCSEOptions)

# A property for accumulating the a global CCSE link line.
set_property(GLOBAL PROPERTY CCSE_LINK_LINE "-L${CMAKE_INSTALL_PREFIX}/lib")

# A property for accumulating CCSE library targets
set_property(GLOBAL PROPERTY CCSE_LIBRARY_TARGETS)
set(CCSE_LINK_LINE_FILE "${CMAKE_INSTALL_PREFIX}/lib/link_line")  # A filename to write link-line to.
include(InstallManager)


# Define the project name
# This command will define
#     CCSE_SOURCE_DIR
#     CCSE_BINARY_DIR
project(CCSE)

# Source files for all binaries and libraries found under src
add_subdirectory(Src)

# Make the exports only _after_ doing the build
create_exports()

ADD_INSTALL_PERLSCRIPT(${CCSE_PERL_DIR}/strip72)

option(ENABLE_Config_Report "Print out detailed information at the end of a configuration")
set(CCSE_CONFIG_LOG "${CCSE_BINARY_DIR}/ccse-config.log"
    CACHE string "CCSE configuration log file")
include(CCSEConfigReport)
if (ENABLE_Config_Report)
    execute_process(COMMAND "${cat_exec}" "${CCSE_CONFIG_LOG}" OUTPUT_VARIABLE config_output)
    message(STATUS "contents of: ${CCSE_CONFIG_LOG}")
    message(STATUS "${config_output}")
endif()    

