#
# Fist, define the library we want to add
# Think of this as a constructor for object "amrex"
# In the following, we will use setters to change
# the properties of this object, like setting the sources,
# setting the compile definitions and so on
#
add_library ( amrex "" )

#
# Next, we implement a wrapper to add sources and headers
# to target object "amrex".
# If headers are detected, this function add the
# header's path to the include paths
# Variable PUBLIC_HEADERS accumulate all the headers for
# final installation
# WARNING: this routine add automatically ${CMAKE_CURRENT_LIST_DIR}
# to the input source file, so always add sources with the relative path
# to the directory where add_sources is called from
set (PUBLIC_HEADERS)
macro ( add_sources )
   foreach ( item IN ITEMS ${ARGV} )
      target_sources ( amrex PRIVATE ${CMAKE_CURRENT_LIST_DIR}/${item})
      # If it's a header, add it the public header group so it gets installed
      get_filename_component ( ext ${item} EXT )
      if ( ( ${ext} STREQUAL ".H" ) OR ( ${ext} STREQUAL ".h" ) )
	 # Do the following command to handle the case when ${item} is not the name
	 # of the source file, but it's in the form <relative-path>/<file-name>
	 # where <relative-path> is the path relative to the directory containing
	 # the list file invoking add_sources
	 get_filename_component ( dir ${CMAKE_CURRENT_LIST_DIR}/${item} DIRECTORY )
	 list ( APPEND PUBLIC_HEADERS ${CMAKE_CURRENT_LIST_DIR}/${item} )
	 target_include_directories ( amrex PUBLIC $<BUILD_INTERFACE:${dir}> )
      endif()      
   endforeach ()
   unset (ext)
   unset (dir)
endmacro ()

#
# Set some properties for target 'amrex', like
# compile definitions and so on
#

# General configuration
include ( AMReX_Config )
configure_amrex ()

#
# Now, one by one, let's add all the source files in the subdirectories
#
include (Base/CMakeLists.txt)
include (Boundary/CMakeLists.txt)
include (AmrCore/CMakeLists.txt)
include (Amr/CMakeLists.txt)

#
# Optional components
#
if (ENABLE_EB)
   include (EB/CMakeLists.txt)
endif ()

if (ENABLE_AMRDATA)
   include(Extern/amrdata/CMakeLists.txt)
endif()

if (ENABLE_PROFPARSER)
   include(Extern/ProfParser/CMakeLists.txt)
endif ()

if (ENABLE_LINEAR_SOLVERS)
   include(LinearSolvers/CMakeLists.txt)
endif ()

if (ENABLE_3D_NODAL_MLMG)
   include(Extern/Algoim/CMakeLists.txt)
endif ()

if (ENABLE_FORTRAN_INTERFACES) 
   include (F_Interfaces/CMakeLists.txt)
endif ()

if ( ENABLE_PARTICLES )
   include (Particle/CMakeLists.txt)
endif ()

if (ENABLE_SENSEI_INSITU)
   include (Extern/SENSEI/CMakeLists.txt)
endif ()


#
# Here we generate AMReX_BuildInfo.cpp
# If Python < 2.7, script won't work so
# we do not include this in library
#
#
# Check for Python >= 2.7: this is needed to build AMReX_BuildInfo.cpp
# This is not a required package. This provides PYTHONINTERP_FOUND and
# PYTHON_VERSION_STRING, PYTHON_VERSION_MAJOR and PYTHON_VERSION_MINOR  
#
find_package (PythonInterp QUIET)

if ( (NOT (PYTHON_VERSION_STRING VERSION_LESS "2.7") ) AND PYTHONINTERP_FOUND )
   add_custom_command(
      COMMAND  ${PROJECT_SOURCE_DIR}/Tools/C_scripts/makebuildinfo_C.py
      --amrex_home "${PROJECT_SOURCE_DIR}"
      --COMP ${CMAKE_C_COMPILER_ID} --COMP_VERSION ${CMAKE_C_COMPILER_VERSION}
      --FCOMP ${CMAKE_Fortran_COMPILER_ID} --FCOMP_VERSION ${CMAKE_C_COMPILER_VERSION}
      #--GIT ". ${AMREX_GIT}"
      OUTPUT AMReX_buildInfo.cpp
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating AMReX_buildInfo.cpp" )
   
   target_sources ( amrex PRIVATE  ${CMAKE_CURRENT_BINARY_DIR}/AMReX_buildInfo.cpp )
endif ()

#
# Set property to export public header
#
set_target_properties ( amrex PROPERTIES PUBLIC_HEADER "${PUBLIC_HEADERS}")


#
# Write and install configure file 
#
include(CMakePackageConfigHelpers)

# File for install tree
configure_package_config_file (
   ${AMREX_CMAKE_MODULES_PATH}/AMReXConfig.cmake.in
   ${CMAKE_BINARY_DIR}/export/AMReXConfig.cmake
   INSTALL_DESTINATION bin/cmake/AMReX )

# # File for build tree
# configure_package_config_file (
#    ${AMREX_CMAKE_MODULES_PATH}/AMReXConfig.cmake.in
#    ${CMAKE_BINARY_DIR}/AMReXConfig.cmake
#    INSTALL_DESTINATION bin/cmake/AMReX )

# 
# Package version is a modified form of AMREX_GIT_VERSION
#
if (AMREX_GIT_VERSION)
   string (FIND "${AMREX_GIT_VERSION}" "-" IDX REVERSE)
   string (SUBSTRING "${AMREX_GIT_VERSION}" 0 "${IDX}" AMREX_VERSION)
   string (FIND "${AMREX_VERSION}" "-" IDX REVERSE)
   string (SUBSTRING "${AMREX_VERSION}" 0 "${IDX}" AMREX_VERSION )
   string (REPLACE "-" "." AMREX_VERSION "${AMREX_VERSION}")
endif ()

write_basic_package_version_file ( ${CMAKE_BINARY_DIR}/AMReXConfigVersion.cmake
   VERSION ${AMREX_VERSION}
   COMPATIBILITY AnyNewerVersion )



install ( FILES
   ${CMAKE_BINARY_DIR}/export/AMReXConfig.cmake
   ${CMAKE_BINARY_DIR}/AMReXConfigVersion.cmake
   DESTINATION lib/cmake/AMReX ) 

#
# Define install  
#
install ( TARGETS amrex
   EXPORT        AMReXTargets
   ARCHIVE       DESTINATION lib 
   LIBRARY       DESTINATION lib 
   INCLUDES      DESTINATION include # Adds proper directory to INTERFACE_INCLUDE_DIRECTORIES 
   PUBLIC_HEADER DESTINATION include )

install ( EXPORT AMReXTargets
   NAMESPACE AMReX::
   DESTINATION lib/cmake/AMReX )

# Install fortran modules
get_target_property (AMREX_Fortran_MODULE_DIR amrex Fortran_MODULE_DIRECTORY )
install ( DIRECTORY ${AMREX_Fortran_MODULE_DIR}/ # Trailing backslash is crucial here!
   DESTINATION include )

# This header in a weird path has to be copied to install includes
install ( FILES ${PROJECT_SOURCE_DIR}/Tools/C_scripts/AMReX_buildInfo.H
   DESTINATION include )

# Install Tools directory
install (DIRECTORY ${PROJECT_SOURCE_DIR}/Tools/
   DESTINATION Tools 
   USE_SOURCE_PERMISSIONS )



# #Export targets in build-tree ( so third-party can use amrex without installing???)
# #Check this

# export ( TARGETS
#          amrex
# 	 NAMESPACE AMReX::
# 	 FILE "${CMAKE_BINARY_DIR}/AMReXTargets.cmake" )

