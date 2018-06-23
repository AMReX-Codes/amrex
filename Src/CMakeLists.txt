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
# to the source file given, so always add sources with the relative path
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
configure_amrex ()


# Include paths 
target_include_directories ( amrex
   PUBLIC
#   "${AMREX_EXTRA_Fortran_INCLUDE_PATH};${AMREX_EXTRA_CXX_INCLUDE_PATH}"
   $<INSTALL_INTERFACE:include> # This is the include dir in installation tree
   $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}> 
   )

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
   include (GeometryShop/CMakeLists.txt)
   include (EB2/CMakeLists.txt)
endif ()

if (ENABLE_PROFPARSER OR ENABLE_AMRDATA)
   include (Extern/CMakeLists.txt)
endif ()

if (ENABLE_LINEAR_SOLVERS)
   include (LinearSolvers/CMakeLists.txt)
endif ()

if (ENABLE_FORTRAN_INTERFACES) 
   include (F_Interfaces/CMakeLists.txt)
endif ()

if ( ENABLE_PARTICLES )
   include (Particle/CMakeLists.txt)
endif ()

#
# Here we generate AMReX_BuildInfo.cpp
# If Python < 2.7, script won't work so
# we do not include this in library
# 
if ( NOT (PYTHON_VERSION_STRING VERSION_LESS "2.7") )
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
set_target_properties ( amrex PROPERTIES PUBLIC_HEADER "${PUBLIC_HEADERS}" )

#
# Define install  
#
install ( TARGETS amrex EXPORT AMReX ARCHIVE DESTINATION ${AMREX_LIBRARY_DIR}
   LIBRARY DESTINATION ${AMREX_LIBRARY_DIR} PUBLIC_HEADER DESTINATION ${AMREX_INCLUDE_DIR} )

install ( EXPORT AMReX NAMESPACE AMReX:: DESTINATION ${AMREX_CMAKE_DIR} )

install ( DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${AMREX_INCLUDE_DIR})


#
# This header in a weird path has to be copied to install includes
# 
install ( FILES ${PROJECT_SOURCE_DIR}/Tools/C_scripts/AMReX_buildInfo.H
   DESTINATION ${AMREX_INCLUDE_DIR} )

#
# Write and install configure file 
# 
configure_file ( ${AMREX_CONFIG_INSTALL_INFILE} ${AMREX_CONFIG_INSTALL_OUTFILE} @ONLY)
install ( FILES ${AMREX_CONFIG_INSTALL_OUTFILE} DESTINATION ${AMREX_CMAKE_DIR})

#
# Install Tools directory
#
install (DIRECTORY ${PROJECT_SOURCE_DIR}/Tools
   DESTINATION ${CMAKE_INSTALL_PREFIX}
   USE_SOURCE_PERMISSIONS )

