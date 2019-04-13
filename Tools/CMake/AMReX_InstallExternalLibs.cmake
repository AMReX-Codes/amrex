#
# Options provided to the user
#
set(BLITZ_INSTALL_PREFIX ""  CACHE PATH "Path to Blitz installation directory")
set(ALGOIM_INSTALL_PREFIX "" CACHE PATH "Path to Algoim installation directory")

if (NOT TARGET amrex)
   message(FATAL_ERROR "Target amrex must be defined before including AMReX_InstallExternalLibs.cmake")
endif ()

#
# Some variables for internal use
#
set(_EXTERNAL_LIBS_PATH "${CMAKE_BINARY_DIR}/external"
   CACHE INTERNAL "Installation path for external libraries (blitz and Algoim)")

file( MAKE_DIRECTORY ${_EXTERNAL_LIBS_PATH} )

#
# Build and install blitz if required.
# If blitz is already installed and the path of the root installation
# is given via BLITZ_INSTALL_DIR, search for blitz
#
if (NOT BLITZ_INSTALL_PREFIX)

   # Set internal variables for blitz installation
   set(_BLITZ_REPO   "https://github.com/mic84/blitz.git" CACHE INTERNAL "Blitz git repo url")
   set(_BLITZ_GIT_TAG master CACHE INTERNAL "Blitz git tag")
   set(_BLITZ_ROOT_DIR  "${_EXTERNAL_LIBS_PATH}/blitz" CACHE INTERNAL "Blitz root directory" )
   set(_BLITZ_BUILD_DIR "${_BLITZ_ROOT_DIR}/builddir"   CACHE INTERNAL "Blitz build directory"  )
   set(_BLITZ_INSTALL_DIR "${_BLITZ_ROOT_DIR}/installdir" CACHE INTERNAL "Blitz install directory")

   if (NOT _BLITZ-INSTALLED)

      if (EXISTS ${PROJECT_SOURCE_DIR}/Src/Extern/blitz)
         file( REMOVE_RECURSE  ${PROJECT_SOURCE_DIR}/Src/Extern/blitz)
      endif ()

      # Clone
      message(STATUS "Cloning Blitz")
      execute_process(
         COMMAND           git clone -q ${_BLITZ_REPO} # -q to have it quiet
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/Src/Extern #${_EXTERNAL_LIBS_PATH}
         RESULT_VARIABLE   _RVAR
         OUTPUT_QUIET
         )

      if (NOT "${_RVAR}" STREQUAL "0")
         message(FATAL_ERROR "Fatal error when cloning Blitz repo")
      endif()

      # Checkout correct version
      message(STATUS "Checking out Blitz ${_BLITZ_GIT_TAG}")
      execute_process(
         COMMAND           git checkout -q ${_BLITZ_GIT_TAG}
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/Src/Extern/blitz #${_BLITZ_ROOT_DIR}
         RESULT_VARIABLE   _RVAR
         OUTPUT_QUIET
         )

      if (NOT "${_RVAR}" STREQUAL "0")
         message(FATAL_ERROR "Fatal error when checking out Blitz ${_BLITZ_GIT_TAG}")
      endif()

      set(_BLITZ-INSTALLED TRUE CACHE INTERNAL "Blitz has been installed")

   endif ()

   add_subdirectory(${PROJECT_SOURCE_DIR}/Src/Extern/blitz)
   target_link_libraries(amrex PUBLIC blitz)

else ()

   # Try to find the package
   find_package(blitz HINTS ${BLITZ_INSTALL_PREFIX})

   if (NOT blitz_FOUND)  # If blitz is too old to have cmake enabled

      set(_NO_BLITZ_DEPENDENCY_NEEDED TRUE CACHE INTERNAL "Weather config file needs find_dependency for blitz")

      # Try package
      find_package(PkgConfig)
      if (PkgConfig_FOUND)
         set(CMAKE_PREFIX_PATH ${BLITZ_INSTALL_PREFIX})
         pkg_check_modules(BLITZ IMPORTED_TARGET blitz)
      endif ()

      if (BLITZ_FOUND) # This is true if both pkgconfig is found and blitz is found
         get_target_property(_BLITZ_LIBRARIES    PkgConfig::BLITZ INTERFACE_LINK_LIBRARIES)
         get_target_property(_BLITZ_INCLUDE_DIR  PkgConfig::BLITZ INTERFACE_INCLUDE_DIRECTORIES)
         get_target_property(_BLITZ_DEFINES      PkgConfig::BLITZ INTERFACE_COMPILE_DEFINITIONS)
      else () # manual search
         message( STATUS "Trying to find blitz \"manually\" ")
         find_library(_BLITZ_LIBRARIES NAMES libblitz.a HINTS ${BLITZ_INSTALL_PREFIX})
         find_path(_BLITZ_INCLUDE_DIR  NAMES blitz      HINTS ${BLITZ_INSTALL_PREFIX})

         include(CheckIncludeFileCXX)
         set(CMAKE_REQUIRED_INCLUDES ${_BLITZ_INCLUDE_DIR})
         check_include_file_cxx("blitz/blitz.h" HAVE_BLITZ_H)
         if (  ${_BLITZ_LIBRARIES}   STREQUAL "_BLITZ_LIBRARIES-NOTFOUND" OR
               ${_BLITZ_INCLUDE_DIR} STREQUAL "_BLITZ_INCLUDE_DIR-NOTFOUND" OR
               NOT HAVE_BLITZ_H)
            message(FATAL_ERROR "Cannot locate blitz")
         endif ()
      endif ()
      target_include_directories(amrex PUBLIC ${_BLITZ_INCLUDE_DIR})
      target_link_libraries(amrex PUBLIC ${_BLITZ_LIBRARIES})
   else()
      target_link_libraries(amrex PUBLIC blitz)
   endif()

endif()

#
# Download algoim if required
#
if ( NOT ALGOIM_INSTALL_PREFIX )

   set(_ALGOIM_REPO  "https://github.com/algoim/algoim.git" CACHE INTERNAL "Algoim git repo url")
   set(_ALGOIM_GIT_TAG  "a3d0b7bb2872cd414f77dbe7e77b25b9e707eaf3" CACHE INTERNAL "Algoim git tag")
   set(_ALGOIM_INSTALL_DIR "${_EXTERNAL_LIBS_PATH}/algoim" CACHE INTERNAL "Algoim install directory")

   if (NOT _ALGOIM-INSTALLED )

      # Clone
      message(STATUS "Cloning Algoim")
      execute_process(
         COMMAND           git clone -q ${_ALGOIM_REPO} # -q to have it quiet
         WORKING_DIRECTORY ${_EXTERNAL_LIBS_PATH}/
         RESULT_VARIABLE   _RVAR
         OUTPUT_QUIET
         )

      if (NOT "${_RVAR}" STREQUAL "0")
         message(FATAL_ERROR "Fatal error when cloning ALGOIM repo")
      endif()

      # Fix source code NOTE: on Mac, BSD sed requires backup suffix for -i option
      # https://stackoverflow.com/questions/7573368/in-place-edits-with-sed-on-os-x
      execute_process(
         COMMAND           sed -i'' /tinyvec-et.h/d  ${_ALGOIM_INSTALL_DIR}/src/algoim_blitzinc.hpp
         WORKING_DIRECTORY ${_ALGOIM_INSTALL_DIR}
         RESULT_VARIABLE   _RVAR
         OUTPUT_QUIET
         )

      if (NOT "${_RVAR}" STREQUAL "0")
         message(FATAL_ERROR "Fatal error when fixing ALGOIM source code.")
      endif()

      set(_ALGOIM-INSTALLED TRUE CACHE INTERNAL "Algoim has been installed")

   endif ()

   target_include_directories(amrex PUBLIC
      $<BUILD_INTERFACE:${_ALGOIM_INSTALL_DIR}/src>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include/algoim>
      )

   file(GLOB_RECURSE _ALGOIM_HEADERS  "${_ALGOIM_INSTALL_DIR}/src/*.hpp" )
   install(FILES ${_ALGOIM_HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/algoim)

else ()

   set_target_properties(algoim PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES ${ALGOIM_INSTALL_PREFIX}/src )

endif ()
