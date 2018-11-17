set(EXTERNAL_LIBS_PATH ${CMAKE_BINARY_DIR}/external)
file(MAKE_DIRECTORY ${EXTERNAL_LIBS_PATH})

#
# Build and install blitz if required
#
if (NOT BLITZ_INSTALL_DIR)

   set(BLITZ_REPO        https://github.com/blitzpp/blitz.git     )
   set(BLITZ_GIT_TAG     0d5d1feaad1b1b31aa9881da00bad68f7036ff63 )
   set(BLITZ_ROOT_DIR    ${EXTERNAL_LIBS_PATH}/blitz              )
   set(BLITZ_INSTALL_DIR ${BLITZ_ROOT_DIR}/installdir CACHE PATH
      "Path to Blitz installation directory")

   # Clone
   message(STATUS "Cloning Blitz")
   execute_process(
      COMMAND           git clone -q ${BLITZ_REPO} # -q to have it quiet
      WORKING_DIRECTORY ${EXTERNAL_LIBS_PATH}
      RESULT_VARIABLE   RVAR
      OUTPUT_QUIET
      )

   if (NOT "${RVAR}" STREQUAL "0")
      message(FATAL_ERROR "Fatal error when cloning BLITZ repo")
   endif()

   # Configure
   message(STATUS "Configuring Blitz")
   execute_process(
      COMMAND           git checkout ${BLITZ_GIT_TAG}
      COMMAND           autoreconf -fiv
      OUTPUT_FILE       ${BLITZ_ROOT_DIR}/cmake_autoreconf.log
      ERROR_FILE        ${BLITZ_ROOT_DIR}/cmake_autoreconf.error
      WORKING_DIRECTORY ${BLITZ_ROOT_DIR}
      RESULT_VARIABLE   RVAR
      )

   if (NOT "${RVAR}" STREQUAL "0")
      message(FATAL_ERROR "Fatal error when running autoreconf for BLITZ ")
   endif()

   if (DDEBUG)
      set(BLITZ_BUILD_TYPE_FLAG --enable-debug)
   else ()
      set(BLITZ_BUILD_TYPE_FLAG --enable-optimize)
   endif ()

   # If would be "more correct" to define FC, CXX, and CC before ./configure
   # and not after. However, this does not work with execute_process.
   # Nonetheless, this should not be a problem since many 'configure' implementation
   # support this approach.
   # Also, we are using a hack to compile on cray machines because blitz configure
   # setup sucks.
   if ("${SITE}" STREQUAL "nersc")
      # Fortran
      if ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
         set(FC ifort)
      elseif ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
         set(FC gfortran)
      elseif ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
         set(FC pgf90)
      endif ()

      # C++
      if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
         set(CXX icpc)
      elseif ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
         set(CXX g++)
      elseif ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "PGI")
         set(CXX pgc++)
      endif ()

      # C
      if ( ${CMAKE_C_COMPILER_ID} STREQUAL "Intel")
         set(CC icc)
      elseif ( ${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
         set(CC gcc)
      elseif ( ${CMAKE_C_COMPILER_ID} STREQUAL "PGI")
         set(CC pgcc)
      endif ()
   else ()
      set(FC  ${CMAKE_Fortran_COMPILER} )
      set(CXX ${CMAKE_CXX_COMPILER} )
      set(CC  ${CMAKE_C_COMPILER} )
   endif ()

   execute_process(
      COMMAND           ${BLITZ_ROOT_DIR}/configure
                        CC=${CC} FC=${FC} CXX=${CXX}
                        --prefix=${BLITZ_INSTALL_DIR} ${BLITZ_BUILD_TYPE_FLAG}
      OUTPUT_FILE       ${BLITZ_ROOT_DIR}/cmake_configure.log
      ERROR_FILE        ${BLITZ_ROOT_DIR}/cmake_configure.error
      WORKING_DIRECTORY ${BLITZ_ROOT_DIR}
      RESULT_VARIABLE   RVAR
      )

   if (NOT "${RVAR}" STREQUAL "0")
      message(FATAL_ERROR "Fatal error when running configure for BLITZ ")
   endif()

   # Install
   message(STATUS "Installing Blitz")
   execute_process(
      COMMAND           make install
      OUTPUT_FILE       ${BLITZ_ROOT_DIR}/cmake_install.log
      ERROR_FILE        ${BLITZ_ROOT_DIR}/cmake_install.error
      WORKING_DIRECTORY ${BLITZ_ROOT_DIR}
      RESULT_VARIABLE   RVAR
      )

   if (NOT "${RVAR}" STREQUAL "0")
      message(FATAL_ERROR "Fatal error when installing BLITZ ")
   endif()

endif()


#
# Download algoim if required
#
if (NOT ALGOIM_INSTALL_DIR)

   set(ALGOIM_REPO        https://github.com/algoim/algoim.git     )
   set(ALGOIM_GIT_TAG     a3d0b7bb2872cd414f77dbe7e77b25b9e707eaf3 )
   set(ALGOIM_INSTALL_DIR ${EXTERNAL_LIBS_PATH}/algoim CACHE PATH
      "Path to Algoim installation directory")

   # Clone
   message(STATUS "Cloning Algoim")
   execute_process(
      COMMAND           git clone -q ${ALGOIM_REPO} # -q to have it quiet
      WORKING_DIRECTORY ${EXTERNAL_LIBS_PATH}/
      RESULT_VARIABLE   RVAR
      OUTPUT_QUIET
      )

   if (NOT "${RVAR}" STREQUAL "0")
      message(FATAL_ERROR "Fatal error when cloning ALGOIM repo")
   endif()

   # Fix source code NOTE: there is a odd problem using the macOS default sed:
   # https://stackoverflow.com/questions/7573368/in-place-edits-with-sed-on-os-x
   execute_process(
      COMMAND           sed -i /tinyvec-et.h/d  ${ALGOIM_INSTALL_DIR}/src/algoim_blitzinc.hpp
      WORKING_DIRECTORY ${ALGOIM_INSTALL_DIR}
      RESULT_VARIABLE   RVAR
      OUTPUT_QUIET
      )

   if (NOT "${RVAR}" STREQUAL "0")
      message(FATAL_ERROR "Fatal error when fixing ALGOIM source code. If you're using macOS, please install gnu-sed via homebrew.")
   endif()

endif ()

   set(  BLITZ_INSTALL_DIR "" CACHE PATH
      "Path to Blitz installation directory (leave empty for superbuild)")
