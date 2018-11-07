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

   # Configure command
   set(BLITZ_CONFIG_COMMAND 
      CC=${CMAKE_C_COMPILER}
      FC=${CMAKE_Fortran_COMPILER}
      CXX=${CMAKE_CXX_COMPILER}
      ./configure --prefix=${BLITZ_INSTALL_DIR}
      )
   # if (DDEBUG)
   #    set(BLITZ_CONFIGURE_COMMAND ${BLITZ_CONFIGURE_COMMAND} --enable-debug)
   # else ()
   #    set(BLITZ_CONFIGURE_COMMAND ${BLITZ_CONFIGURE_COMMAND} --enable-optimize)
   # endif ()

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
   
   execute_process(
      COMMAND           ${BLITZ_ROOT_DIR}/configure --prefix=${BLITZ_INSTALL_DIR}
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

   # Fix source code
   execute_process( 
      COMMAND           sed -i /tinyvec-et.h/d  ${ALGOIM_INSTALL_DIR}/src/algoim_blitzinc.hpp 
      WORKING_DIRECTORY ${ALGOIM_INSTALL_DIR}
      RESULT_VARIABLE   RVAR
      OUTPUT_QUIET
      )
 
   if (NOT "${RVAR}" STREQUAL "0")
      message(FATAL_ERROR "Fatal error when fixing ALGOIM source code")
   endif()

endif ()   

   set(  BLITZ_INSTALL_DIR "" CACHE PATH
      "Path to Blitz installation directory (leave empty for superbuild)")
