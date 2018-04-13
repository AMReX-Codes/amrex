#
# This file provides the following variables
#
#   TPP_CFLAGS
#   TPP_FFLAGS
#   TPP_CXXFLAGS
#   TPP_C_INCLUDE_PATH
#   TPP_Fortran_INCLUDE_PATH
#   TPP_CXX_INCLUDE_PATH
#   TPP_C_LINK_LINE
#   TPP_Fortran_LINK_LINE
#   TPP_CXX_LINK_LINE
# 

#
# Check if AMReX_Options.cmake  and AMREX_CMakeVariables.cmake
# have been already processed
#
if ( NOT (DEFINED __AMREX_OPTIONS__ ) )
   message ( FATAL_ERROR "AMReX_Options.cmake must be\
included before AMReX_ThirdPartyProfilers.cmake" )
endif ()

if ( NOT (DEFINED __AMREX_CMAKEVARIABLES__ ) )
   message ( FATAL_ERROR "AMReX_CMakeVariables.cmake must be\
included before AMReX_ThirdPartyProfilers.cmake" )
endif ()


set ( TPP_CFLAGS )
set ( TPP_FFLAGS ) 
set ( TPP_CXXFLAGS )
set ( TPP_C_INCLUDE_PATH )
set ( TPP_Fortran_INCLUDE_PATH )
set ( TPP_CXX_INCLUDE_PATH )
set ( TPP_C_LINK_LINE )
set ( TPP_Fortran_LINK_LINE )
set ( TPP_CXX_LINK_LINE )


if ( ${TP_PROFILE} MATCHES "CRAYPAT" )

   if ( ${SITE} MATCHES "nersc" )
      list (APPEND TPP_Fortran_INCLUDE_PATH "$ENV{CRAYPAT_ROOT}/include")
      list (APPEND TPP_C_INCLUDE_PATH "$ENV{CRAYPAT_ROOT}/include")
      list (APPEND TPP_CXX_INCLUDE_PATH "$ENV{CRAYPAT_ROOT}/include")
   endif ()


elseif ( ${TP_PROFILE} MATCHES "FORGE" )

   if ( ${SITE} MATCHES "nersc" )
      list (APPEND TPP_Fortran_INCLUDE_PATH
	 "$ENV{ALLINEA_TOOLS_DIR}/$ENV{ALLINEA_TOOLS_VERSION}/map/wrapper")
      list (APPEND TPP_C_INCLUDE_PATH 
	 "$ENV{ALLINEA_TOOLS_DIR}/$ENV{ALLINEA_TOOLS_VERSION}/map/wrapper")
      list (APPEND TPP_CXX_INCLUDE_PATH 
	 "$ENV{ALLINEA_TOOLS_DIR}/$ENV{ALLINEA_TOOLS_VERSION}/map/wrapper")
   endif ()


elseif ( ${TP_PROFILE} MATCHES "VTUNE" )

   set ( TPP_CFLAGS "-debug inline-debug-info -parallel-source-info=2")
   set ( TPP_FFLAGS "-debug inline-debug-info -parallel-source-info=2")
   set ( TPP_CXXFLAGS "-debug inline-debug-info -parallel-source-info=2")

   if ( ${SITE} MATCHES "nersc" )
      set ( TPP_CFLAGS "-dynamic ${TPP_CFLAGS}")
      set ( TPP_FFLAGS "-dynamic ${TPP_FFLAGS}")
      set ( TPP_CXXFLAGS "-dynamic ${TPP_CXXFLAGS}")

      list ( APPEND TPP_Fortran_INCLUDE_PATH "$ENV{VTUNE_AMPLIFIER_XE_2018_DIR}/include" )
      list ( APPEND TPP_C_INCLUDE_PATH "$ENV{VTUNE_AMPLIFIER_XE_2018_DIR}/include" )
      list ( APPEND TPP_CXX_INCLUDE_PATH "$ENV{VTUNE_AMPLIFIER_XE_2018_DIR}/include" )

      list ( APPEND TPP_Fortran_LINK_LINE "$ENV{VTUNE_AMPLIFIER_XE_2018_DIR}/lib64 -littnotify" )
      list ( APPEND TPP_C_LINK_LINE "$ENV{VTUNE_AMPLIFIER_XE_2018_DIR}/lib64 -littnotify" )
      list ( APPEND TPP_CXX_LINK_LINE "$ENV{VTUNE_AMPLIFIER_XE_2018_DIR}/lib64 -littnotify" )
   endif ()

endif ()
