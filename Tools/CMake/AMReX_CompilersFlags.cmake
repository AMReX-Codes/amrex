#
# This module provides the following macros
#
#    load_default_cxxflags: set CMAKE_CXX_FLAGS to AMReX defaults
#
#    load_default_fflags:   set CMAKE_Fortran_FLAGS to AMReX defaults
#



#
# Set CMAKE_CXX_FLAGS to AMReX defaults
# 
macro ( load_default_cxxflags )

   if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" )  # GNU compiler specific flags

      if ( DEBUG )
	 set ( CMAKE_CXX_FLAGS "-O0 -fno-inline -ggdb -Wall -Wno-sign-compare")
      else ()
	 set ( CMAKE_CXX_FLAGS "" )
      endif ()
      
   elseif ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel" )

      if ( DEBUG )
	 set ( CMAKE_CXX_FLAGS  "-O0 -traceback -Wcheck" )
      else ()
	 set ( CMAKE_CXX_FLAGS  "-ip -qopt-report=5 -qopt-report-phase=vec" )
      endif ()
      
   elseif ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "PGI")

      if ( DEBUG )
	 set ( CMAKE_CXX_FLAGS "-O0 -Mbounds")
      else ()
	 set ( CMAKE_CXX_FLAGS "-gopt -fast")
      endif ()
      
   elseif ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Cray" )

      if ( DEBUG )
	 set ( CMAKE_CXX_FLAGS "-O0")
      else ()
	 set ( CMAKE_CXX_FLAGS "")
      endif ()
      
   endif ()

endmacro ()

#
# Set CMAKE_Fortran_FLAGS to AMReX defaults
# 
macro ( load_default_fflags )

   if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" )  # GNU compiler specific flags

      if ( DEBUG )
	 set ( CMAKE_Fortran_FLAGS "-O0 -ggdb -fbounds-check -fbacktrace -Wuninitialized" )
	 set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wunused -finit-real=snan" )
	 set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -finit-integer=2147483647" )
      else ()
	 set ( CMAKE_Fortran_FLAGS "" )
      endif ()
      
   elseif ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel" )

      if ( DEBUG )
	 set ( CMAKE_Fortran_FLAGS  "-O0 -traceback -check bounds,uninit,pointers" )
      else ()
	 set ( CMAKE_Fortran_FLAGS  "-ip -qopt-report=5 -qopt-report-phase=vec" )
      endif ()
      
   elseif ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "PGI")

      if ( DEBUG )
	 set ( CMAKE_Fortran_FLAGS "-O0 -Mbounds -Mchkptr")
      else ()
	 set ( CMAKE_Fortran_FLAGS "-gopt -fast")
      endif ()
      
   elseif ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Cray" )

      if ( DEBUG )
	 set ( CMAKE_Fortran_FLAGS "-O0 -e i")
      else ()
	 set ( CMAKE_Fortran_FLAGS "")
      endif ()
      
   endif ()

endmacro ()


#
# Append required flags to CMAKE_CXX_FLAGS
# Required flags include those related to FPE
# 
macro ( append_required_cxxflags ) 

   if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" )  

      if ( ENABLE_FPE )
	 set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftrapv " )
      endif ()
      
   elseif ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel" )

      set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
      # if ( ENABLE_FPE )
      # 	 set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" )
      # endif ()
      
   elseif ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "PGI")

      # if ( ENABLE_FPE )
      # 	 set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" )
      # endif ()
      
   elseif ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Cray" )

      set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -h std=c++11 -h list=a")
      # if ( ENABLE_FPE )
      # 	 set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" )
      # endif ()
      
   endif ()

endmacro ()


#
# Append required flags to CMAKE_Fortran_FLAGS
# 
macro ( append_required_fflags ) 

   if ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" )
      
      set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-line-length-none" )
      set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none" )
      set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-range-check -fno-second-underscore")

      if ( ENABLE_FPE )
	 set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffpe-trap=invalid,zero -ftrapv" )
      endif ()
      
   elseif ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" )
      
      set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -extend_source")
      # if ( ENABLE_FPE )
      # 	 set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" )
      # endif ()

   elseif ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
      
      set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -extend")
      if ( ENABLE_FPE )
      	 set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ktrap=divz,inv" )
      endif ()
      
   elseif ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray" )

      set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -N 255 -h list=a")
      # if ( ENABLE_FPE )
      # 	 set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" )
      # endif ()

   endif ()

   # Add required Fortran-only definition
   set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${AMREX_Fortran_DEFINITIONS}" )
   
endmacro ()









