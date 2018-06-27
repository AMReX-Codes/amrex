#
#
# FUNCTION: configure_amrex
#
# Set all the properties (except sources and installation-related)
# required to build target "amrex".
# Target "amrex" must exist befroe this function is called.
#
# Author: Michele Rosso
# Date  : June 26, 2018
#
#
function (configure_amrex)

   # 
   # Check if target "amrex" has been defined before
   # calling this macro
   #
   if (NOT TARGET amrex)
      message (FATAL_ERROR "Target 'amrex' must be defined before calling function 'configure_amrex'" )     
   endif ()

   #
   # Check that needed options have already been defined
   # 
   if ( ( NOT ( DEFINED ENABLE_MPI ) ) OR ( NOT (DEFINED ENABLE_OMP) ) 
	 OR ( NOT (DEFINED ENABLE_PIC ) ) )
      message ( AUTHOR_WARNING "Required options are not defined" )
   endif ()

   #
   # Include the required modules
   # 
   include ( AMReX_ThirdPartyProfilers )
   include ( AMReX_Defines )
   include ( AMReX_Compilers )
   
   # 
   # Set compile definitions
   # 
   set_amrex_defines ()

   #
   # Set compiler flags
   #
   set_amrex_compilers ()

   #
   # Decide whether or not to use PIC 
   #
   if ( ENABLE_PIC OR BUILD_SHARED_LIBS )
      set (CMAKE_POSITION_INDEPENDENT_CODE TRUE)
   endif ()
   
   #
   # Include paths for Fortran modules
   #
   target_include_directories ( amrex
      PUBLIC "$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>" )
   
   #
   # Setup MPI (CMake 3.11 allows to import MPI as an imported TARGET)
   #  
   if (ENABLE_MPI)
      find_package (MPI REQUIRED)

      # Includes
      target_include_directories ( amrex PUBLIC
	 $<$<COMPILE_LANGUAGE:Fortran>:${MPI_Fortran_INCLUDE_PATH}>
	 $<$<COMPILE_LANGUAGE:C>:${MPI_C_INCLUDE_PATH}>
	 $<$<COMPILE_LANGUAGE:CXX>:${MPI_CXX_INCLUDE_PATH}>
	 )

      # Additional compiler flags
      target_compile_options ( amrex PUBLIC
	 $<$<COMPILE_LANGUAGE:Fortran>:${MPI_Fortran_COMPILE_FLAGS}>
	 $<$<COMPILE_LANGUAGE:C>:${MPI_C_COMPILE_FLAGS}>
	 $<$<COMPILE_LANGUAGE:CXX>:${MPI_CXX_COMPILE_FLAGS}>
	 )
      
      # Libraries to link to + linker flags
      # CMake doc suggests not to use absolute paths for portability reasons
      # (and we ignore it for now)
      if (MPI_CXX_LINK_FLAGS)
      	 string ( STRIP ${MPI_CXX_LINK_FLAGS} MPI_CXX_LINK_FLAGS )
      endif ()
      target_link_libraries ( amrex PUBLIC ${MPI_CXX_LINK_FLAGS}
	 ${MPI_Fortran_LIBRARIES} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
      
   endif ()

   #
   # Setup OpenMP
   # 
   if (ENABLE_OMP)
      find_package (OpenMP REQUIRED)

      # Additional compiler flags
      target_compile_options ( amrex PUBLIC
	 $<$<COMPILE_LANGUAGE:Fortran>:${OpenMP_Fortran_FLAGS}>
	 $<$<COMPILE_LANGUAGE:C>:${OpenMP_C_FLAGS}>
	 $<$<COMPILE_LANGUAGE:CXX>:${OpenMP_CXX_FLAGS}>
	 )

      # The openMP flags are required also during linking phase
      # The NON-clean way of achiving this is by adding a flag to the linker
      # as follows. Notice that in more recent version of CMake, openMP becomes
      # an imported target, thus provinding a clean way to do this
      if (OpenMP_CXX_FLAGS)
	 string ( STRIP "${OpenMP_CXX_FLAGS}" OpenMP_CXX_FLAGS )
      endif()    
      target_link_libraries ( amrex PUBLIC ${OpenMP_CXX_FLAGS} )

   else ()
      # Cray compiler has OMP turned on by default
      target_compile_options ( amrex PUBLIC
	 "$<$<C_COMPILER_ID:Cray>:-h noomp>" )
   endif()

   #
   # Setup third-party profilers
   # 
   set_amrex_profilers ()
   
   #
   # Print out summary
   # 
   print_amrex_configuration_summary ()
   
endfunction ()

# 
#
# Prints out configuration details
# 
# 
function (print_amrex_configuration_summary)

   # 
   # Check if target "amrex" has been defined before
   # calling this macro
   #
   if (NOT TARGET amrex)
      message (FATAL_ERROR "Target 'amrex' must be defined before calling function 'set_amrex_defines'" )
   endif ()

   
   #
   # Retrieve defines
   #
   get_target_property ( AMREX_DEFINES amrex COMPILE_DEFINITIONS )
   extract_by_language ( AMREX_DEFINES AMREX_GENERAL_DEFINES AMREX_Fortran_DEFINES )
   string (REPLACE " " " -D" AMREX_GENERAL_DEFINES "-D${AMREX_GENERAL_DEFINES}" )
   string (REPLACE ";" " -D" AMREX_GENERAL_DEFINES "${AMREX_GENERAL_DEFINES}" )
   string (REPLACE " " " -D" AMREX_Fortran_DEFINES "-D${AMREX_Fortran_DEFINES}" )
   string (REPLACE ";" " -D" AMREX_Fortran_DEFINES "${AMREX_Fortran_DEFINES}" )
   
   #
   # Retrieve compile flags
   #
   get_target_property ( AMREX_FLAGS   amrex COMPILE_OPTIONS)
   extract_by_language ( AMREX_FLAGS AMREX_CXX_FLAGS AMREX_Fortran_FLAGS)
   string ( REPLACE ";" " " AMREX_CXX_FLAGS "${AMREX_CXX_FLAGS}" )
   string ( REPLACE ";" " " AMREX_Fortran_FLAGS "${AMREX_Fortran_FLAGS}" )

   #
   # Include paths
   #
   get_target_property ( AMREX_INCLUDE_PATHS amrex INTERFACE_INCLUDE_DIRECTORIES )
   extract_by_language ( AMREX_INCLUDE_PATHS AMREX_CXX_INCLUDE_PATHS AMREX_Fortran_INCLUDE_PATHS)

   #
   # Link libraries
   # 
   get_target_property ( AMREX_LINK_LINE amrex LINK_LIBRARIES )
   if (NOT AMREX_LINK_LINE) # LINK_LIBRARIES property can return "NOT_FOUND"
      set (AMREX_LINK_LINE "")
   endif ()   
   list ( REMOVE_DUPLICATES AMREX_LINK_LINE )
   string ( REPLACE ";" " " AMREX_LINK_LINE "${AMREX_LINK_LINE}" )
   
   
   #
   # Config summary
   #
   message( STATUS "AMReX configuration summary: ")
   message( STATUS "   Build type               = ${CMAKE_BUILD_TYPE}")
   message( STATUS "   Install directory        = ${CMAKE_INSTALL_PREFIX}")
   message( STATUS "   General defines          = ${AMREX_GENERAL_DEFINES}")
   message( STATUS "   Fortran-specific defines = ${AMREX_Fortran_DEFINES}")
   message( STATUS "   C++ compiler             = ${CMAKE_CXX_COMPILER}")
   message( STATUS "   Fortran compiler         = ${CMAKE_Fortran_COMPILER}")
   message( STATUS "   C++ flags                = ${AMREX_CXX_FLAGS}")
   message( STATUS "   Fortran flags            = ${AMREX_Fortran_FLAGS}")
   message( STATUS "   C++ include paths        = ${AMREX_CXX_INCLUDE_PATHS}")  
   message( STATUS "   Fortran include paths    = ${AMREX_Fortran_INCLUDE_PATHS}")
   message( STATUS "   Link line                = ${AMREX_LINK_LINE}") 

endfunction ()


#
#  Helper macro to differentiate configuration options
#  by language
# 
macro ( extract_by_language list cxx_list fortran_list )

   #
   # Re-position list separators to reshape list 
   # 

   # Replace all ; with a place holder (*)
   string ( REPLACE ";" "*" ${list} "${${list}}" )

   # Add list delimiter only where it suits us
   string ( REPLACE ">*" ">;" ${list} "${${list}}" )
   string ( REPLACE "*$" ";$" ${list} "${${list}}" )
   string ( REPLACE "*/" ";/" ${list} "${${list}}" )
   string ( REPLACE "*" " "   ${list} "${${list}}" )
   
   #
   # First remove entries related to a compiler other than
   # the one currently in use or BUILD_INTERFACES keyword
   # 
   set (tmp_list)
   foreach ( item IN ITEMS ${${list}} )
      string (REPLACE "$<" "" item ${item} )
      string (REPLACE ">" "" item ${item} )
      string (REPLACE ":" "" item ${item} )

      # Skip build interface generator expressions 
      string ( FIND ${item} "BUILD_INTERFACE" idx )
      if ( ${idx} GREATER -1 )
	 continue ()
      endif()
      
      string ( FIND ${item} "C_COMPILER_ID" idx1 )
      if ( ${idx1} GREATER -1 )
	 string ( FIND ${item} "${CMAKE_C_COMPILER_ID}" idx2 )
	 if ( ${idx2} GREATER -1 )
	    string (REPLACE "C_COMPILER_ID" "" item ${item} )
	    string (REPLACE "${CMAKE_C_COMPILER_ID}" "" item ${item} )
	    list (APPEND tmp_list ${item})
	 endif ()
      else ()
	 list (APPEND tmp_list ${item})
      endif ()
   endforeach ()
   
   set (${fortran_list})
   set (${cxx_list})
   
   foreach ( item IN ITEMS ${tmp_list} )

      string ( FIND ${item} "COMPILE_LANGUAGEFortran" idx1 )
      string ( FIND ${item} "COMPILE_LANGUAGECXX"     idx2 )
      string ( FIND ${item} "COMPILE_LANGUAGEC"       idx3 )
      
      if ( ${idx1} GREATER -1 )
	 string (REPLACE "COMPILE_LANGUAGEFortran" "" item ${item} )
	 if ( NOT ( item STREQUAL "") )
   	    list ( APPEND ${fortran_list} ${item} )
	 endif ()
      elseif ( ${idx2} GREATER -1)
	 string (REPLACE "COMPILE_LANGUAGECXX" "" item ${item} )
	 if ( NOT ( item STREQUAL "") )
   	    list ( APPEND ${cxx_list} ${item} )
	 endif ()
      elseif ( ${idx3} GREATER -1 )
      else ()
	 list ( APPEND ${cxx_list} ${item} )	    
      endif () 
     
   endforeach ()

   if (${cxx_list})
      list ( REMOVE_DUPLICATES ${cxx_list} )
   endif ()
   if (${fortran_list})
      list ( REMOVE_DUPLICATES ${fortran_list} )
   endif ()
   
   unset (tmp_list)
   
endmacro ()
