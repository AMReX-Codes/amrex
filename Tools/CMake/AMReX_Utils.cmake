###################################################################
# Check if dir or file given by path exists and issue a warning or
#  error if not
##################################################################
function ( check_path  path  message_type )
  if ( EXISTS ${path} )
  else ()
    message(${message_type} ${path} " does not exist!")
  endif ( EXISTS ${path} )
endfunction ()


#
# This function turns a list into a string
# 
function ( list_to_string list )
  string (REPLACE ";" " " tmp "${${list}}")
  set ( ${list} "${tmp}" PARENT_SCOPE)
endfunction ()

#
# Create list of all include directories
# cmake must be re-run if new dirs with Headers are introduced
#
# Arguments:
#
#  dirlist  = the list of subdir
#  ROOT     = top level directory from where to start search
#             If not given, default is CMAKE_CURRENT_LIST_DIR
#  EXCLUDE  = list of path to exclude from search
#
function ( find_include_paths dirlist )

  cmake_parse_arguments ( ARG "" "ROOT" "EXCLUDE"  ${ARGN} )

  if (NOT ARG_ROOT)
    set ( ARG_ROOT ${CMAKE_CURRENT_LIST_DIR} )
  endif ()

  # Check that root and exclude dir exist
  set ( alldirs ${ARG_ROOT} )
  if ( ARG_EXCLUDE )
    list (APPEND alldirs ${ARG_EXCLUDE} )
  endif ()

  foreach ( dir ${alldirs} )
    if ( NOT EXISTS ${dir} )
      message (WARNING "${dir} is not a valid path")
    endif ()
  endforeach ()

  # This list all the directories containing headers 
  file ( GLOB_RECURSE includes LIST_DIRECTORIES true 
    ${ARG_ROOT}/*.h ${ARG_ROOT}/*.H )
  
  
  foreach (item ${includes})

    get_filename_component ( path ${item} PATH )
    
    if (IS_DIRECTORY ${path})

      # Check first if it is a valid path
      set (path_is_valid "YES")
      
      foreach ( exclude ${ARG_EXCLUDE})
	string (FIND ${path} ${exclude} out )
	if ( NOT (${out} EQUAL -1) ) 
	  set (path_is_valid "NO")
	endif ()
      endforeach ()

      
      if ( NOT (${path} IN_LIST tmp ) AND path_is_valid )	   	    
	list ( APPEND tmp ${path} )
      endif ()
      
    endif ()
    
  endforeach ()

  
  
  set ( ${dirlist} ${tmp} PARENT_SCOPE )
  
endfunction ()


#
# Append new_var to all_var
# 
function ( append new_var all_var )
  if ( ${new_var} )
    set ( tmp  "${${all_var}} ${${new_var}}" )

    # Since I am OCD, remove the double spaces.
    string ( REPLACE "  " " " tmp ${tmp} )
    set ( ${all_var}  ${tmp} PARENT_SCOPE )
  endif ()
endfunction ()


#
# Print variable (useful for debugging)
#
function ( print var )
  message (" ${var} = ${${var}}" )
endfunction ()

#
# Print list
#
function ( print_list list )
  
  list ( LENGTH ${list} len )

  if ( ${len} GREATER 0 )
    message ("")
    message ( STATUS " LIST NAME:  ${list}")
    foreach ( item ${${list}})
      message ( STATUS "  ${item}")
    endforeach ()
    message ("")
  endif ()
  
endfunction ()


#
# Function to append to link line
#
function ( append_to_link_line libs link_line )

  if ( ${ARGC} EQUAL 3 )  # Third one is optional flags
    set ( flags  ${ARGV2} )
  else ()
    set ( flags )
  endif ()
  
  set ( tmp "${${link_line}} ${${flags}} ${${libs}} " )
  string ( STRIP "${tmp}" tmp )
  set ( ${link_line} ${tmp} PARENT_SCOPE )
  
endfunction ()


#
# Function to install include files
#
function ( install_include_files )
  foreach (file ${ARGV})
    install ( FILES ${file} DESTINATION include )
  endforeach()
endfunction ()

#
# Function to prepend path to list items
#
function (prepend list prefix)

  set ( tmp "" )
  foreach (item ${${list}})
    set ( name   ${prefix}/${item} )
    string ( REPLACE "//" "/" name ${name})
    list ( APPEND tmp ${name} )
  endforeach ()

  set ( ${list} ${tmp}  PARENT_SCOPE )

endfunction ()



#
#  USE AT YOUR OWN RISK
#
function (scan_for_sources f90src f77src cxxsrc allheaders)

  cmake_parse_arguments ( ARG "" "ROOT" ""  ${ARGN} )

  if (NOT (ARG_ROOT))
    set (ARG_ROOT ${CMAKE_CURRENT_LIST_DIR})
  endif ()
  
  file (GLOB_RECURSE tmp  "${ARG_ROOT}/*.f90"
    "${ARG_ROOT}/*.F90")
  set (${f90src} ${tmp} PARENT_SCOPE)

  
  file (GLOB_RECURSE f77src  "${ARG_ROOT}/*.f"
    "${ARG_ROOT}/*.F")

  file (GLOB_RECURSE cxxsrc  "${ARG_ROOT}/*.cpp" )

  file (GLOB_RECURSE allheaders  "${ARG_ROOT}/*.H"
    "${ARG_ROOT}/*.H")

  set (f90src)

endfunction ()


#
# Find all source files ( Fortran, C++, Headers )
# in CMAKE_CURRENT_LIST_DIR.
#
# Arguments:
#
#  source_list  = the list of sources (prefixed with their absolute path)
#  ROOT         = directory to search.
#                 If not given, default is CMAKE_CURRENT_LIST_DIR
#  RECURSE      = if given, enables search for subdirectories
#
# This macro returns a list of files with their absolute paths
#
# WARNING: it is dangerous and definetely discouraged to use
#          GLOBBING to find sources. Use at your own risk.
# 
macro ( find_all_sources sources_list include_paths )

  cmake_parse_arguments ( ARG "RECURSE" "ROOT" ""  ${ARGN} )

  if ( NOT (ARG_ROOT) )
    set (ARG_ROOT ${CMAKE_CURRENT_LIST_DIR})
  endif ()

  if (ARG_RECURSE)
    set ( search_type GLOB_RECURSE )
  else ()
    set ( search_type GLOB )
  endif ()
  
  file ( ${search_type} ${sources_list}
    "${ARG_ROOT}/*.f90"
    "${ARG_ROOT}/*.F90"
    "${ARG_ROOT}/*.F"
    "${ARG_ROOT}/*.cpp"
    "${ARG_ROOT}/*.H"
    "${ARG_ROOT}/*.h"
    )
  
  unset (search_type)

  # Find include paths
  if ( ${sources_list} )
    
    set ( include_paths )
    
    foreach ( file IN LISTS ${sources_list} )

      get_filename_component ( ext ${file} EXT )

      if ( ("${ext}" STREQUAL ".h") OR ("${ext}" STREQUAL ".H") )

	get_filename_component (path ${file} DIRECTORY)
	list ( APPEND ${include_paths} ${path} ) 

      endif()
      
    endforeach ()

    unset (path)

    if ( ${include_paths} )
      list ( REMOVE_DUPLICATES ${include_paths} )
    endif ()
    
  endif ()

endmacro ()



#
# This sets CMake_<LANG>_FLAGS_<CONFIG> to default values
# if the variable is empty
#
macro ( set_default_config_flags )
  
  if ( NOT CMAKE_Fortran_FLAGS_DEBUG )
    set (CMAKE_Fortran_FLAGS_DEBUG "-g")
  endif ()

  if ( NOT CMAKE_Fortran_FLAGS_RELEASE )
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  endif ()

  if ( NOT CMAKE_CXX_FLAGS_DEBUG )
    set (CMAKE_CXX_FLAGS_DEBUG "-g")
  endif ()
  
  if ( NOT CMAKE_CXX_FLAGS_RELEASE )
    set (CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
  endif ()
  
endmacro ()





#
#  Helper macro to replace genex in list 
# 
macro (replace_genex input_list output_list )

  cmake_parse_arguments ( ARG "" "LANGUAGE" "" ${ARGN} )

  set (tmp_list ${${input_list}})
  
  # Replace all ; with a place holder (*)
  string ( REPLACE ";" "*" tmp_list "${tmp_list}" )

  # Add tmp_list delimiter only where it suits us
  string ( REPLACE ">*" ">;" tmp_list "${tmp_list}" )
  string ( REPLACE "*$" ";$" tmp_list "${tmp_list}" )
  string ( REPLACE "*/" ";/" tmp_list "${tmp_list}" )
  string ( REPLACE "*" " "   tmp_list "${tmp_list}" )
  
  #
  # First remove entries related to:
  # 1) a compiler other than the one currently in use
  # 2) a build type other than the current one
  # 
  foreach ( item IN ITEMS ${tmp_list} )
    string (REPLACE "$<" "" item ${item} )
    string (REPLACE ">" "" item ${item} )
    string (REPLACE ":" "" item ${item} )

    # Accept build interface generator expressions 
    string (REPLACE "BUILD_INTERFACE" "" item ${item})

    # Skip genex for compilers other than the one in use
    string ( FIND ${item} "C_COMPILER_ID" idx1 )      
    if ( ${idx1} GREATER -1 )
      string ( FIND ${item} "C_COMPILER_ID${CMAKE_C_COMPILER_ID}" idx2 )
      if ( ${idx2} GREATER -1 )
   	string (REPLACE "C_COMPILER_ID${CMAKE_C_COMPILER_ID}" "" item ${item} )
      else ()
   	continue ()
      endif ()
    endif ()

    string ( FIND ${item} "STREQUAL" idx1 )      
    if ( ${idx1} GREATER -1 )
      string ( FIND ${item} "\"${CMAKE_Fortran_COMPILER_ID}\",\"${CMAKE_Fortran_COMPILER_ID}\"" idx2 )
      if ( ${idx2} GREATER -1 )
   	string (REPLACE "STREQUAL\"${CMAKE_Fortran_COMPILER_ID}\",\"${CMAKE_Fortran_COMPILER_ID}\"" "" item ${item} )
      else ()
   	continue ()
      endif ()
    endif ()

    string ( FIND ${item} "CXX_COMPILER_ID" idx1 )
    if ( ${idx1} GREATER -1 )
      string ( FIND ${item} "CXX_COMPILER_ID${CMAKE_CXX_COMPILER_ID}" idx2 )
      if ( ${idx2} GREATER -1 )
	string (REPLACE "CXX_COMPILER_ID${CMAKE_CXX_COMPILER_ID}" "" item ${item} )
      else ()
	continue ()
      endif ()
    endif ()
    
    string (FIND ${item} "CONFIG" idx3 )
    if ( ${idx3} GREATER -1 )
      string (FIND ${item} "${CMAKE_BUILD_TYPE}" idx4)
      if ( ${idx4} GREATER -1 )
   	string (REPLACE "CONFIG${CMAKE_BUILD_TYPE}" "" item ${item} )
      else ()
   	continue ()
      endif () 
    endif ()

    # Extract by Language part
    if ( ARG_LANGUAGE )
      string ( FIND ${item} "COMPILE_LANGUAGE" idx1 )
      if (${idx1} GREATER -1)
	if (${ARG_LANGUAGE} STREQUAL Fortran )
	  string ( FIND ${item} "Fortran" idx2 )
	  if ( ${idx2} GREATER -1)
	    string (REPLACE "COMPILE_LANGUAGEFortran" "" item ${item} )
	  else()
	    continue ()
	  endif ()
	elseif (${ARG_LANGUAGE} STREQUAL CXX)
	  string ( FIND ${item} "CXX" idx2 )
	  if ( ${idx2} GREATER -1)
	    string (REPLACE "COMPILE_LANGUAGECXX" "" item ${item} )
	  else()
	    continue ()
	  endif ()
	endif ()
      endif ()	    
    endif ()
    
    # Now item should be ok to be added to final list
    list ( APPEND ${output_list} ${item})
    
  endforeach ()

  if (${output_list})
    list (REMOVE_DUPLICATES ${output_list} )
  endif ()
  
endmacro ()


#
# Strip string from trailing and leading whitespace
# after veryfing it is not empty
#
macro (strip var)
  if (${var})
    string ( STRIP "${${var}}" ${var} )
  endif ()  
endmacro ()
