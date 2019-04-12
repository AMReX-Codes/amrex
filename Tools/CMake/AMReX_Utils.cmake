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
# WARNING: it is dangerous and definitely discouraged to use
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
# Strip string from trailing and leading whitespace
# after veryfing it is not empty
#
macro (strip var)
  if (${var})
    string ( STRIP "${${var}}" ${var} )
  endif ()  
endmacro ()


# 
#
# FUNCTION: add_amrex_define
#
# Add definition to target "amrex" compile definitions.
#
# Arguments:
#
#    new_define    = variable containing the definition to add.
#                    The new define can be any string with no "-D" prepended.
#                    If the new define is in the form AMREX_SOMENAME,
#                    this function also adds BL_SOMENAME to the list,
#                    unless NO_LEGACY is specified (see below)
# 
#    NO_LEGACY     = if specified, the legacy version of a new_define given in the
#                    form AMREX_SOMENAME will not be added.
#                     
#    IF <cond-var> = new_define is added only if <cond-var> is true
#
# Author: Michele Rosso
# Date  : June 26, 2018
# 
function ( add_amrex_define new_define )

   # 
   # Check if target "amrex" has been defined before
   # calling this macro
   #
   if (NOT TARGET amrex)
      message(FATAL_ERROR "Target 'amrex' must be defined before calling function 'add_amrex_define'" )
   endif ()
   
   cmake_parse_arguments( DEFINE "NO_LEGACY" "IF" ""  ${ARGN} )

   set( condition  1 )

   if (DEFINE_IF)
      set( condition ${${DEFINE_IF}} )
   endif()

   # Return if flags does not need to be included
   if (NOT condition)
      return()
   endif ()

   target_compile_definitions( amrex PUBLIC $<BUILD_INTERFACE:${new_define}> )

   if ( NOT DEFINE_NO_LEGACY )
      # Add legacy definition
      string( FIND ${new_define} "AMREX_" out )

      if (${out} GREATER -1 )
	 string(REPLACE "AMREX_" "BL_" legacy_define ${new_define})
	 target_compile_definitions( amrex PUBLIC $<BUILD_INTERFACE:${legacy_define}> )
      endif ()
   endif () 
   
endfunction ()



#
#
# FUNCTION: evaluate_genex
#
# Takes a list as input and evaluate the generator expressions
# according to user-defined criteria (see arguments below).
#
# Optional arguments:
#   LANG       <compile-language>
#   COMP       <compiler-id>
#   CONFIG     <config-type>
#   INTERFACE  <interface-type>  (BUILD is default)
#
# Additionally, passing the keyword STRING will return a string
# rather than a list (default).
# 
# Example:
#
# evaluate_genex(IN OUT COMP GNU LANG CXX STRING)
#
# Returns in OUT all the elements in IN that do not contain genex-es
# + all the elements in IN contaning genex-es which test for true
# when $<COMPILE_LANGUAGE:CXX> and the compiler id is GNU.
# Also, OUT will be returned as a string, i.e. without ";"
#
# Author: Michele Rosso
# Date  : Jan 30, 2019
#
#
function( evaluate_genex input output )

   set(option_arg STRING)
   set(one_value_args LANG COMP CONFIG INTERFACE)
   cmake_parse_arguments( ARG "${option_arg}" "${one_value_args}" "" ${ARGN} )

   
   # Split input into two lists: one containing genex, the other with
   # no genex inside.
   # Since ";" can be used as both item separator in a list and
   # semi-column, we replace ";" with placeholders to indicate which
   # of these two function ";" is performing.
   set(list_sep "$SEP$" ) # Use this to indicate a item separator in list
   set(semicol  "$SC$"  ) # Use this to indicate a simple semi-column 
   set(genex_list ${${input}})
   string(GENEX_STRIP "${genex_list}" nogenex_list)
   if (nogenex_list)
      list(REMOVE_ITEM genex_list ${nogenex_list})
   endif ()

   if (genex_list)
      # The genex_list is a list of genex-es in the form
      # $<SOME_GENEX:something;something;something;...>.
      # In order to avoid genex-ex containg ";" to be split over multiple
      # items in the list - ";" is also the list separator - we replace
      # all the ";" meant to be regular semicolumns with a placeholder
      # to be converted back to ";" later on.
      string(REPLACE ">;" "${list_sep}" genex_list "${genex_list}")
      string(REPLACE ";" "${semicol}" genex_list "${genex_list}" )
      string(REPLACE  "${list_sep}" ";" genex_list "${genex_list}" )

      # Now we can tag wih true or false the genex-es which satisfy
      # or not our criteria  
      set(true ".true.")
      set(false ".false.")

      # Tag by language
      if (ARG_LANG)
         string(REPLACE "$<COMPILE_LANGUAGE:${ARG_LANG}>" "${true}"
            genex_list "${genex_list}")
      endif ()
      string(REPLACE "$<COMPILE_LANGUAGE:" "${false}"
         genex_list "${genex_list}")

      # Tag by compiler
      if (ARG_COMP)
         string(REPLACE "$<C_COMPILER_ID:${ARG_COMP}>" "${true}"
            genex_list "${genex_list}")
         string(REPLACE "$<CXX_COMPILER_ID:${ARG_COMP}>" "${true}"
            genex_list "${genex_list}")
         string(REPLACE
            "$<STREQUAL:\"${CMAKE_Fortran_COMPILER_ID}\",\"${ARG_COMP}\">"
            "${true}" genex_list "${genex_list}")
      endif ()
      string(REPLACE "$<C_COMPILER_ID:" "${false}" genex_list "${genex_list}")
      string(REPLACE "$<CXX_COMPILER_ID:" "${false}" genex_list "${genex_list}")
      string(REPLACE "$<STREQUAL:"        "${false}" genex_list "${genex_list}")

      # Tag by configuration
      if (ARG_COMP)
         string(REPLACE "$<CONFIG:${ARG_CONFIG}>" "${true}"
            genex_list "${genex_list}")
      endif ()
      string(REPLACE "$<CONFIG:" "${false}" genex_list "${genex_list}")

      # Tag by interface
      if (ARG_INTERFACE)
         if ( "${ARG_INTERFACE}" STREQUAL "BUILD")
            string(REPLACE "$<BUILD_INTERFACE:" "${true}"
               genex_list "${genex_list}")
         elseif ("${ARG_INTERFACE}" STREQUAL "INSTALL")
            string(REPLACE "$<INSTALL_INTERFACE:" "${true}"
               genex_list "${genex_list}")         
         endif ()
      endif ()
      # Build interface is default
      string(REPLACE "$<BUILD_INTERFACE:" "${true}"  genex_list "${genex_list}")
      string(REPLACE "_INTERFACE:"        "${false}" genex_list "${genex_list}")   
      
      # Remove genex-es not satisfying our criteria
      set(valid_genex_list ${genex_list})
      if(valid_genex_list)
         list(FILTER valid_genex_list EXCLUDE REGEX "${false}")
      endif()

      # Replace back regular semi-columns
      string(REPLACE  "${semicol}" ";" valid_genex_list "${valid_genex_list}" )

      #
      # remove all remaining placeholders and genex characters
      # OR operation not supported for now
      #
      string(REPLACE "$<AND:"    "" valid_genex_list "${valid_genex_list}")
      # string(REPLACE "$<OR:"     "" valid_genex_list "${valid_genex_list}")
      string(REPLACE "${true},"  "" valid_genex_list "${valid_genex_list}")
      string(REPLACE "${true}>:" "" valid_genex_list "${valid_genex_list}")
      string(REPLACE "${true}"   "" valid_genex_list "${valid_genex_list}")
      string(REPLACE ">:"        "" valid_genex_list "${valid_genex_list}")
      string(REPLACE "<:"        "" valid_genex_list "${valid_genex_list}")
      string(REGEX REPLACE "\\$|<|>" "" valid_genex_list "${valid_genex_list}")
   else()
     set(valid_genex_list "") 
   endif ()

   if (ARG_STRING)
      string(REPLACE ";" " " tmp "${nogenex_list};${valid_genex_list}" )
      string(STRIP "${tmp}" tmp)
      set(${output} ${tmp} PARENT_SCOPE)
   else ()
      set(${output} ${nogenex_list} ${valid_genex_list} PARENT_SCOPE)
   endif ()
   
endfunction ()
