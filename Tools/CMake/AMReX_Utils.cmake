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
# Find Git Version
#

#
# Find Git Version
#
function (find_git_version version)

   set (output "")
   
   # Check whether .git is present and git installed
   if (EXISTS ${CMAKE_SOURCE_DIR}/.git AND ${GIT_FOUND})
      execute_process ( COMMAND git describe --abbrev=12 --dirty --always --tags
	 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
	 OUTPUT_VARIABLE output )
      string (STRIP ${output} output)
   endif ()

   set ( ${version} ${output} PARENT_SCOPE )

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
