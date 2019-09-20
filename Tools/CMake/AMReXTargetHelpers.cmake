#[=======================================================================[
AMReXTargetHelpers
----------------

Provides tools to perform non-standard operations on CMake targets
#]=======================================================================]

#
# Recursively finds the properties associated to _target via all its dependecies.
# Basically, it returns all ( including transitive ) includes, defines, flags,
# and link libraries associated with _target 
# This is essentially a wrapper for get_target_prop_recursive below and should
# be used in its stead. 
# 
function (get_target_properties_flattened _target _includes _defines _flags _link_line )

   if (NOT TARGET ${_target})
      message(FATAL_ERROR "Target ${_target} does not exist!")
   endif ()

   # Init local variables
   set(_local_includes)
   set(_local_defines)
   set(_local_flags)
   set(_local_link_line)

   get_target_prop_recursive(${_target} _local_includes _local_defines _local_flags _local_link_line)
   # Set results in parent scope
   set(${_includes}  "${_local_includes}"   PARENT_SCOPE)
   set(${_defines}   "${_local_defines}"    PARENT_SCOPE)
   set(${_flags}     "${_local_flags}"      PARENT_SCOPE)
   set(${_link_line} "${_local_link_line}"  PARENT_SCOPE)
      
endfunction ()


#
# Recursively finds the properties associated to _target via all its dependecies
# NOT RECOMMENDED to use this function directly.
# Use get_target_properties_flattened INSTEAD!
# 
function (get_target_prop_recursive _target _lincludes _ldefines _lflags _llink_line )

   if (NOT TARGET ${_target})
      return()
   endif ()

   # defines
   get_target_property(_interface_defines ${_target} INTERFACE_COMPILE_DEFINITIONS)
   if (_interface_defines)
      list(APPEND ${_ldefines} ${_interface_defines})
   endif ()
   
   # includes
   get_target_property(_interface_includes ${_target} INTERFACE_INCLUDE_DIRECTORIES)
   if (_interface_includes)
      list(APPEND ${_lincludes} ${_interface_includes})
   endif ()
   
   # compile options
   get_target_property(_interface_flags ${_target} INTERFACE_COMPILE_OPTIONS)
   if (_interface_flags)
      list(APPEND ${_lflags} ${_interface_flags})
   endif ()

   # libraries
   get_target_property(_interface_link_libraries ${_target} INTERFACE_LINK_LIBRARIES)

   # INTERFACE_LINK_LIBRARIES may contain targets in the form target-name::@<hexadecimal-code>
   # Remove the @<hexadecimal-code> bit to get the target name.
   # check CMake doc regarding INTERFACE_LINK_LIBRARIES for an in-depth explanation of this
   string(REGEX REPLACE ":*@+<[A-Za-z0-9]+>" "" _interface_link_libraries
      "${_interface_link_libraries}")
   
   # Remove INTERFACE genex: choose build
   include(AMReXGenexHelpers)
   evaluate_genex(_interface_link_libraries
      _interface_link_libraries
      CONFIG ${CMAKE_BUILD_TYPE}
      INTERFACE BUILD)
   
   if (_interface_link_libraries)
      foreach(_item IN LISTS _interface_link_libraries )
         if ( NOT TARGET ${_item} )
            list(APPEND ${_llink_line} ${_item})
         else ()
            get_target_prop_recursive(${_item} ${_lincludes} ${_ldefines} ${_lflags} ${_llink_line})
         endif ()
      endforeach ()
   endif ()

   # Set output variables to parent scope
   foreach (_var _lincludes _ldefines _lflags _llink_line)
      set(${${_var}} "${${${_var}}}" PARENT_SCOPE)
   endforeach ()
   
endfunction ()
