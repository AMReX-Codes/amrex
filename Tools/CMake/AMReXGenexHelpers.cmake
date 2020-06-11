function ( eval_conditional_expressions _out _in)

   # Replace each ";" with a place holder so that we can use REGEX MATCHALL
   # which returns a list of matches without getting messed up by
   # ";" which are not delimiters in a list
   # This is need it here because the genex below may contain ";" as part of
   # their argument or because CMake introduced it
   set(_semicol "$SC$")
   string(REPLACE ";" "${_semicol}" _in "${_in}")

   #
   # Genex in the form $<bool:true_string> where bool=<0|1>
   #
   string(REGEX REPLACE "\\$<0:[^>]*>" "" _in "${_in}")     # false
   string(REGEX MATCHALL "\\$<1:[^>]*>" _matches "${_in}")  # true
   foreach(_match IN LISTS _matches)
      string(REGEX REPLACE "(\\$<1:)|>$" "" _true_string "${_match}")
      string(REPLACE "${_match}" "${_true_string}" _in "${_in}")
   endforeach ()

   #
   # Genex in the form $<IF:condition,true_string,false_string> where bool=<0|1>
   #
   string(REGEX REPLACE "\\$<IF:0[^>]*>"  "" _in "${_in}")       # false
   string(REGEX MATCHALL "\\$<IF:1[^>]*>" _matches "${_in}")  # true
   foreach(_match IN LISTS _matches)
      string(REGEX MATCH ",.*," _true_string "${_match}")
      string(REPLACE "," "" _true_string "${_true_string}")
      string(REPLACE "${_match}" "${_true_string}" _in "${_in}")
   endforeach()

   # Put back the ";"s
   string(REPLACE "${_semicol}" ";" _in "${_in}")
   set(${_out} "${_in}" PARENT_SCOPE)

endfunction ()


#
# Does not support yet BOOL
#
function (eval_logical_operators _out _in)

   #
   # Genex in the form $<AND:conditions> where conditions=<0|1>
   #
   string(REGEX REPLACE "\\$<AND:[1,]*>" "1" _in "${_in}")  # Deal with "true" first
   string(REGEX REPLACE "\\$<AND:[^>]*>" "0" _in "${_in}")  # The remaining are "false"

   #
   # Genex in the form $<OR:conditions> where conditions=<0|1>
   #
   string(REGEX REPLACE "\\$<OR:[0,]*>" "0" _in "${_in}")  # Deal with "false" first
   string(REGEX REPLACE "\\$<OR:[^>]*>" "1" _in "${_in}")  # The remaining are "true"

   #
   # Genex in the form $<NOT:condition> where condition=<0|1>
   #
   string(REGEX REPLACE "\\$<NOT:1>" "0" _in "${_in}")
   string(REGEX REPLACE "\\$<NOT:0>" "1" _in "${_in}")

   set(${_out} "${_in}" PARENT_SCOPE)

endfunction ()


#
# Does not support yet STREQUAL, EQUAL, IN_LIST
#
function (eval_string_comparisons _out _in)

   #
   # Genex in the form $<ops:v1,v2> where ops=VERSION_* and v1,v2 = version strings
   #
   string(REGEX MATCHALL "\\$<VERSION[^>]*>" _matches "${_in}")

   foreach(_match IN LISTS _matches)
      string(REGEX MATCHALL "[0-9][^,>]*" _versions "${_match}")
      string(REGEX MATCHALL "[A-Z][^\$<:,>]*" _ops "${_match}")
      list(GET _versions 0 _v1)
      list(GET _versions 1 _v2)
      if (_v1 ${_ops} _v2)
         string(REPLACE "${_match}" "1" _in "${_in}")
      else ()
         string(REPLACE "${_match}" "0" _in "${_in}")
      endif ()
   endforeach()

   set(${_out} "${_in}" PARENT_SCOPE)

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
         foreach(_lang C CXX Fortran CUDA)
            if ("${_lang}" STREQUAL "${ARG_LANG}")
               string(REPLACE "$<COMPILE_LANGUAGE:${_lang}>" "${true}"
                  genex_list "${genex_list}")
            else ()
               string(REPLACE "$<COMPILE_LANGUAGE:${_lang}>" "${false}"
                  genex_list "${genex_list}")
            endif ()
         endforeach()
      endif ()

      string(REPLACE "$<COMPILE_LANGUAGE:" "${false}"
         genex_list "${genex_list}")

      # Tag by compiler
      if (ARG_COMP)
         string(REPLACE "$<C_COMPILER_ID:${ARG_COMP}>" "${true}"
            genex_list "${genex_list}")
         string(REPLACE "$<CXX_COMPILER_ID:${ARG_COMP}>" "${true}"
            genex_list "${genex_list}")
         string(REPLACE "$<Fortran_COMPILER_ID:${ARG_COMP}>" "${true}"
            genex_list "${genex_list}")
      endif ()
      string(REPLACE "$<C_COMPILER_ID:" "${false}" genex_list "${genex_list}")
      string(REPLACE "$<CXX_COMPILER_ID:" "${false}" genex_list "${genex_list}")
      string(REPLACE "$<Fortran_COMPILER_ID:" "${false}" genex_list "${genex_list}")


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

      # Take care of logical NOT
      string(REPLACE "NOT:${true}"  "${false}"  genex_list "${genex_list}")
      string(REPLACE "NOT:${false}" "${true}"   genex_list "${genex_list}")

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


#  _in must be a string
function ( eval_genex _out _in )

   set(_option_arg STRING)
   set(_one_value_args LANG COMP COMP_VERSION CONFIG INTERFACE)
   cmake_parse_arguments( ARG "${_option_arg}" "${_one_value_args}" "" ${ARGN} )

   # Check if there are any genex-es in the input. If not, return _in as is
   # and exit right away
   string(GENEX_STRIP "${_in}" _nogenex_string)
   if (_nogenex_string STREQUAL _in)
      set(${_out} "${_in}" PARENT_SCOPE)
      return()
   endif ()


   #
   #  Replace GENEX variable queries
   #

   # Genex in the form $<CONFIG:cfg>
   # Ignore if no CONFIG arg is given
   if (ARG_CONFIG)
      string(REGEX REPLACE "\\$<CONFIG:${ARG_CONFIG}>" "1"  _in "${_in}")
   endif ()
   string(REGEX REPLACE "\\$<CONFIG:[A-Za-z]*>" "0"  _in "${_in}")

   # Genex in the form $<*_COMPILER_ID:compiler_id>
   # Ignore if no COMP arg is given
   if (ARG_COMP)
      string(REGEX REPLACE "\\$<[A-Za-z]*_COMPILER_ID:${ARG_COMP}>" "1"  _in "${_in}")
   endif ()
   string(REGEX REPLACE "\\$<[A-Za-z]*_COMPILER_ID:[A-Za-z]*>" "0"  _in "${_in}")

   # Genex in the form $<*_COMPILER_VERSION:version>
   # Ignore if no COMP_VERSION arg is given
   if (ARG_COMP_VERSION)
      string(REGEX REPLACE "\\$<[A-Za-z]*_COMPILER_VERSION:${ARG_COMP_VERSION}>" "1"  _in "${_in}")
   endif ()
   string(REGEX REPLACE "\\$<[A-Za-z]*_COMPILER_ID:[^>]*>" "0"  _in "${_in}")

   # Genex in the form $<COMPILE_LANGUAGE:language>
   # Ignore if no LANG arg is given
   if (ARG_LANG)
      string(REGEX REPLACE "\\$<COMPILE_LANGUAGE:${ARG_LANG}>" "1"  _in "${_in}")
   endif ()
   string(REGEX REPLACE "\\$<COMPILE_LANGUAGE:[A-Za-z]*>" "0"  _in "${_in}")

   # Logical operators
   eval_logical_operators(_in "${_in}")

   # Conditional expressions
   eval_conditional_expressions(_in "${_in}")

   # Remove empty elements, i.e. the ";;" resulting from
   # deleting false genex.
   # In this case, _in is treated as a list
   list(REMOVE_ITEM _in "")

   # More clean up
   list(TRANSFORM _in STRIP)

   # If string is required
   if (ARG_STRING)
      string(REPLACE ";" " " _in "${_in}" )
      string(STRIP "${tmp}" tmp)
   endif ()

   set(${_out} "${_in}" PARENT_SCOPE)

endfunction ()
