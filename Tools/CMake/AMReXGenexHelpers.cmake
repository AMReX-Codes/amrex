function ( eval_conditional_expressions _out _in)

   # IMPORTANT: the group "\\", i.e. quoted double backslash, is parsed as a single
   # backslash. Therefore, "\\$" is parsed as just \$.

   #
   # Genex in the form $<condition:true_string> where condition=<0|1>
   #
   string(REGEX REPLACE "\\$<1:([^\$>]*)>" "\\1" _in "${_in}")  # true (\\1 is first subexpression in the match)
   string(REGEX REPLACE "\\$<0:[^\$>]*>"   ""    _in "${_in}")  # false

   #
   # Genex in the form $<IF:condition,true_string,false_string> where condition=<0|1>
   #
   string(REGEX REPLACE "\\$<IF:0,[^\$,]*,([^\$>]*)>"   "\\1" _in "${_in}")  # false
   string(REGEX REPLACE "(\\$<IF):1,([^\$,]*),[^\$>]*>" "\\2" _in "${_in}")  # true

   set(${_out}  "${_in}" PARENT_SCOPE)

endfunction ()


#
#
#
function (eval_logical_operators _out _in)

   #
   # Genex in the form $<BOOL:string>
   #
   set(_semicol "$SC$")   # Use to replace ; with a custom mark so the following foreach doesn't get confused
   string(REPLACE ";" "${_semicol}" _in "${_in}")
   string(REGEX MATCHALL "\\$<BOOL:[^\$>]*>" _matches "${_in}" )
   foreach (_match IN LISTS _matches)
      string(REGEX MATCHALL "\\$<BOOL:([^>]*)>" _tmp "${_match}" )
      if (CMAKE_MATCH_1)
         string(REPLACE "${_match}" "1" _in "${_in}")
      else()
         string(REPLACE "${_match}" "0" _in "${_in}")
      endif ()
   endforeach ()
   string(REPLACE "${_semicol}" ";" _in "${_in}")

   #
   # Genex in the form $<AND:conditions> where conditions=<0|1>
   #
   string(REGEX REPLACE "\\$<AND:[1,]*>"           "1" _in "${_in}")  # True
   string(REGEX REPLACE "\\$<AND:[01\,]*0[01\,]*>" "0" _in "${_in}")  # False

   #
   # Genex in the form $<OR:conditions> where conditions=<0|1>
   #
   string(REGEX REPLACE "\\$<OR:[0,]*>"           "0" _in "${_in}")  # True
   string(REGEX REPLACE "\\$<OR:[01\,]*1[01\,]*>" "1" _in "${_in}")  # False

   #
   # Genex in the form $<NOT:condition> where condition=<0|1>
   #
   string(REGEX REPLACE "\\$<NOT:1>" "0" _in "${_in}")
   string(REGEX REPLACE "\\$<NOT:0>" "1" _in "${_in}")

   set(${_out} "${_in}" PARENT_SCOPE)

endfunction ()


#
#
#
function (eval_string_comparisons _out _in)

   # To avoid confusing ";" with list separators when using REGEX MATCHALL
   set(_semicol "$SC$")
   string(REPLACE ";" "${_semicol}" _in "${_in}")

   #
   # Genex in the form $<ops:a1,a2> where ops=VERSION_*,EQUAL,STREQUAL,IN_LIST
   #
   string(REGEX MATCHALL "\\$<(VERSION|EQUAL|STREQUAL|IN_LIST)[^\$>]*>" _matches "${_in}")
   foreach(_match IN LISTS _matches)
      string(REGEX MATCH "([^\$<:]*):([^,]*),([^>]*)" _tmp "${_match}")
      set(_ops ${CMAKE_MATCH_1})
      set(_a1  ${CMAKE_MATCH_2})
      set(_a2  ${CMAKE_MATCH_3})
      if (_a1 ${_ops} _a2)
         string(REPLACE "${_match}" "1" _in "${_in}")
      else ()
         string(REPLACE "${_match}" "0" _in "${_in}")
      endif ()
   endforeach()

   string(REPLACE "${_semicol}" ";" _in "${_in}")
   set(${_out} "${_in}" PARENT_SCOPE)

endfunction ()

#
#
#
function (eval_string_transformations _out _in)

   # To avoid confusing ";" with list separators when using REGEX MATCHALL
   set(_semicol "$SC$")
   string(REPLACE ";" "${_semicol}" _in "${_in}")

   # Genex that deal with strings
   string(REGEX MATCHALL "\\$<(LOWER_CASE|UPPER_CASE)[^\$>]*>" _matches "${_in}")
   foreach(_match IN LISTS _matches)
      string(REGEX MATCH "([^\$<:]*):([^>]*)" _tmp "${_match}")
      set(_ops  ${CMAKE_MATCH_1})
      set(_args  ${CMAKE_MATCH_2})
      string(REPLACE "UPPER_CASE" "TOUPPER" _ops ${_ops})
      string(REPLACE "LOWER_CASE" "TOLOWER" _ops ${_ops})
      string(${_ops} ${CMAKE_MATCH_2} CMAKE_MATCH_2)
      string(REPLACE "${_match}" "${CMAKE_MATCH_2}" _in "${_in}")
   endforeach()

   # Genex that deal with lists
   string(REGEX MATCHALL "\\$<(JOIN|REMOVE_DUPLICATES|FILTER)[^\$>]*>" _matches "${_in}")
   foreach(_match IN LISTS _matches)
      string(REGEX MATCH "([^\$<:]*):([^>]*)" _tmp "${_match}")
      set(_ops  ${CMAKE_MATCH_1})
      set(_args  ${CMAKE_MATCH_2})
      string(REPLACE "$SC$" ";" _args "${_args}")

      if (_ops MATCHES "REMOVE_DUPLICATES")
         list(REMOVE_DUPLICATES _args)
         string(REPLACE "${_match}" "${_args}" _in "${_in}")
      elseif (_ops MATCHES "JOIN")
         string(REGEX MATCH "([^,]*),(.*)" _tmp "${_args}")
         list(JOIN CMAKE_MATCH_1 ${CMAKE_MATCH_2} _args)
         string(REPLACE "${_match}" "${_args}" _in "${_in}")
      elseif (_ops MATCHES "FILTER")
         string(REGEX MATCH "([^,]*),(INCLUDE|EXCLUDE),(.*)" _tmp "${_args}")
         list(FILTER CMAKE_MATCH_1 ${CMAKE_MATCH_2} REGEX "${CMAKE_MATCH_3}")
         string(REPLACE "${_match}" "${CMAKE_MATCH_1}" _in "${_in}")
      endif ()
   endforeach()

   string(REPLACE "${_semicol}" ";" _in "${_in}")
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


#  _in must be a string  ### TODO: INTERFACE
function ( eval_genex _list _lang _comp )

   #
   # Optional arguments
   #
   set(_option_arg STRING)
   set(_one_value_args COMP_VERSION CONFIG INTERFACE)
   cmake_parse_arguments( ARG "${_option_arg}" "${_one_value_args}" "" ${ARGN} )

   # Auxiliary list
   set(_in "${${_list}}")

   #
   # Loop to deal with nested genex if any are present
   #
   foreach ( iter RANGE 10 )

      # Check if there are any genex-es in the input. If not, exit foreach
      string(GENEX_STRIP "${_in}" _nogenex_string)
      if (_nogenex_string STREQUAL _in)
         break ()
      endif ()

      #
      # String variable queries
      #

      # Genex in the form $<CONFIG>
      # Ignore if no CONFIG arg is given? NOTE: should I keep it or ignore it?
      if (ARG_CONFIG)
         string(REGEX REPLACE "\\$<CONFIG>" "${ARG_CONFIG}"  _in "${_in}")
      endif ()

      # Genex in the form $<*_COMPILER_ID>
      string(REGEX REPLACE "\\$<${_lang}_COMPILER_ID>" "${_comp}"  _in "${_in}")
      string(REGEX REPLACE "\\$<[A-Za-z]*_COMPILER_ID>" ""  _in "${_in}")

      # Genex in the form $<*_COMPILER_VERSION>
      if (ARG_COMP_VERSION)
         string(REGEX REPLACE "\\$<${_lang}_COMPILER_VERSION>" "${ARG_COMP_VERSION}>" "1"  _in "${_in}")
      endif ()
      string(REGEX REPLACE "\\$<[A-Za-z]*_COMPILER_VERSION>" "0"  _in "${_in}")

      # Genex in the form $<COMPILE_LANGUAGE>
      string(REGEX REPLACE "\\$<COMPILE_LANGUAGE>" "${_lang}"  _in "${_in}")

      # Genex in the form $<COMPILE_LANGUAGE>
      string(REGEX REPLACE "\\$<LINK_LANGUAGE>" "${_lang}"  _in "${_in}")


      #
      #  Boolean variable queries
      #

      # Genex in the form $<CONFIG:cfg>
      # Ignore if no CONFIG arg is given? NOTE: should I keep it or ignore it?
      if (ARG_CONFIG)
         string(REGEX REPLACE "\\$<CONFIG:${ARG_CONFIG}>" "1"  _in "${_in}")
      endif ()
      string(REGEX REPLACE "\\$<CONFIG:[A-Za-z]*>" "0"  _in "${_in}")

      # Genex in the form $<*_COMPILER_ID:compiler_ids>
      string(REGEX REPLACE "\\$<${_lang}_COMPILER_ID:[^>]*${_comp}[^>]*>" "1"  _in "${_in}")
      string(REGEX REPLACE "\\$<[A-Za-z]*_COMPILER_ID:[A-Za-z]*>" "0"  _in "${_in}")

      # Genex in the form $<*_COMPILER_VERSION:version>
      # Ignore if no COMP_VERSION arg is given
      if (ARG_COMP_VERSION)
         string(REGEX REPLACE "\\$<${_lang}_COMPILER_VERSION:${ARG_COMP_VERSION}>" "1"  _in "${_in}")
      endif ()
      string(REGEX REPLACE "\\$<[A-Za-z]*_COMPILER_VERSION:[^\$>]*>" "0"  _in "${_in}")

      # Genex in the form $<COMPILE_LANGUAGE_AND_ID:language,compiler_ids>
      string(REGEX REPLACE "\\$<COMPILE_LANGUAGE_AND_ID:${_lang},[^\$>]*${_comp}[^\$>]*>" "1"  _in "${_in}")
      string(REGEX REPLACE "\\$<COMPILE_LANGUAGE_AND_ID:[A-Za-z,]*>" "0"  _in "${_in}")

      # Genex in the form $<COMPILE_LANGUAGE:languages>
      string(REGEX REPLACE "\\$<COMPILE_LANGUAGE:[^\$>]*${_lang}[^\$>]*>" "1"  _in "${_in}")
      string(REGEX REPLACE "\\$<COMPILE_LANGUAGE:[A-Za-z]*>" "0"  _in "${_in}")

      # Genex in the form $<LINK_LANGUAGE_AND_ID:language,compiler_ids>
      string(REGEX REPLACE "\\$<LINK_LANGUAGE_AND_ID:${_lang},[^\$>]*${_comp}[^\$>]*>" "1"  _in "${_in}")
      string(REGEX REPLACE "\\$<LINK_LANGUAGE_AND_ID:[A-Za-z,]*>" "0"  _in "${_in}")

      # Genex in the form $<LINK_LANGUAGE:languages>
      string(REGEX REPLACE "\\$<LINK_LANGUAGE:[^\$>]*${_lang}[^\$>]*>" "1"  _in "${_in}")
      string(REGEX REPLACE "\\$<LINK_LANGUAGE:[A-Za-z]*>" "0"  _in "${_in}")

      # String transformation
      eval_string_transformations( _in "${_in}")

      # String comparisons
      eval_string_comparisons(_in "${_in}")

      # Logical operators
      eval_logical_operators(_in "${_in}")

      # Conditional expressions
      eval_conditional_expressions(_in "${_in}")

      # Deal with interface if present
      if (ARG_INTERFACE)
         string(REGEX REPLACE "\\$<${ARG_INTERFACE}_INTERFACE:([^\$>]*)>" "\\1"  _in "${_in}")
      endif ()
      string(REGEX REPLACE "\\$<[A-Za-z]*_INTERFACE:[^\$>]*>" ""  _in "${_in}")

      # Remove empty elements, i.e. the ";;" resulting from
      # deleting false genex.
      # In this case, _in is treated as a list
      list(REMOVE_ITEM _in "")

      # More clean up
      list(TRANSFORM _in STRIP)

   endforeach ()

   # If string is required
   if (ARG_STRING)
      string(REPLACE ";" " " _in "${_in}" )
      string(STRIP "${_in}" _in)
   endif ()

   set(${_list} "${_in}" PARENT_SCOPE)

endfunction ()
