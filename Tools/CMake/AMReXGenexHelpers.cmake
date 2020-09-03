#[=======================================================================[
AMReXGenexHelpers
-----------------

Provides functions to evaluate generator expressions.

The most important of these functions is

  * eval_genex( _list _lang _comp )

eval_genex(...) evaluates the generator expression in _list
according to user-specified parameters such as compiler language (_lang)
and compiler id (_comp).

Other functions provided by this module are:

  * eval_conditional_expressions( _out _in )
  * eval_logical_operators( _out _in )
  * eval_string_comparisons( _out _in )
  * eval_string_transformations( _out _in )

Each of these functions evaluate a specific type of genex in a list
and therefore are rarely needed outside of this module.

IMPORTANT NOTE:
when specifying regex, be aware that the group "\\", i.e. quoted
double backslash, is parsed as a single  backslash.Therefore, "\\$" is
parsed as just \$.

#]=======================================================================]

#
# Evaluate string-valued generator expressions that depend on a boolean condition
# that must be 0 or 1.
#
function ( eval_conditional_expressions _list)

   set(_in ${${_list}})

   # Genex in the form $<condition:true_string> where condition=<0|1>
   #
   string(REGEX REPLACE "\\$<1:([^\$>]*)>" "\\1" _in "${_in}")  # true (\\1 is first subexpression in the match)
   string(REGEX REPLACE "\\$<0:[^\$>]*>"   ""    _in "${_in}")  # false

   #
   # Genex in the form $<IF:condition,true_string,false_string> where condition=<0|1>
   #
   string(REGEX REPLACE "\\$<IF:0,[^\$,]*,([^\$>]*)>"   "\\1" _in "${_in}")  # false
   string(REGEX REPLACE "(\\$<IF):1,([^\$,]*),[^\$>]*>" "\\2" _in "${_in}")  # true

   set(${_list}  "${_in}" PARENT_SCOPE)

endfunction ()



#
# Evaluate boolean generator expressions that perform logical operations
#
function (eval_logical_operators _list)

   set(_in ${${_list}})

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

   set(${_list} "${_in}" PARENT_SCOPE)

endfunction ()

#
# Evaluate boolean generator expressions that perform string operations
#
function (eval_string_comparisons _list)

   set(_in ${${_list}})

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
   set(${_list} "${_in}" PARENT_SCOPE)

endfunction ()

#
# Evaluate string-valued generator expressions that perform string
# transformations
#
function (eval_string_transformations _list)

   set(_in ${${_list}})

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
   set(${_list} "${_in}" PARENT_SCOPE)

endfunction ()

#
#
# FUNCTION: eval_genex
#
# Evaluate the generator expressions in an input list
# according to user-defined criteria (see arguments below).
#
# Arguments:
#   _list         the input list
#   _lang         the compile language
#   _comp         the compiler id
#
# Optional one-value arguments:
#   COMP_VERSION  <compiler-version>
#   CONFIG        <config-type>
#   INTERFACE     <interface-type>
#
# Optional no-value arguments:
#   STRING        Returns a string rather than a list
#
# Additionally, passing the keyword STRING will return a string
# rather than a list
# If an optional argument if not specify, any genex requiring that
# argument will be stripped from _list.
#
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
      eval_string_transformations(_in)

      # String comparisons
      eval_string_comparisons(_in)

      # Logical operators
      eval_logical_operators(_in)

      # Conditional expressions
      eval_conditional_expressions(_in)

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
