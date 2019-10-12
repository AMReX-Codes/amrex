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
