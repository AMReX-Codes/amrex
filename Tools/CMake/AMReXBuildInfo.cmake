include(AMReXTargetHelpers)

macro (generate_buildinfo _target )
   
   cmake_parse_arguments("_arg" "REQUIRED" "" "" ${ARGN})

   # 
   # check if _target is a valid target
   # 
   if (NOT TARGET ${_target} )
      message(FATAL_ERROR "Target ${_target} does not exists")
   endif ()
   
   #
   # Look for python
   # 
   find_package(Python)

   # If Python >= 2.7 is available, generate AMReX_BuildInfo.cpp
   # If Python is not available, do not include AMReX_BuildInfo.cpp in library
   # AMReX_Buildinfo.cpp is optional, not required.
   #
   if ( Python_Interpreter_FOUND AND (NOT (Python_VERSION VERSION_LESS "2.7") ) )


      get_target_properties_flattened(${_target} _includes _defines _flags _link_line)

      set(_lang CXX Fortran)
      set(_prop _includes _defines _flags _link_line )

      foreach( _p IN LISTS _prop)
         foreach(_l IN LISTS _lang)

            string(TOLOWER ${_l} _ll)
            
            evaluate_genex(${_p} _${_ll}${_p}
               LANG ${_l}
               COMP ${CMAKE_${_l}_COMPILER_ID}
               CONFIG ${CMAKE_BUILD_TYPE}
               INTERFACE BUILD)

            if (_${_ll}${_p})

               list(REMOVE_DUPLICATES _${_ll}${_p})
               
               if ("${_p}" STREQUAL "_defines")
                  string(REPLACE ";" " -D" _${_ll}${_p} "-D${_${_ll}${_p}}")
               elseif ("${_p}" STREQUAL "_includes")
                  string(REPLACE ";" " -I" _${_ll}${_p} "-I${_${_ll}${_p}}")
               else()
                  string(REPLACE ";" " " _${_ll}${_p} "${_${_ll}${_p}}")
               endif ()              

                  print(_${_ll}${_p})

            endif ()

            unset(_ll)

            
         endforeach()
      endforeach ()

      string(TOUPPER ${CMAKE_BUILD_TYPE} _ubuild_type)
      set(_cxx_flags "${CMAKE_CXX_FLAGS_${_ubuild_type}} ${CMAKE_CXX_FLAGS} ${_cxx_flags}")
      set(_fortran_flags "${CMAKE_Fortran_FLAGS_${_ubuild_type}} ${CMAKE_Fortran_FLAGS} ${_fortran_flags}")
      
      # Determine whether we are using AMReX as a library or as a sub-project
      if (TARGET amrex)   # amrex included as subproject via add_subdirectory
         
      elseif (TARGET AMReX::amrex) # amrex used as library
         string(REPLACE "CMake" "" _c_scripts_dir "${AMReX_MODULES_PATH}")
         string(REPLACE "/Tools" "" _amrex_home    "${_c_scripts_dir}")
         set(_c_scripts_dir "${_c_scripts_dir}C_scripts/")
      else ()

      endif ()

      set(_script "${_c_scripts_dir}makebuildinfo_C.py")
      



      add_custom_command(
         COMMAND ${Python_EXECUTABLE} ${_script}
                 --amrex_home ${_amrex_home}    
                 #  CXX
                 --COMP ${CMAKE_CXX_COMPILER_ID}
                 --COMP_VERSION ${CMAKE_CXX_COMPILER_VERSION}
                 --CXX_comp_name ${CMAKE_CXX_COMPILER}
                 --CXX_flags "${_ccx_defines} ${_cxx_includes} ${_cxx_flags}"
                 # Fortran
                 --FCOMP ${CMAKE_Fortran_COMPILER_ID}
                 --FCOMP_VERSION ${CMAKE_Fortran_COMPILER_VERSION}
                 --F_comp_name ${CMAKE_Fortran_COMPILER}
                 --F_flags "${_fortran_defines} ${_fortran_includes} ${_fortran_flags}"
                 #--link_flags
                 --libraries "${_cxx_link_line}"
         OUTPUT AMReX_buildInfo.cpp
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
         COMMENT "Generating AMReX_buildInfo.cpp" )
      
      target_sources( ${_target}
         PRIVATE
         ${CMAKE_CURRENT_BINARY_DIR}/AMReX_buildInfo.cpp
         )

      target_sources( ${_target}
         PRIVATE
         ${_c_scripts_dir}AMReX_buildInfo.H
         )

      target_include_directories( ${_target}
         PUBLIC
         $<BUILD_INTERFACE:${_c_scripts_dir}>
         )


      #
      # Clean-up
      #
      unset(_ubuild_type)
      unset(_amrex_home)
      unset(_c_scripts_dir)
      unset(_script)

      foreach( _p IN LISTS _prop)
         foreach(_l IN LISTS _lang)
            unset(_${_ll}${_p})           
         endforeach()
      endforeach ()

      unset(_prop)
      unset(_lang)
      
   else ()
      
      if (_arg_REQUIRED)
         message(FATAL_ERROR "Cannot find Python > 2.7")
      else ()
         message(WARNING "Cannot find Python > 2.7: skipping generation of build info")
      endif ()
      
   endif ()

   unset(_arg_REQUIRED)
   
endmacro ()
