#[=======================================================================[
AMReXTypecheck
----------------

Provides the function add_typecheck_target(_target).

add_typecheck_target(_target) allows to check consistency of data types between
C++ signature of Fortran routines and actual implementation of these same
routines.

If AMReX is included in the build-tree of another project as a subproject,
add_typecheck_target(_target) typechecks the source code of both ${_target}
and amrex. If AMReX is an imported target, i.e. it has been included via
a find_package(AMReX) call, add_typecheck_target(_target) typechecks the
source code of ${_target} alone.

add_typecheck_target(_target) assumes that all the C++ interfaces to
Fortran routines are located in header files with name ending in "_f.H"
or "_F.H". It also assumes that the Fortran implementations are stored
in files with extentions ".f", ".f90", ".F", and ".F90". 


Typechecking is performed in 4 steps:

1) All fortran sources are compiled. This is necessary so that all
".mod" files are available for the next steps.

2) Preprocess C++ header files with name ending in "_f.H" or "_F.H".
Store the preprocessed header in "<header_name>-cppd.h".

3) Output the internal parse tree for each Fortran source.
This steps requires the module file generated at step 1).
Internal parse tree for file "<fortran-source>" is stored in
"<fortran-source>.orig"

4) Call AMReX "typechecker.py" to compare interfaces defined in "-cppd.h"
files (step 2) and internal parse tree stored in ".orig" files (step 3)


add_typecheck_target(_target) adds two targets:

* typecheckobjs_${_target}: perform step 1)
* typecheck_${_target}: perform step 2), 3), and 4)
#]=======================================================================]
function( add_typecheck_target _target)

   # 
   # Check if we have all we need to define the typecheck target
   # 
   if (  NOT (CMAKE_Fortran_COMPILER_ID MATCHES GNU) OR
         NOT (CMAKE_C_COMPILER_ID MATCHES GNU)       OR
         NOT (CMAKE_Fortran_COMPILER_ID MATCHES GNU) )
      message(WARNING "Typecheck disabled because compiler ID is not GNU")
      return ()
   endif ()

   
   # We need either "AMReX::amrex", if imported target, or "amrex", if included as subproject.
   if ( (NOT TARGET AMReX::amrex ) AND (NOT TARGET amrex) )
      message(WARNING "No valid amrex target exists: typecheck is disabled")
      return()
   endif ()

   if ( NOT TARGET ${_target} )
      message(WARNING "Skipping typecheck for target ${_target} bacause it does not exist")
      return ()
   endif ()

   find_package(Python3 COMPONENTS Interpreter Development QUIET)
   if (NOT Python3_FOUND)
      message(WARNING "Typecheck disabled because Python 3 was not found")
      return ()
   endif ()

   #
   # Set directory for typecheck
   #
   set( _typecheck_dir  "${CMAKE_CURRENT_BINARY_DIR}/TypeCheckTemp/${_target}" )

   #
   # Define which form of amrex target we are dealing with 
   #
   if (TARGET amrex)
      # AMReX is included as subproject
      set(_amrex_target amrex)
   elseif (TARGET AMReX::amrex)
      # AMReX is an imported target
      set(_amrex_target AMReX::amrex)
   endif ()

   #
   # Get the fortran sources and the fortran-interfaces headers from _target
   #
   get_target_property( _sources ${_target} SOURCES )

   # If amrex is included as a subproject, get its Fortran sources too
   if (TARGET amrex)
      get_target_property( _amrex_sources amrex SOURCES )
      list(APPEND _sources ${_amrex_sources})
   endif ()
   
   set(_fsources ${_sources})
   list(FILTER _fsources INCLUDE REGEX "(\.f90|\.f|\.F90|\.F)$" )

   set(_fheaders ${_sources})
   list(FILTER _fheaders INCLUDE REGEX "(\_f\.H|\_F\.H)$")

   #
   # Find includes and defines required for those sources
   # Must be done manually since we will use them in a custom command
   #
   set(_includes ${_typecheck_dir})
   get_target_property( _target_includes ${_target}        INCLUDE_DIRECTORIES )
   get_target_property( _amrex_includes  ${_amrex_target}  INTERFACE_INCLUDE_DIRECTORIES )
   get_target_property( _amrex_mods      ${_amrex_target}  Fortran_MODULE_DIRECTORY )  

   foreach (_item IN LISTS _target_includes _amrex_includes _amrex_mods)
      if (_item)
         list(APPEND _includes  ${_item})
      endif ()
   endforeach ()

   list( REMOVE_DUPLICATES _includes )

   
   set(_defines AMREX_TYPECHECK) # this needs to be added
   get_target_property( _target_defines ${_target}        COMPILE_DEFINITIONS )
   get_target_property( _amrex_defines  ${_amrex_target}  INTERFACE_COMPILE_DEFINITIONS )

   foreach (_item IN LISTS _target_defines _amrex_defines)
      if (_item)
         list(APPEND _defines  ${_item})
      endif ()
   endforeach ()

   list( REMOVE_DUPLICATES _defines )
  
   # If amrex is included as a subproject, get its required compiler flags
   # If amrex is a subproject, fortran sources need Flags_Fortran_REQUIRED
   set(_amrex_flags)
   if (TARGET amrex)
      get_target_property( _amrex_flags Flags_Fortran_REQUIRED INTERFACE_COMPILE_OPTIONS )
   endif ()

   
   # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # STEP 1: create fortran modules
   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   #
   # Define the typecheckobjs library.
   # Basically, we only want to generate the module files
   # associated with fortran otherwise the other two steps
   # of the type check (see below) will fail because
   # it will not be able to include symbols
   # if modules are not there.
   #

   # ---------------------------------------------------------------------------------
   # NOTE:
   # ---------------------------------------------------------------------------------
   # We could skip generating the objs file and create only the modules
   # by setting the following compile option:
   # target_compile_options(typecheckobjs PRIVATE -fsyntax-only)
   # However, since this would only generate ".mod" files and no ".o" files, the target
   # would always be out of date and repeat itself.
   # A work around would be to create the module and the orig file at the same time.
   # this could be achieved with a add_custom_command which is aware of dependecies info.
   # To this end, we could use IMPLICIT_DEPENDS. Unfortunately this option is supported
   # for C and C++ only for now.
   set(_typecheckobjs  typecheckobjs_${_target})

   add_library( ${_typecheckobjs} OBJECT EXCLUDE_FROM_ALL )

   set_target_properties( ${_typecheckobjs}
      PROPERTIES
      SOURCES                   "${_fsources}"
      INCLUDE_DIRECTORIES       "${_includes}"
      Fortran_MODULE_DIRECTORY  "${_typecheck_dir}"
      COMPILE_DEFINITIONS       "${_defines}"
      COMPILE_OPTIONS           "${_amrex_flags}"
      )

   # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # STEP 2: create CPPD files from C++ headers
   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   # Manually setup includes and defines to use in custom commands
   if (_includes)
      string(REPLACE ";" ";-I" _includes "-I${_includes}")
   endif ()

   # expand genex-es
   include(AMReXGenexHelpers)
   
   evaluate_genex(_defines _cxx_defines
      LANG CXX
      COMP GNU
      CONFIG ${CMAKE_BUILD_TYPE}
      INTERFACE BUILD)
   
   evaluate_genex(_defines _fortran_defines
      LANG Fortran
      COMP GNU
      CONFIG ${CMAKE_BUILD_TYPE}
      INTERFACE BUILD)

   evaluate_genex(_amrex_flags _amrex_fortran_flags
      LANG Fortran
      COMP GNU
      CONFIG ${CMAKE_BUILD_TYPE}
      INTERFACE BUILD) 
   
   if (_cxx_defines)
      string(REPLACE ";" ";-D" _cxx_defines "-D${_cxx_defines}")
   endif ()

   if (_fortran_defines)
      string(REPLACE ";" ";-D" _fortran_defines "-D${_fortran_defines}")
   endif ()

   # Set proper replacements for real types
   set(AMREX_REAL double)
   if ( ("-DAMREX_USE_FLOAT" IN_LIST _cxx_defines)   OR
        (DEFINED AMReX_DP_FOUND AND NOT AMReX_DP_FOUND) )
      set(AMREX_REAL float)
   endif ()

   set(AMREX_PARTICLE_REAL double)
   if (  ("-DAMREX_SINGLE_PRECISION_PARTICLES" IN_LIST _cxx_defines) OR
         (DEFINED AMReX_DPARTICLES_FOUND AND NOT AMReX_DPARTICLES_FOUND) )
      set(AMREX_PARTICLE_REAL float)
   endif ()

   set (_cppd)
   foreach ( _file IN LISTS _fheaders )
      get_filename_component( _fname    ${_file} NAME )     # This strips away the path
      get_filename_component( _fullname ${_file} ABSOLUTE ) # This add the absolute path to fname
      set( _cppd_file ${_typecheck_dir}/${_fname}-cppd.h )
      add_custom_command(
         OUTPUT  ${_cppd_file}
         COMMAND ${CMAKE_C_COMPILER}
         ARGS    ${_cxx_defines} ${_includes} -E -P -x c -std=c99 ${_fullname} > ${_cppd_file}
         COMMAND sed
         ARGS -i -e 's/amrex::Real/${AMREX_REAL}/g' ${_cppd_file}
         COMMAND sed
         ARGS -i -e 's/amrex_real/${AMREX_REAL}/g' ${_cppd_file}
         COMMAND sed
         ARGS -i -e 's/amrex_particle_real/${AMREX_PARTICLE_REAL}/g' ${_cppd_file}
         COMMAND sed
         ARGS -i -e '/typedef\\s*${AMREX_REAL}/d' ${_cppd_file}
         COMMAND sed
         ARGS -i -e 's/AMREX_GPU_DEVICE/ /g' ${_cppd_file}
         COMMAND sed
         ARGS -i -e 's/AMREX_GPU_HOST_DEVICE/ /g' ${_cppd_file}
         COMMAND sed
         ARGS -i -e 's/\\&/*/g' ${_cppd_file}
         DEPENDS ${_file}
         WORKING_DIRECTORY ${_typecheck_dir}
         COMMENT "Generating ${_fname}-cppd.h" )
      list(APPEND _cppd ${_cppd_file})
   endforeach ()


   # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # STEP 3: generate orig files from fortran sources
   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   set(_orig)
   set(_fortran_flags -fsyntax-only -fdump-fortran-original ${_amrex_fortran_flags})
   foreach ( _file IN LISTS _fsources )
      get_filename_component( _fname    ${_file} NAME )     # This strips away the path
      get_filename_component( _fullname ${_file} ABSOLUTE ) # This add the absolute path to fname
      set( _orig_file ${_typecheck_dir}/${_fname}.orig )     # Use fullpath otherwise it rebuilds everytime
      add_custom_command(
         OUTPUT   ${_orig_file}
         COMMAND  ${CMAKE_Fortran_COMPILER}
         ARGS     ${_fortran_defines} ${_includes} ${_fortran_flags}  ${_fullname} > ${_orig_file}
         DEPENDS  ${_file} ${_typecheckobjs}
         WORKING_DIRECTORY    ${_typecheck_dir}
         COMMENT  "Generating ${_fname}.orig" )
      list(APPEND _orig ${_orig_file})
   endforeach ()

   #
   # Add typecheck target
   #
   set(_outfile  "${_typecheck_dir}/${_target}_typecheck.ou" )

   # Find typechecker 
   find_file(_typechecker "typechecker.py"
      HINTS ${AMReX_SOURCE_DIR} ${AMReX_ROOT} ENV AMReX_ROOT PATH_SUFFIXES Tools/typechecker)

   add_custom_target( typecheck_${_target}
      COMMAND python3  ${_typechecker}
      --workdir ${_typecheck_dir} --output ${_outfile}
      DEPENDS ${_cppd} ${_orig}
      WORKING_DIRECTORY ${_typecheck_dir}
      )

   # Add rules to remove typecheck outfile when cleaning
   set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${_outfile})

endfunction ()
