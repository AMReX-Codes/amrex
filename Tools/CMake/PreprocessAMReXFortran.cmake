# -*- mode: cmake -*-
# - Pass a list of Fortran files through a BoxLib-specific preprocessing 
#
# PREPROCESS_BOXLIB_FORTRAN( OUTVAR source1 ... sourceN )
#
#  OUTVAR  A list containing all the output file names, suitable
#          to be passed to add_executable or add_library.
#
# The output files are placed in the same relative location to
# CMAKE_CURRENT_BINARY_DIR as they are to CMAKE_CURRENT_SOURCE_DIR.
#
# Example:
#  preprocess_boxlib_fortran( SRCS src/test1.F src/test2.F )
#  add_executable( test ${SRCS} )

if(__PREPROCESS_BOXLIB_FORTRAN_INCLUDED)
  return()
endif()
set(__PREPROCESS_BOXLIB_FORTRAN_INCLUDED)

find_package(Perl REQUIRED)

function(preprocess_boxlib_fortran OUTVAR)

  if(MSVC)
    set(BL_DEFINE_FLAG "/D")
  else(MSVC)
    set(BL_DEFINE_FLAG "-D")
  endif(MSVC)

  get_directory_property(cmpdefs COMPILE_DEFINITIONS)
  set(defflags "${BL_DEFINE_FLAG}BL_LANG_FORT")
  foreach (d ${cmpdefs})
    list(APPEND defflags "${BL_DEFINE_FLAG}${d}")
  endforeach (d ${compdefs})

  get_directory_property(incldirs INCLUDE_DIRECTORIES)
  set(inclflags)
  foreach (i ${incldirs})
    list(APPEND inclflags "${CMAKE_INCLUDE_FLAG_C}${i}${CMAKE_INCLUDE_FLAG_C_SEP}")
  endforeach (i ${incldirs})

  set(outfiles)
  foreach( f ${ARGN} )

    if(NOT IS_ABSOLUTE "${f}")
      get_filename_component(f "${f}" ABSOLUTE)
    endif()
    file(RELATIVE_PATH r "${CMAKE_CURRENT_SOURCE_DIR}" "${f}")
    get_filename_component(n "${r}" NAME_WE)
    get_filename_component(p "${r}" PATH)
    get_filename_component(e "${r}" EXT)
    if(MSVC)
      set(of "${CMAKE_CURRENT_BINARY_DIR}\${p}${n}${e}.f")
      add_custom_command(
        OUTPUT ${of}
        COMMAND fpp /m /ansi ${inclflags} ${defflags} ${f} |
                 ${PERL_EXECUTABLE} ${CCSE_PERL_DIR}/strip72 -c > ${of}
        DEPENDS ${f}
        COMMENT "Preprocessing ${f}..."
      )
    else(MSVC)
      set(of "${CMAKE_CURRENT_BINARY_DIR}/${p}${n}${e}.f")
      set(PREPROCESS_FLAGS)
      if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        list(APPEND PREPROCESS_FLAGS "-traditional")
      elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        list(APPEND PREPROCESS_FLAGS "-traditional")
      endif()

      add_custom_command(
        OUTPUT ${of}
        COMMAND ${CMAKE_C_COMPILER} -E ${PREPROCESS_FLAGS} ${inclflags} ${defflags} ${f} |
                 ${PERL_EXECUTABLE} ${CCSE_PERL_DIR}/strip72 -c > ${of}
        DEPENDS ${f}
        COMMENT "Preprocessing ${f}..."
      )
    endif(MSVC)
    set_source_files_properties(${of} PROPERTIES GENERATED ON)
    list(APPEND outfiles "${of}")
  endforeach(f)
  set(${OUTVAR} ${outfiles} PARENT_SCOPE)

endfunction(preprocess_boxlib_fortran)
