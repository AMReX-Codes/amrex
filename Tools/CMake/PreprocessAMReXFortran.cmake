function(preprocess_boxlib_fortran OUTVAR)

   set_source_files_properties(${ARGN} PROPERTIES COMPILE_DEFINITIONS "BL_LANG_FORT")
   set_source_files_properties(${ARGN} PROPERTIES COMPILE_FLAGS "-ffixed-line-length-none")
   set(${OUTVAR} ${ARGN} PARENT_SCOPE)

endfunction(preprocess_boxlib_fortran)
