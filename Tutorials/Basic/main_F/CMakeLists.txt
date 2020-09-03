if (NOT ENABLE_FORTRAN_INTERFACES)
   return()
endif ()

set(_sources main.F90)
set(_include_dirs)
set(_input_files)

setup_tutorial(_sources _include_dirs _input_files)

unset(_sources)
unset(_include_dirs)
unset(_input_files)
