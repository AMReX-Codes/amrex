if ( (NOT AMReX_EB) OR NOT (AMReX_SPACEDIM EQUAL 3))
   return()
endif ()

set(_sources main.cpp)
set(_input_files inputs_3d)

setup_tutorial(_sources _input_files)

unset(_sources)
unset(_input_files)
