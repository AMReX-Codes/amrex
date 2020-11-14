if (NOT (DIM EQUAL 3))
   return()
endif ()

set(_sources MyKernel_F.H main.cpp )
set(_input_files)

setup_tutorial(_sources _input_files)

unset(_sources)
unset(_input_files)
