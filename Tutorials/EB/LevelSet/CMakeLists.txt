if (NOT DIM EQUAL 3)
   return()
endif ()

set(_sources main.cpp make_cylinder.cpp make_eb.cpp make_poly.cpp make_shapes.H)
list(TRANSFORM _sources PREPEND Src/)

set(_input_files Exec/inputs Exec/inputs_eb2)

setup_tutorial(_sources _input_files)

unset(_sources)
unset(_input_files)
