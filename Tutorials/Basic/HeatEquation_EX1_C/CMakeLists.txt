set(_sources main.cpp myfunc.cpp myfunc.H mykernel.H)
list(TRANSFORM _sources PREPEND "${CMAKE_CURRENT_LIST_DIR}/Source/")

set(_include_dirs ${CMAKE_CURRENT_LIST_DIR}/Source)

file( GLOB_RECURSE _input_files LIST_DIRECTORIES false ${CMAKE_CURRENT_LIST_DIR}/Exec/input* )

setup_tutorial(_sources _include_dirs _input_files)

unset( _sources )
unset( _include_dirs  )
unset( _input_files   )
