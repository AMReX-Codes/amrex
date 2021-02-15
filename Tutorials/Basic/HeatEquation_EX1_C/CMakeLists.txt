if (AMReX_SPACEDIM EQUAL 1)
   return()
endif ()

# List of source files
set(_sources main.cpp myfunc.cpp myfunc.H mykernel.H)
list(TRANSFORM _sources PREPEND "Source/")

# List of input files
file( GLOB_RECURSE _input_files LIST_DIRECTORIES false Exec/input* )

setup_tutorial(_sources _input_files)

unset( _sources )
unset( _input_files   )
