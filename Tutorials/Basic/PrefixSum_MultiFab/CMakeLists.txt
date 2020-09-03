# List of source files
set(_sources main.cpp)

# List of include directories
set(_include_dirs)

# List of input files
set(_input_files inputs Parallel-Prefix-Sum.ipynb)

setup_tutorial(_sources _include_dirs _input_files)

unset( _sources )
unset( _include_dirs  )
unset( _input_files   )
