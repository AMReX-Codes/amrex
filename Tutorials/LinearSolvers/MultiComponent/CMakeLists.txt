set(_sources main.cpp MCNodalLinOp.cpp MCNodalLinOp.H)

file( GLOB_RECURSE _input_files LIST_DIRECTORIES false inputs)

setup_tutorial(_sources _input_files)

unset(_sources)
unset(_input_files)
