set(_sources main.cpp MyTest.cpp MyTest.H MyTestPlotfile.cpp)

file( GLOB_RECURSE _input_files LIST_DIRECTORIES false ${CMAKE_CURRENT_LIST_DIR}/input*)

setup_tutorial(_sources _input_files)

unset(_sources)
unset(_input_files)
