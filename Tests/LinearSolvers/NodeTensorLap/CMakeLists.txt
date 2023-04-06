if (NOT AMReX_CUDA AND NOT AMReX_SPACEDIM EQUAL 1)
set(_sources main.cpp MyTest.cpp MyTest.H MyTestPlotfile.cpp)
set(_input_files)

setup_test(_sources _input_files )

unset(_sources)
unset(_input_files)
endif ()
