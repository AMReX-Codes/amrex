set(_sources     main.cpp)
set(_input_files inputs.rt)

setup_test(_sources _input_files NTASKS 2)

unset(_sources)
unset(_input_files)
