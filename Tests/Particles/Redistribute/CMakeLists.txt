set(_sources     main.cpp)
set(_input_files inputs.rt  )  # There are others but we use only this one for now

setup_test(_sources _input_files NTASKS 2)

unset(_sources)
unset(_input_files)
