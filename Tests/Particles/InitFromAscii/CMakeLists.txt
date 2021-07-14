set(_sources     main.cpp)
set(_input_files inputs  )

setup_test(_sources _input_files)

file(COPY particles.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

unset(_sources)
unset(_input_files)
