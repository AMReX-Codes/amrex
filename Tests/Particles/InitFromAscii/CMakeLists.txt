if ( NOT (3 IN_LIST AMReX_SPACEDIM) OR AMReX_MPI )
   return ()
endif ()

set(_sources     main.cpp)
set(_input_files inputs  )

setup_test(3 _sources _input_files)

file(COPY particles.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/3d)

unset(_sources)
unset(_input_files)
