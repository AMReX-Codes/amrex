set(_sources     main.cpp)
if (AMReX_CUDA)
  set(_input_files inputs.rt.cuda  )
else ()
  set(_input_files inputs.rt  )
endif ()

setup_test(_sources _input_files)

unset(_sources)
unset(_input_files)
