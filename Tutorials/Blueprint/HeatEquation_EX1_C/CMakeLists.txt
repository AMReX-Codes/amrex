if (AMReX_SPACEDIM EQUAL 1 OR NOT CMAKE_Fortran_COMPILER_LOADED)
   return()
endif ()

set(_sources advance_${AMReX_SPACEDIM}d.f90 advance.cpp init_phi_${AMReX_SPACEDIM}d.f90 main.cpp myfunc_F.H  myfunc.H)
list(TRANSFORM _sources PREPEND Source/)

set(_input_files Exec/inputs_${AMReX_SPACEDIM}d)

setup_tutorial(_sources _input_files)

unset(_sources)
unset(_input_files)
