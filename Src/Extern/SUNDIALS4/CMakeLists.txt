set(_sundials_components nvecserial;cvode;arkode)

if (ENABLE_OMP)
   list(APPEND _sundials_components nvecopenmp)
endif ()

if (ENABLE_CUDA)
   list(APPEND _sundials_components nveccuda)
endif ()

if (ENABLE_FORTRAN_INTERFACES)
   add_sources(fsundials_types_mod.f90)
   add_sources(fnvector_mod.f90 fnvector_serial_mod.f90)
   add_sources(fsunmatrix_mod.f90 fsunmatrix_dense_mod.f90)
   add_sources(fsunlinsol_mod.f90 fsunlinsol_dense_mod.f90)

   # ARKode interfaces
   add_sources(arkode_interface.f90 farkode_mod.f90 farkode_arkstep_mod.f90)

   # CVODE interfaces
   add_sources(cvode_interface.f90 fcvode_mod.f90)
endif ()

#
# We link to libraries and always include nvecserial (in case app code needs it)
#
find_package(SUNDIALS REQUIRED COMPONENTS ${_sundials_components})
foreach(_comp ${_sundials_components})
   target_link_libraries(amrex PUBLIC SUNDIALS::${_comp})
endforeach ()
target_compile_definitions(amrex PUBLIC AMREX_USE_SUNDIALS_3x4x)



