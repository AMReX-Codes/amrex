if (ENABLE_FORTRAN_INTERFACES)
   target_sources( amrex
      PRIVATE
      fsundials_types_mod.f90
      fnvector_mod.f90
      fnvector_serial_mod.f90
      fsunmatrix_mod.f90
      fsunmatrix_dense_mod.f90
      fsunlinsol_mod.f90
      fsunlinsol_dense_mod.f90
      # ARKode interfaces
      arkode_interface.f90
      farkode_mod.f90
      farkode_arkstep_mod.f90
      # CVODE interfaces
      cvode_interface.f90
      fcvode_mod.f90
      )
endif ()

#
# We link to libraries and always include nvecserial (in case app code needs it)
#
set(_sundials_components nvecserial;cvode;arkode)

if (ENABLE_OMP)
   list(APPEND _sundials_components nvecopenmp)
endif ()

if (ENABLE_CUDA)
   list(APPEND _sundials_components nveccuda)
endif ()

find_package(SUNDIALS 4 REQUIRED COMPONENTS ${_sundials_components})

foreach (_comp ${_sundials_components})
   target_link_libraries(amrex PUBLIC SUNDIALS::${_comp})   
endforeach ()

target_compile_definitions(amrex PUBLIC AMREX_USE_SUNDIALS_3x4x)
