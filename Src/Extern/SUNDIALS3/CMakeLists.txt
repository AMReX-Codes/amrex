if (ENABLE_FORTRAN_INTERFACES)
   target_sources( amrex
      PRIVATE
      fnvector_serial.f90
      fnvector_serial_fprefix.f90
      fsunlinsol_dense.f90
      fsunmat_dense.f90
      # Arkcode interfaces
      arkode_interface.f90
      farkode.f90
      # CVODE interfaces
      cvode_interface.f90
      fcvode.f90
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

find_package(SUNDIALS REQUIRED COMPONENTS ${_sundials_components})

foreach(_comp ${_sundials_components})
   target_link_libraries(amrex PUBLIC SUNDIALS::${_comp})   
endforeach ()

target_compile_definitions(amrex PUBLIC AMREX_USE_SUNDIALS_3x4x)
   


