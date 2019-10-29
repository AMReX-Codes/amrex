if (ENABLE_FORTRAN_INTERFACES)
   target_sources( amrex
      PRIVATE
      # ARKode interfaces
      arkode_interface.f90
      # CVODE interfaces
      cvode_interface.f90
      )

   list(APPEND _sundials_components farkode_mod fcvode_mod fnvecserial_mod)
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

foreach (_comp ${_sundials_components})
   target_link_libraries(amrex PUBLIC SUNDIALS::${_comp})
endforeach ()

target_compile_definitions(amrex PUBLIC AMREX_USE_SUNDIALS)
