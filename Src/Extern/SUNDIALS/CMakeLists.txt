target_sources( amrex
   PRIVATE
   arkode_interface.f90
   cvode_interface.f90
   )

target_compile_definitions(amrex PUBLIC AMREX_USE_SUNDIALS)
