#
# Fortran interfaces to SUNDIALS routines
#
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

target_compile_definitions(amrex PUBLIC AMREX_USE_SUNDIALS_3x4x)
