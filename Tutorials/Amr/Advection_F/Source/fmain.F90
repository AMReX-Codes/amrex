
program main

  use amrex_amr_module

  use my_amr_module
  use initdata_module
  use evolve_module
  use amrex_amrtracerparticlecontainer_module
  
  implicit none

  type(c_ptr) :: amrcore = c_null_ptr
  
  call amrex_init()
  call amrex_amrcore_init()

  call my_amr_init()

  call initdata()

  amrcore = amrex_get_amrcore()
  call amrex_amrtracerparticlecontainer_init(amrcore)
  call amrex_init_particles_one_per_cell()
  
  call evolve()

  call amrex_amrtracerparticlecontainer_finalize()
  
  call my_amr_finalize()

  call amrex_amrcore_finalize()
  call amrex_finalize()

end program main
