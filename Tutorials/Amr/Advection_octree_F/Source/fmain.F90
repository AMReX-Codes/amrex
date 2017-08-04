
program main

  use amrex_amr_module
  use amrex_octree_module

  use my_amr_module
  use initdata_module
  use evolve_module

  implicit none

  call amrex_init()
  call amrex_octree_init()
  call amrex_amrcore_init()

  call my_amr_init()

  call initdata()

  call evolve()

  call my_amr_finalize()

  call amrex_amrcore_finalize()
  call amrex_octree_finalize()
  call amrex_finalize()

end program main
