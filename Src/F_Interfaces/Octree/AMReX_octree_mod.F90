module amrex_octree_module

  use iso_c_binding

  implicit none

  private

  public :: amrex_octree_init, amrex_octree_finalize

  interface
     subroutine amrex_fi_init_octree () bind(c)
     end subroutine amrex_fi_init_octree
  end interface

contains

  subroutine amrex_octree_init ()
    call amrex_fi_init_octree()
  end subroutine amrex_octree_init

  subroutine amrex_octree_finalize ()
  end subroutine amrex_octree_finalize

end module amrex_octree_module
