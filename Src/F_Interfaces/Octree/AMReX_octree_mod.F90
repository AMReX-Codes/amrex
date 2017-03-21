module amrex_octree_module

  use iso_c_binding
  use amrex_base_module
  use amrex_amrcore_module

  implicit none

  private

  public :: amrex_octree_init, amrex_octree_finalize

  type, bind(c) :: node
     integer :: level, grid
  end type node

  type(node), allocatable :: leaf_nodes(:)

  interface
     subroutine amrex_fi_init_octree () bind(c)
     end subroutine amrex_fi_init_octree

     subroutine amrex_fi_build_octree_leaves (amr, nleaves, leaves) bind(c)
       import
       implicit none
       integer, intent(out) :: nleaves
       type(c_ptr), intent(in), value :: amr
       type(c_ptr), intent(out) :: leaves
     end subroutine amrex_fi_build_octree_leaves

     subroutine amrex_fi_copy_octree_leaves (leaves, leaves_copy) bind(c)
       import
       implicit none
       type(c_ptr), value :: leaves
       type(node), intent(inout) :: leaves_copy(*)
     end subroutine amrex_fi_copy_octree_leaves
  end interface

contains

  subroutine amrex_octree_init ()
    call amrex_fi_init_octree()
    call amrex_init_post_regrid_function(amrex_octree_post_regrid)
  end subroutine amrex_octree_init

  subroutine amrex_octree_finalize ()
    if (allocated(leaf_nodes)) deallocate(leaf_nodes)
  end subroutine amrex_octree_finalize

  subroutine amrex_octree_post_regrid ()
    type(c_ptr) :: amrcore, leaves
    integer :: nleaves
    amrcore = amrex_get_amrcore()
    call amrex_fi_build_octree_leaves(amrcore, nleaves, leaves)
    if (allocated(leaf_nodes)) deallocate(leaf_nodes)
    allocate(leaf_nodes(nleaves));
    call amrex_fi_copy_octree_leaves(leaves, leaf_nodes);
  end subroutine amrex_octree_post_regrid

end module amrex_octree_module
