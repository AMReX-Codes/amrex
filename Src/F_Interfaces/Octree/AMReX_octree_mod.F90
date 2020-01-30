module amrex_octree_module

  use iso_c_binding
  use amrex_base_module
  use amrex_amrcore_module

  implicit none

  private

  public :: amrex_octree_init, amrex_octree_finalize, amrex_octree_iter_build, &
       amrex_octree_iter_destroy, amrex_octree_average_down_leaves

  type, public :: amrex_octree_iter
     integer, private :: begin_index    = 0
     integer, private :: end_index      = 0 ! exclusive
     integer, private :: current_index  = 0
   contains
     procedure :: clear      => amrex_octree_iter_clear
     procedure :: next       => amrex_octree_iter_next
     procedure :: level      => amrex_octree_iter_level
     procedure :: grid_index => amrex_octree_iter_grid_index
     procedure :: box        => amrex_octree_iter_box
  end type amrex_octree_iter

  type, bind(c) :: node
     integer :: level, grid
  end type node

  type(node), allocatable :: leaf_nodes(:)

  interface amrex_octree_average_down_leaves
     module procedure amrex_octree_average_down_leaves_all
     module procedure amrex_octree_average_down_leaves_all_comp
     module procedure amrex_octree_average_down_leaves_level
     module procedure amrex_octree_average_down_leaves_level_comp
  end interface amrex_octree_average_down_leaves

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

     subroutine amrex_fi_octree_average_down_level &
          (amr, flev, fine, coarse, scomp, ncomp) bind(c)
       import
       implicit none
       integer, intent(in), value :: flev, scomp, ncomp
       type(c_ptr), intent(in), value :: amr, fine, coarse
     end subroutine amrex_fi_octree_average_down_level
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

  subroutine amrex_octree_iter_build (oti)
    type(amrex_octree_iter) :: oti
    integer :: nnodes, tid, nthreads, n_per_thread, nlft
    nnodes = size(leaf_nodes)
    tid = omp_get_thread_num()
    nthreads = omp_get_num_threads()
    if (nnodes < nthreads) then  ! there are more threads than nodes
       if (tid < nnodes) then
          oti%begin_index   = tid
          oti%end_index     = tid+1
       else
          oti%begin_index   = 0
          oti%end_index     = 0
       end if
    else
       n_per_thread = nnodes / nthreads
       nlft = nnodes - n_per_thread*nthreads
       if (tid < nlft) then
          oti%begin_index = tid*(n_per_thread+1)
          oti%end_index   = oti%begin_index + n_per_thread+1
       else
          oti%begin_index = tid*n_per_thread + nlft
          oti%end_index   = oti%begin_index + n_per_thread
       end if
    end if
    oti%begin_index   = oti%begin_index + 1  ! because this is Fortran
    oti%end_index     = oti%end_index   + 1
    oti%current_index = oti%begin_index - 1
  end subroutine amrex_octree_iter_build

  subroutine amrex_octree_iter_destroy (oti)
    type(amrex_octree_iter) :: oti
  end subroutine amrex_octree_iter_destroy

  subroutine amrex_octree_iter_clear (oti)
    class(amrex_octree_iter) :: oti
    oti%current_index = oti%begin_index - 1 ! as when just built
  end subroutine amrex_octree_iter_clear

  logical function amrex_octree_iter_next (this)
    class(amrex_octree_iter) :: this
    this%current_index = this%current_index+1
    amrex_octree_iter_next = this%current_index < this%end_index
  end function amrex_octree_iter_next

  integer function amrex_octree_iter_level (this)
    class(amrex_octree_iter) :: this
    amrex_octree_iter_level = leaf_nodes(this%current_index)%level
  end function amrex_octree_iter_level

  integer function amrex_octree_iter_grid_index (this)
    class(amrex_octree_iter) :: this
    amrex_octree_iter_grid_index = leaf_nodes(this%current_index)%grid
  end function amrex_octree_iter_grid_index

  type(amrex_box) function amrex_octree_iter_box (this)
    class(amrex_octree_iter) :: this
    type(amrex_boxarray) :: ba
    ba = amrex_get_boxarray(leaf_nodes(this%current_index)%level)
    amrex_octree_iter_box = ba%get_box(leaf_nodes(this%current_index)%grid)
  end function amrex_octree_iter_box

  ! average down all components for all levels
  subroutine amrex_octree_average_down_leaves_all (mfs)
    type(amrex_multifab), intent(inout) :: mfs(0:)
    call amrex_octree_average_down_leaves_all_comp(mfs, 1, mfs(0)%ncomp())
  end subroutine amrex_octree_average_down_leaves_all

  ! average down ncomp components starting from component scomp for all levels
  subroutine amrex_octree_average_down_leaves_all_comp(mfs, scomp, ncomp)
    type(amrex_multifab), intent(inout) :: mfs(0:)
    integer, intent(in) :: scomp, ncomp
    integer :: ilev
    do ilev = amrex_max_level, 1, -1
       call amrex_octree_average_down_leaves_level_comp(ilev, &
            mfs(ilev), mfs(ilev-1), scomp, ncomp)
    end do
  end subroutine amrex_octree_average_down_leaves_all_comp

  ! average down flev to flev-1 for all component
  subroutine amrex_octree_average_down_leaves_level (flev, fine, coarse)
    integer, intent(in) :: flev
    type(amrex_multifab), intent(in) :: fine
    type(amrex_multifab), intent(inout) :: coarse
    call amrex_octree_average_down_leaves_level_comp(flev, &
         fine, coarse, 1, fine%ncomp())
  end subroutine amrex_octree_average_down_leaves_level

  ! average down flev to flev-1 for components [scomp, scomp+ncomp-1]
  subroutine amrex_octree_average_down_leaves_level_comp &
       (flev, fine, coarse, scomp, ncomp)
    integer, intent(in) :: flev, scomp, ncomp
    type(amrex_multifab), intent(in) :: fine
    type(amrex_multifab), intent(inout) :: coarse
    type(c_ptr) :: amrcore
    amrcore = amrex_get_amrcore()
    call amrex_fi_octree_average_down_level &
         (amrcore, flev, fine%p, coarse%p, scomp-1, ncomp) ! -1 for C++
  end subroutine amrex_octree_average_down_leaves_level_comp

end module amrex_octree_module
