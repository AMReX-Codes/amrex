
module amrex_fab_module

  use iso_c_binding
  use amrex_fort_module
  use amrex_box_module
  use mempool_module

  implicit none
  private

  public :: amrex_fab_build, amrex_fab_destroy

  type, public :: amrex_fab
     type(amrex_box) :: bx
     integer         :: nc = 0
     logical, private :: owner = .false.
     type(c_ptr), private :: cp = c_null_ptr
     real(amrex_real), private, pointer, dimension(:,:,:,:) :: fp => null()
   contains
!     generic   :: assignment(=) => amrex_fab_assign  ! shallow copy
!     procedure :: dataptr       => amrex_fab_dataptr
!     procedure :: resize        => amrex_fab_resize
!     procedure, private :: amrex_fab_assign, amrex_fab_resize
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_fab_destroy
#endif
  end type amrex_fab

contains

  subroutine amrex_fab_build (fab, bx, ncomp)
    type(amrex_fab), intent(inout) :: fab
    type(amrex_box), intent(in   ) :: bx
    integer, intent(in) :: ncomp
    call amrex_fab_destroy(fab)
    fab%owner = .true.
    call amrex_allocate(fab%fp, bx%lo, bx%hi, ncomp)
    fab%cp = c_loc(fab%fp(bx%lo(1),bx%lo(2),bx%lo(3),1))
  end subroutine amrex_fab_build

  subroutine amrex_fab_destroy (fab)
    type(amrex_fab), intent(inout) :: fab
    if (fab%owner) then
       if (associated(fab%fp)) then
          call amrex_deallocate(fab%fp)
       end if
    end if
    fab%owner = .false.
    fab%cp    = c_null_ptr
    fab%fp    => null()
  end subroutine amrex_fab_destroy

end module amrex_fab_module

