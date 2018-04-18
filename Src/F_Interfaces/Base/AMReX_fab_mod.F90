
module amrex_fab_module

  use iso_c_binding
  use amrex_fort_module
  use amrex_box_module
  use amrex_mempool_module

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
     generic   :: assignment(=) => amrex_fab_assign  ! shallow copy
     procedure :: dataptr       => amrex_fab_dataptr
     procedure :: resize        => amrex_fab_resize
     procedure :: norminf       => amrex_fab_norminf
     procedure, private :: amrex_fab_assign
     procedure, private :: amrex_fab_resize
     procedure, private :: amrex_fab_norminf
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_fab_destroy
#endif
  end type amrex_fab

contains

  subroutine amrex_fab_build (fab, bx, nc)
    type(amrex_fab), intent(inout) :: fab
    type(amrex_box), intent(in   ) :: bx
    integer, intent(in) :: nc
    call amrex_fab_destroy(fab)
    fab%bx    = bx
    fab%nc    = nc
    fab%owner = .true.
    call amrex_allocate(fab%fp, bx%lo, bx%hi, nc)
    fab%cp = c_loc(fab%fp(bx%lo(1),bx%lo(2),bx%lo(3),1))
  end subroutine amrex_fab_build

  impure elemental subroutine amrex_fab_destroy (fab)
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

  subroutine amrex_fab_assign (dst, src)
    class(amrex_fab), intent(inout) :: dst
    type(amrex_fab), intent(in) :: src
    call amrex_fab_destroy(dst)
    dst%bx    =  src%bx
    dst%nc    =  src%nc
    dst%owner =  .false.
    dst%cp    =  src%cp
    dst%fp    => src%fp
  end subroutine amrex_fab_assign

  function amrex_fab_dataptr (this) result(dp)
    class(amrex_fab), intent(in) :: this
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: dp
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: tmp
    call c_f_pointer(this%cp, tmp, shape=[this%bx%hi(1)-this%bx%lo(1)+1, &
         &                                this%bx%hi(2)-this%bx%lo(2)+1, &
         &                                this%bx%hi(3)-this%bx%lo(3)+1, &
         &                                this%nc])
    dp(this%bx%lo(1):,this%bx%lo(2):,this%bx%lo(3):,1:) => tmp(:,:,:,:)
  end function amrex_fab_dataptr

  subroutine amrex_fab_resize (this, bx, nc)
    class(amrex_fab), intent(inout) :: this
    type(amrex_box), intent(in) :: bx
    integer, intent(in) :: nc
    if (.not.associated(this%fp) .or. &
         &         bx%numpts()*int(     nc,c_long) &
         .gt. this%bx%numpts()*int(this%nc,c_long)) then
       call amrex_fab_build(this, bx, nc)
    else
       this%bx = bx
       this%nc = nc
    end if
  end subroutine amrex_fab_resize

  real(amrex_real) function amrex_fab_norminf (this, ic, nc) result(r)
    class(amrex_fab), intent(in) :: this
    integer, intent(in) :: ic, nc
    real(amrex_real), contiguous, pointer :: data(:,:,:,:)
    integer :: i,j,k,n
    data => this%dataptr()
    r = -1._amrex_real
    do n = ic, ic+nc-1
       do       k = this%bx%lo(3), this%bx%hi(3)
          do    j = this%bx%lo(2), this%bx%hi(2)
             do i = this%bx%lo(1), this%bx%hi(1)
                r = max(r, abs(data(i,j,k,n)))
             end do
          end do
       end do
    end do
  end function amrex_fab_norminf

end module amrex_fab_module

