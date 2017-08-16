
module amrex_box_module

  use iso_c_binding
  use amrex_fort_module, only : ndims => amrex_spacedim, amrex_real

  implicit none

  private

  type, public :: amrex_box
     integer, dimension(3) :: lo    = 1
     integer, dimension(3) :: hi    = 1
     logical, dimension(3) :: nodal = .false.
   contains
     procedure :: numpts      => amrex_box_numpts
     procedure :: nodalize    => amrex_box_nodalize
     procedure :: cellize     => amrex_box_cellize
     procedure :: convert     => amrex_box_convert
     generic   :: grow        => amrex_box_grow_s, amrex_box_grow_v
     generic   :: intersects  => amrex_box_intersects_box, amrex_box_intersects_fp
     generic   :: contains    => amrex_box_contains_box, amrex_box_contains_fp, amrex_box_contains_pt
     procedure, private :: amrex_box_grow_s
     procedure, private :: amrex_box_grow_v
     procedure, private :: amrex_box_intersects_box
     procedure, private :: amrex_box_intersects_fp
     procedure, private :: amrex_box_contains_box
     procedure, private :: amrex_box_contains_fp
     procedure, private :: amrex_box_contains_pt
  end type amrex_box

  public :: amrex_print

  interface amrex_box
     module procedure amrex_box_build
  end interface amrex_box

  interface amrex_print
     module procedure amrex_box_print
  end interface amrex_print

  interface
     subroutine amrex_fi_print_box(lo,hi,nodal) bind(c)
       import
       implicit none
       integer(c_int), intent(in) :: lo(*), hi(*), nodal(*)
     end subroutine amrex_fi_print_box
  end interface

contains

  function amrex_box_build (lo, hi, nodal) result(bx)
    integer, intent(in) :: lo(ndims), hi(ndims)
    logical, intent(in), optional :: nodal(ndims)
    type(amrex_box) :: bx
    bx%lo(1:ndims) = lo(1:ndims)
    bx%hi(1:ndims) = hi(1:ndims)
    if (present(nodal)) bx%nodal(1:ndims) = nodal(1:ndims)
  end function amrex_box_build

  subroutine amrex_box_print (bx)
    type(amrex_box), intent(in) :: bx
    integer :: inodal(3)
    inodal = 0
    where (bx%nodal) inodal = 1
    call amrex_fi_print_box(bx%lo, bx%hi, inodal)
  end subroutine amrex_box_print

  function amrex_box_numpts (this) result(npts)
    integer(c_long) :: npts
    class(amrex_box), intent(in) :: this
    npts = (int(this%hi(1),c_long)-int(this%lo(1),c_long)+1_c_long) &
         * (int(this%hi(2),c_long)-int(this%lo(2),c_long)+1_c_long) &
         * (int(this%hi(3),c_long)-int(this%lo(3),c_long)+1_c_long)
  end function amrex_box_numpts

  subroutine amrex_box_nodalize (this, dir)
    class(amrex_box), intent(inout) :: this
    integer, intent(in) :: dir
    if (.not.this%nodal(dir)) then
       this%hi(dir) = this%hi(dir) + 1
       this%nodal(dir) = .true.
    end if
  end subroutine amrex_box_nodalize

  subroutine amrex_box_cellize (this, dir)
    class(amrex_box), intent(inout) :: this
    integer, intent(in) :: dir
    if (this%nodal(dir)) then
       this%hi(dir) = this%hi(dir) - 1
       this%nodal(dir) = .false.
    end if
  end subroutine amrex_box_cellize

  subroutine amrex_box_convert (this, flag)
    class(amrex_box), intent(inout) :: this
    logical, intent(in) :: flag(ndims)
    integer :: dir
    do dir = 1, ndims
       if (flag(dir) .and. .not.this%nodal(dir)) then
          this%hi(dir) = this%hi(dir) + 1
          this%nodal(dir) = .true.
       else if (.not.flag(dir) .and. this%nodal(dir)) then
          this%hi(dir) = this%hi(dir) - 1
          this%nodal(dir) = .false.
       end if
    end do
  end subroutine amrex_box_convert

  subroutine amrex_box_grow_s (this, i)
    class(amrex_box), intent(inout) :: this
    integer, intent(in) :: i
    this%lo = this%lo - i
    this%hi = this%hi + i
  end subroutine amrex_box_grow_s

  subroutine amrex_box_grow_v (this, i)
    class(amrex_box), intent(inout) :: this
    integer, intent(in) :: i(ndims)
    this%lo(1:ndims) = this%lo(1:ndims) - i(1:ndims)
    this%hi(1:ndims) = this%hi(1:ndims) + i(1:ndims)
  end subroutine amrex_box_grow_v

  function amrex_box_intersects_box (this, bx) result(r)
    class(amrex_box), intent(in) :: this
    type(amrex_box), intent(in) :: bx
    logical :: r
    integer :: seclo(3), sechi(3)
    seclo = max(this%lo, bx%lo)
    sechi = min(this%hi, bx%hi)
    r = all(seclo(1:ndims) .le. sechi(1:ndims))
  end function amrex_box_intersects_box

  function amrex_box_intersects_fp (this, p) result(r)
    class(amrex_box), intent(in) :: this
    real(amrex_real), pointer, intent(in) :: p(:,:,:,:)
    logical :: r
    integer :: seclo(3), sechi(3), plo(4), phi(4)
    plo = lbound(p)
    phi = ubound(p)
    seclo = max(this%lo, plo(1:3))
    sechi = min(this%hi, phi(1:3))
    r = all(seclo(1:ndims) .le. sechi(1:ndims))    
  end function amrex_box_intersects_fp

  function amrex_box_contains_box (this, bx) result(r)
    class(amrex_box), intent(in) :: this
    type(amrex_box), intent(in) :: bx
    logical :: r
    r = all(this%lo(1:ndims) .le. bx%lo(1:ndims)) .and. all(this%hi(1:ndims) .ge. bx%hi(1:ndims))
  end function amrex_box_contains_box

  function amrex_box_contains_fp (this, p) result(r)
    class(amrex_box), intent(in) :: this
    real(amrex_real), pointer, intent(in) :: p(:,:,:,:)
    logical :: r
    integer :: plo(4), phi(4)
    plo = lbound(p)
    phi = ubound(p)
    r = all(this%lo(1:ndims) .le. plo(1:ndims)) .and. all(this%hi(1:ndims) .ge. phi(1:ndims))
  end function amrex_box_contains_fp

  function amrex_box_contains_pt (this, pt) result(r)
    class(amrex_box), intent(in) :: this
    integer, intent(in) :: pt(ndims)
    logical :: r
    r = all(this%lo(1:ndims) .le. pt) .and. all(this%hi(1:ndims) .ge. pt)
  end function amrex_box_contains_pt

end module amrex_box_module

