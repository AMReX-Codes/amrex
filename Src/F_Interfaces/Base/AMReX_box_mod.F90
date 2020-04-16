
module amrex_box_module

  use iso_c_binding
  use amrex_fort_module, only : ndims => amrex_spacedim, amrex_real, amrex_long, amrex_coarsen_intvect
  use amrex_error_module, only : amrex_error

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
     procedure :: refine      => amrex_box_refine
     procedure :: coarsen     => amrex_box_coarsen
     generic   :: grow        => amrex_box_grow_s, amrex_box_grow_v
     generic   :: intersects  => amrex_box_intersects_box, amrex_box_intersects_fp
     generic   :: contains    => amrex_box_contains_box, amrex_box_contains_fp, amrex_box_contains_pt
     procedure, private :: amrex_box_refine
     procedure, private :: amrex_box_coarsen
     procedure, private :: amrex_box_grow_s
     procedure, private :: amrex_box_grow_v
     procedure, private :: amrex_box_intersects_box
     procedure, private :: amrex_box_intersects_fp
     procedure, private :: amrex_box_contains_box
     procedure, private :: amrex_box_contains_fp
     procedure, private :: amrex_box_contains_pt
  end type amrex_box

  public :: amrex_print
  public :: amrex_intersection

  interface amrex_box
     module procedure amrex_box_build
  end interface amrex_box

  interface amrex_print
     module procedure amrex_box_print
  end interface amrex_print

  interface amrex_intersection
     module procedure amrex_box_intersection
  end interface amrex_intersection

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
    integer(amrex_long) :: npts
    class(amrex_box), intent(in) :: this
    npts = (int(this%hi(1),amrex_long)-int(this%lo(1),amrex_long)+1_amrex_long) &
         * (int(this%hi(2),amrex_long)-int(this%lo(2),amrex_long)+1_amrex_long) &
         * (int(this%hi(3),amrex_long)-int(this%lo(3),amrex_long)+1_amrex_long)
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

  subroutine amrex_box_refine (this, rr)
    class(amrex_box), intent(inout) :: this
    integer, intent(in) :: rr
    integer :: i
    this%lo(1:ndims) = this%lo(1:ndims) * rr
    do i = 1, ndims
       if (this%nodal(i)) then
          this%hi(i) = this%hi(i)*rr
       else
          this%hi(i) = (this%hi(i)+1)*rr-1
       end if
    end do
  end subroutine amrex_box_refine

  subroutine amrex_box_coarsen (this, rr)
    class(amrex_box), intent(inout) :: this
    integer, intent(in) :: rr
    integer :: i, off(ndims)
    this%lo(1:ndims) = amrex_coarsen_intvect(ndims, this%lo(1:ndims), rr)
    off = 0
    do i = 1, ndims
       if (this%nodal(i)) then
          if (mod(this%hi(i),rr) .ne. 0) then
             off(i) = 1
          end if
       end if
    end do
    this%hi(1:ndims) = amrex_coarsen_intvect(ndims, this%hi(1:ndims), rr) + off
  end subroutine amrex_box_coarsen

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

  function amrex_box_intersection(box1, box2) result(r)
    type(amrex_box), intent(in) :: box1
    type(amrex_box), intent(in) :: box2
    type(amrex_box) :: r

    if ( any(box1%nodal .neqv. box2%nodal) ) then
       call amrex_error("amrex_intersection: Boxes must be from same index space")
    end if

    r%nodal = box1%nodal
    r%lo(1:ndims) = max(box1%lo(1:ndims), box2%lo(1:ndims))
    r%hi(1:ndims) = min(box1%hi(1:ndims), box2%hi(1:ndims))
  end function amrex_box_intersection

end module amrex_box_module

