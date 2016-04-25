
module multifab_module

  use iso_c_binding
  use bl_space_module, only : ndims => bl_num_dims
  use box_module
  use boxarray_module
  use geometry_module
  use omp_module

  implicit none

  private

  public :: multifab_build, multifab_swap, mfiter_build

  type, public   :: MultiFab
     integer(c_int) :: nc  =  0
     integer(c_int) :: ng  =  0
     type(BoxArray) :: ba
     type   (c_ptr) :: p   =  c_null_ptr
   contains
     procedure :: ncomp         => multifab_ncomp
     procedure :: nghost        => multifab_nghost
     procedure :: dataPtr       => multifab_dataptr
     procedure :: min           => multifab_min
     procedure :: max           => multifab_max
     procedure :: norm0         => multifab_norm0
     procedure :: norm1         => multifab_norm1
     procedure :: norm2         => multifab_norm2
     generic   :: fill_boundary => multifab_fill_boundary, multifab_fill_boundary_c
     procedure, private :: multifab_fill_boundary, multifab_fill_boundary_c
     final :: multifab_destroy
  end type MultiFab

  type, public :: MFIter
     integer ,private :: counter = -1 
     type(c_ptr)      :: p       = c_null_ptr
   contains
     procedure :: clear   => mfiter_clear
     procedure :: next    => mfiter_next
     procedure :: tilebox => mfiter_tilebox
     final :: mfiter_destroy
  end type MFIter

  interface multifab_build
     module procedure multifab_build
  end interface multifab_build

  interface mfiter_build
     module procedure mfiter_build
  end interface mfiter_build

  ! interfaces to c++ functions

  interface 
     subroutine fi_new_multifab (mf,bao,bai,nc,ng,nodal) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: mf, bao
       type(c_ptr), intent(in), value :: bai 
       integer(c_int), intent(in), value :: nc, ng
       integer(c_int), intent(in) :: nodal(3)
     end subroutine fi_new_multifab
     
     subroutine fi_delete_multifab (mf) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: mf
     end subroutine fi_delete_multifab

     subroutine fi_multifab_dataptr (mf, mfi, dp, lo, hi) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: mf, mfi
       type(c_ptr) :: dp
       integer(c_int) :: lo(3), hi(3)
     end subroutine fi_multifab_dataptr

     real(c_double) function fi_multifab_min (mf, comp, nghost) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: mf
       integer(c_int), value, intent(in) :: comp, nghost
     end function fi_multifab_min

     real(c_double) function fi_multifab_max (mf, comp, nghost) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: mf
       integer(c_int), value, intent(in) :: comp, nghost
     end function fi_multifab_max

     real(c_double) function fi_multifab_norm0 (mf, comp) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: mf
       integer(c_int), value, intent(in) :: comp
     end function fi_multifab_norm0

     real(c_double) function fi_multifab_norm1 (mf, comp) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: mf
       integer(c_int), value, intent(in) :: comp
     end function fi_multifab_norm1

     real(c_double) function fi_multifab_norm2 (mf, comp) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: mf
       integer(c_int), value, intent(in) :: comp
     end function fi_multifab_norm2

     subroutine fi_multifab_fill_boundary (mf, geom, c, nc, cross) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: mf, geom
       integer(c_int), value :: c, nc, cross
     end subroutine fi_multifab_fill_boundary
  end interface

  interface
     subroutine fi_new_mfiter (mfi, mf, tiling) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: mfi
       type(c_ptr), intent(in), value :: mf
       integer(c_int), intent(in), value :: tiling
     end subroutine fi_new_mfiter

     subroutine fi_delete_mfiter (p) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: p
     end subroutine fi_delete_mfiter

     subroutine fi_increment_mfiter (p, iv) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: p
       integer(c_int) :: iv
     end subroutine fi_increment_mfiter

     subroutine fi_mfiter_is_valid (p, iv) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: p
       integer(c_int) :: iv
     end subroutine fi_mfiter_is_valid

     subroutine fi_mfiter_tilebox (p, lo, hi) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: p
       integer(c_int) :: lo(3), hi(3)
     end subroutine fi_mfiter_tilebox
  end interface

contains

  subroutine multifab_build (mf, ba, nc, ng, nodal)
    type(MultiFab) :: mf
    type(BoxArray), intent(in ) :: ba
    integer, intent(in) :: nc, ng
    logical, intent(in), optional :: nodal(:)
    integer :: lnodal(3)
    mf%nc = nc
    mf%ng = ng
    lnodal = 0 
    if (present(nodal)) then
       where (nodal .eqv. .true.) lnodal = 1
    end if
    call fi_new_multifab(mf%p, mf%ba%p, ba%p, mf%nc, mf%ng, lnodal)
  end subroutine multifab_build

  subroutine multifab_destroy (this)
    type(MultiFab) :: this
    if (c_associated(this%p)) then
       call fi_delete_multifab(this%p)
       this%p = c_null_ptr
    end if
  end subroutine multifab_destroy

  subroutine multifab_swap(mf1, mf2)
    type(MultiFab), intent(inout) :: mf1, mf2
    integer :: itmp
    type(BoxArray) :: batmp
    type(c_ptr) :: ptmp
    itmp = mf1%nc;  mf1%nc = mf2%nc;  mf2%nc = itmp
    itmp = mf1%ng;  mf1%ng = mf2%ng;  mf2%ng = itmp
    batmp= mf1%ba;  mf1%ba = mf2%ba;  mf2%ba = batmp
    ptmp = mf1%p;   mf1%p  = mf2%p;   mf2%p  = ptmp
    ! the code below causes internal compiler error for gfortran 5.2
    ! tmp = mf1
    ! mf1 = mf2
    ! mf2 = tmp
    ! tmp%p = c_null_ptr
  end subroutine multifab_swap

  pure integer function multifab_ncomp (this)
    class(MultiFab), intent(in) :: this
    multifab_ncomp = this%nc
  end function multifab_ncomp

  pure integer function multifab_nghost (this)
    class(MultiFab), intent(in) :: this
    multifab_nghost = this%ng
  end function multifab_nghost

  function multifab_dataPtr (this, mfi) result(dp)
    class(MultiFab) :: this
    type(MFIter), intent(in) :: mfi
    double precision, pointer, dimension(:,:,:,:) :: dp
    type(c_ptr) :: cp
    double precision, pointer :: fp(:,:,:,:)
    integer(c_int) :: lo(3), hi(3), n(4)
    lo = 1;  hi = 1
    call fi_multifab_dataptr(this%p, mfi%p, cp, lo, hi)
    n(1:3) = hi - lo + 1
    n(4)   = this%ncomp()
    call c_f_pointer(cp, fp, shape=n)
    dp(lo(1):,lo(2):,lo(3):,1:) => fp
  end function multifab_dataPtr

  function multifab_min (this, comp, nghost) result(r)
    class(MultiFab), intent(in) :: this
    integer(c_int), intent(in) :: comp
    integer(c_int), intent(in), optional :: nghost
    double precision :: r
    if (present(nghost)) then
       r = fi_multifab_min(this%p, comp-1, nghost)
    else
       r = fi_multifab_min(this%p, comp-1, 0)
    end if
  end function multifab_min

  function multifab_max (this, comp, nghost) result(r)
    class(MultiFab), intent(in) :: this
    integer(c_int), intent(in) :: comp
    integer(c_int), intent(in), optional :: nghost
    double precision :: r
    if (present(nghost)) then
       r = fi_multifab_max(this%p, comp-1, nghost)
    else
       r = fi_multifab_max(this%p, comp-1, 0)
    end if
  end function multifab_max

  function multifab_norm0 (this, comp) result(r)
    class(MultiFab), intent(in) :: this
    integer(c_int), intent(in), optional :: comp
    double precision :: r
    if (present(comp)) then
       r = fi_multifab_norm0(this%p, comp-1)
    else
       r = fi_multifab_norm0(this%p, 0)
    end if
  end function multifab_norm0

  function multifab_norm1 (this, comp) result(r)
    class(MultiFab), intent(in) :: this
    integer(c_int), intent(in), optional :: comp
    double precision :: r
    if (present(comp)) then
       r = fi_multifab_norm1(this%p, comp-1)
    else
       r = fi_multifab_norm1(this%p, 0)
    end if
  end function multifab_norm1

  function multifab_norm2 (this, comp) result(r)
    class(MultiFab), intent(in) :: this
    integer(c_int), intent(in), optional :: comp
    double precision :: r
    if (present(comp)) then
       r = fi_multifab_norm2(this%p, comp-1)
    else
       r = fi_multifab_norm2(this%p, 0)
    end if
  end function multifab_norm2

  subroutine multifab_fill_boundary (this, geom, cross)
    class(MultiFab) :: this
    type(Geometry), intent(in) :: geom
    logical, intent(in), optional :: cross
    call this%multifab_fill_boundary_c(geom, 1, this%nc, cross)
  end subroutine multifab_fill_boundary

  subroutine multifab_fill_boundary_c (this, geom, c, nc, cross)
    class(MultiFab) :: this
    type(Geometry), intent(in) :: geom
    integer, intent(in) :: c, nc
    logical, intent(in), optional :: cross
    integer :: lcross
    lcross = 0  
    if (present(cross)) then
       if (cross) then
          lcross = 1
       else
          lcross = 0
       end if
    end if
    call fi_multifab_fill_boundary (this%p, geom%p, c-1, nc, lcross)
  end subroutine multifab_fill_boundary_c

!------ MFIter routines ------!

  subroutine mfiter_build (mfi, mf, tiling)
    type(MFIter) :: mfi
    type(MultiFab), intent(in ) :: mf
    logical, intent(in), optional :: tiling
    logical :: ltiling
    integer(c_int) :: t
    ltiling = .false.;  if (present(tiling)) ltiling = tiling
    if (ltiling) then
       t = 1
    else
       t = 0
    end if
    mfi%counter = 0
    call fi_new_mfiter(mfi%p, mf%p, t)
  end subroutine mfiter_build

  subroutine mfiter_destroy (this)
    type(MFIter) :: this
    call this%clear()
  end subroutine mfiter_destroy

  subroutine mfiter_clear (this)
    class(MFIter) :: this
    this%counter = -1
    if (c_associated(this%p)) then
       call fi_delete_mfiter(this%p)
       this%p = c_null_ptr
    end if
  end subroutine mfiter_clear

  logical function mfiter_next (this)
    class(MFIter) :: this
    integer(c_int) :: isvalid
    this%counter = this%counter + 1
    if (this%counter == 1) then
       call fi_mfiter_is_valid(this%p, isvalid)
    else 
       call fi_increment_mfiter(this%p, isvalid)
    end if
    if (isvalid .eq. 1) then
       mfiter_next = .true.
    else
       call this%clear()
       mfiter_next = .false.
    end if
  end function mfiter_next

  type(Box) function mfiter_tilebox (this)
    class(MFIter), intent(in) :: this
    call fi_mfiter_tilebox(this%p, mfiter_tilebox%lo, mfiter_tilebox%hi)
  end function mfiter_tilebox

end module multifab_module

