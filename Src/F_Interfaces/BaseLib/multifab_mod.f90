
module multifab_module

  use iso_c_binding
  use bl_space_module, only : ndims => bl_num_dims
  use box_module
  use boxarray_module

  implicit none

  private

  type, public   :: MultiFab
     integer(c_int) :: nc  =  0
     integer(c_int) :: ng  =  0
     type   (c_ptr) :: p   =  c_null_ptr
   contains
     procedure :: ncomp    => multifab_ncomp
     procedure :: nghost   => multifab_nghost
     procedure :: dataPtr
     final :: destroy_multifab
  end type MultiFab

  type, public :: MFIter
     integer ,private :: counter = -1 
     type(c_ptr)      :: p       = c_null_ptr
   contains
     procedure :: clear   => mfiter_clear
     procedure :: next    => mfiter_next
     procedure :: tilebox
     final :: destroy_mfiter
  end type MFIter

  interface MultiFab
     module procedure build_multifab
  end interface MultiFab

  interface MFIter
     module procedure build_mfiter
  end interface MFIter

  interface 
     subroutine fi_new_multifab (mf,ba,nc,ng) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), intent(out) :: mf
       type(c_ptr), intent(in), value :: ba 
       integer(c_int), intent(in), value :: nc, ng
     end subroutine fi_new_multifab
     
     subroutine fi_delete_multifab (p) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: p
     end subroutine fi_delete_multifab

     subroutine fi_multifab_dataptr (mf, mfi, dp, lo, hi) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: mf, mfi
       type(c_ptr), intent(out) :: dp
       integer(c_int), intent(inout) :: lo(3), hi(3)
     end subroutine fi_multifab_dataptr
  end interface

  interface
     subroutine fi_new_mfiter (mfi, mf, tiling) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), intent(out) :: mfi
       type(c_ptr), intent(in), value :: mf
       integer(c_int), intent(in), value :: tiling
     end subroutine fi_new_mfiter

     subroutine fi_delete_mfiter (p) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: p
     end subroutine fi_delete_mfiter

     subroutine fi_increment_mfiter (p, iv) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: p
       integer(c_int), intent(out) :: iv
     end subroutine fi_increment_mfiter

     subroutine fi_mfiter_is_valid (p, iv) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: p
       integer(c_int), intent(out) :: iv
     end subroutine fi_mfiter_is_valid

     subroutine fi_mfiter_tilebox (p, lo, hi) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: p
       integer(c_int), intent(inout) :: lo(3), hi(3)
     end subroutine fi_mfiter_tilebox
  end interface

contains

  function build_multifab (ba, nc, ng) result(mf)
    type(BoxArray), intent(in) :: ba
    integer, intent(in) :: nc, ng
    type(MultiFab) :: mf
    mf%nc = nc
    mf%ng = ng
    call fi_new_multifab(mf%p, ba%p, mf%nc, mf%ng)
  end function build_multifab

  subroutine destroy_multifab (this)
    type(MultiFab) :: this
    if (c_associated(this%p)) then
       call fi_delete_multifab(this%p)
       this%p = c_null_ptr
    end if
  end subroutine destroy_multifab

  pure integer function multifab_ncomp (this)
    class(MultiFab), intent(in) :: this
    multifab_ncomp = this%nc
  end function multifab_ncomp

  pure integer function multifab_nghost (this)
    class(MultiFab), intent(in) :: this
    multifab_nghost = this%ng
  end function multifab_nghost

!------ MFIter routines ------!

  function build_mfiter (mf, tiling) result(mfi)
    type(MultiFab), intent(in) :: mf
    logical, intent(in), optional :: tiling
    type(MFIter) :: mfi
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
  end function build_mfiter

  subroutine destroy_mfiter (this)
    type(MFIter), intent(inout) :: this
    call this%clear()
  end subroutine destroy_mfiter

  subroutine mfiter_clear (this)
    class(MFIter), intent(inout) :: this
    this%counter = -1
    if (c_associated(this%p)) then
       call fi_delete_mfiter(this%p)
       this%p = c_null_ptr
    end if
  end subroutine mfiter_clear

  logical function mfiter_next (this)
    class(MFIter), intent(inout) :: this
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

  type(Box) function tilebox (this)
    class(MFIter), intent(in) :: this
    call fi_mfiter_tilebox(this%p, tilebox%lo, tilebox%hi)
  end function tilebox

  function dataPtr (this, mfi) result(dp)
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
  end function dataPtr

end module multifab_module

