
module amrex_multifab_module

  use iso_c_binding
  use amrex_error_module
  use amrex_fort_module, only : ndims => amrex_spacedim, amrex_real
  use amrex_box_module
  use amrex_boxarray_module
  use amrex_distromap_module
  use amrex_geometry_module
  use amrex_string_module
  use amrex_omp_module

  implicit none

  private

  public :: amrex_multifab_build, amrex_multifab_swap, amrex_multifab_destroy, amrex_multifab_write
  public :: amrex_imultifab_build_owner_mask
  public :: amrex_imultifab_build, amrex_imultifab_destroy
  public :: amrex_mfiter_build, amrex_mfiter_destroy

  type, public   :: amrex_multifab
     logical               :: owner = .false.
     type   (c_ptr)        :: p     =  c_null_ptr
     integer(c_int)        :: nc    =  0
     integer(c_int)        :: ng    =  0
     type(amrex_boxarray)  :: ba
     type(amrex_distromap) :: dm
   contains
     generic   :: assignment(=) => amrex_multifab_assign, amrex_multifab_install  ! shallow copy
     procedure :: move          => amrex_multifab_move     ! transfer ownership
     procedure :: ncomp         => amrex_multifab_ncomp
     procedure :: nghost        => amrex_multifab_nghost
     generic   :: dataPtr       => amrex_multifab_dataptr_iter, amrex_multifab_dataptr_int
     procedure :: min           => amrex_multifab_min
     procedure :: max           => amrex_multifab_max
     procedure :: norm0         => amrex_multifab_norm0
     procedure :: norm1         => amrex_multifab_norm1
     procedure :: norm2         => amrex_multifab_norm2
     procedure :: setval        => amrex_multifab_setval
     procedure :: copy          => amrex_multifab_copy     ! This copies the data
     generic   :: parallel_copy => amrex_multifab_parallel_copy, amrex_multifab_parallel_copy_c, &
          amrex_multifab_parallel_copy_cg
     generic   :: fill_boundary => amrex_multifab_fill_boundary, amrex_multifab_fill_boundary_c
     generic   :: override_sync => amrex_multifab_override_sync, amrex_multifab_override_sync_mask
     generic   :: sum_boundary  => amrex_multifab_sum_boundary, amrex_multifab_sum_boundary_c
     procedure :: average_sync  => amrex_multifab_average_sync
     procedure, private :: amrex_multifab_fill_boundary
     procedure, private :: amrex_multifab_fill_boundary_c
     procedure, private :: amrex_multifab_parallel_copy
     procedure, private :: amrex_multifab_parallel_copy_c
     procedure, private :: amrex_multifab_parallel_copy_cg
     procedure, private :: amrex_multifab_assign
     procedure, private :: amrex_multifab_install
     procedure, private :: amrex_multifab_dataptr_iter
     procedure, private :: amrex_multifab_dataptr_int
     procedure, private :: amrex_multifab_override_sync
     procedure, private :: amrex_multifab_override_sync_mask
     procedure, private :: amrex_multifab_sum_boundary
     procedure, private :: amrex_multifab_sum_boundary_c
     procedure, private :: amrex_multifab_average_sync
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_multifab_destroy
#endif
  end type amrex_multifab

  type, public   :: amrex_imultifab
     logical               :: owner = .false.
     type   (c_ptr)        :: p     =  c_null_ptr
     integer(c_int)        :: nc    =  0
     integer(c_int)        :: ng    =  0
     type(amrex_boxarray)  :: ba
     type(amrex_distromap) :: dm
   contains
     generic   :: assignment(=) => amrex_imultifab_assign   ! shallow copy
     procedure :: ncomp         => amrex_imultifab_ncomp
     procedure :: nghost        => amrex_imultifab_nghost
     procedure :: dataPtr       => amrex_imultifab_dataptr
     procedure :: setval        => amrex_imultifab_setval
     procedure, private :: amrex_imultifab_assign
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_imultifab_destroy
#endif
  end type amrex_imultifab

  type, public :: amrex_mfiter
     type(c_ptr)      :: p       = c_null_ptr
     integer ,private :: counter = -1 
   contains
     generic   :: assignment(=)    => amrex_mfiter_assign  ! will abort if called
     procedure :: clear            => amrex_mfiter_clear
     procedure :: next             => amrex_mfiter_next
     procedure :: grid_index       => amrex_mfiter_grid_index
     procedure :: tilebox          => amrex_mfiter_tilebox
     procedure :: nodaltilebox     => amrex_mfiter_nodaltilebox
     procedure :: fabbox           => amrex_mfiter_fabbox
     procedure, private :: amrex_mfiter_assign
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_mfiter_destroy
#endif
  end type amrex_mfiter

  interface amrex_mfiter_build
     module procedure amrex_mfiter_build_r
     module procedure amrex_mfiter_build_i
  end interface amrex_mfiter_build

  ! interfaces to c++ functions

  interface 
     subroutine amrex_fi_new_multifab (mf,ba,dm,nc,ng,nodal) bind(c)
       import
       implicit none
       type(c_ptr) :: mf, ba, dm
       integer(c_int), value :: nc, ng
       integer(c_int), intent(in) :: nodal(3)
     end subroutine amrex_fi_new_multifab
     
     subroutine amrex_fi_delete_multifab (mf) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf
     end subroutine amrex_fi_delete_multifab

     integer(c_int) function amrex_fi_multifab_ncomp (mf) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf
     end function amrex_fi_multifab_ncomp

     integer(c_int) function amrex_fi_multifab_ngrow (mf) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf
     end function amrex_fi_multifab_ngrow

     type(c_ptr) function amrex_fi_multifab_boxarray (mf) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf
     end function amrex_fi_multifab_boxarray

     type(c_ptr) function amrex_fi_multifab_distromap (mf) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf
     end function amrex_fi_multifab_distromap

     subroutine amrex_fi_multifab_dataptr_iter (mf, mfi, dp, lo, hi) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf, mfi
       type(c_ptr) :: dp
       integer(c_int) :: lo(3), hi(3)
     end subroutine amrex_fi_multifab_dataptr_iter

     subroutine amrex_fi_multifab_dataptr_int (mf, igrd, dp, lo, hi) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf
       integer, value :: igrd
       type(c_ptr) :: dp
       integer(c_int) :: lo(3), hi(3)
     end subroutine amrex_fi_multifab_dataptr_int

     function amrex_fi_multifab_min (mf, comp, nghost) bind(c)
       import
       implicit none
       real(amrex_real) :: amrex_fi_multifab_min
       type(c_ptr), value :: mf
       integer(c_int), value :: comp, nghost
     end function amrex_fi_multifab_min

     function amrex_fi_multifab_max (mf, comp, nghost) bind(c)
       import
       implicit none
       real(amrex_real) :: amrex_fi_multifab_max
       type(c_ptr), value :: mf
       integer(c_int), value :: comp, nghost
     end function amrex_fi_multifab_max

     function amrex_fi_multifab_norm0 (mf, comp) bind(c)
       import
       implicit none
       real(amrex_real) :: amrex_fi_multifab_norm0
       type(c_ptr), value :: mf
       integer(c_int), value :: comp
     end function amrex_fi_multifab_norm0

     function amrex_fi_multifab_norm1 (mf, comp) bind(c)
       import
       implicit none
       real(amrex_real) :: amrex_fi_multifab_norm1
       type(c_ptr), value :: mf
       integer(c_int), value :: comp
     end function amrex_fi_multifab_norm1

     function amrex_fi_multifab_norm2 (mf, comp) bind(c)
       import
       implicit none
       real(amrex_real) :: amrex_fi_multifab_norm2
       type(c_ptr), value :: mf
       integer(c_int), value :: comp
     end function amrex_fi_multifab_norm2

     subroutine amrex_fi_multifab_setval (mf, val, ic, nc, ng) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf
       real(amrex_real), value :: val
       integer(c_int), value :: ic, nc, ng
     end subroutine amrex_fi_multifab_setval

     subroutine amrex_fi_multifab_copy (dstmf, srcmf, srccomp, dstcomp, nc, ng) bind(c)
       import
       implicit none
       type(c_ptr), value :: dstmf, srcmf
       integer(c_int), value :: srccomp, dstcomp, nc, ng
     end subroutine amrex_fi_multifab_copy

     subroutine amrex_fi_multifab_parallelcopy(dstmf, srcmf, srccomp, dstcomp, nc,&
          srcng, dstng, geom) bind(c)
       import
       implicit none
       type(c_ptr), value :: dstmf, srcmf, geom
       integer(c_int), value :: srccomp, dstcomp, nc, srcng, dstng
     end subroutine amrex_fi_multifab_parallelcopy

     subroutine amrex_fi_multifab_fill_boundary (mf, geom, c, nc, cross) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf, geom
       integer(c_int), value :: c, nc, cross
     end subroutine amrex_fi_multifab_fill_boundary

     subroutine amrex_fi_write_multifab (mf, name) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf
       character(kind=c_char), intent(in) :: name(*)
     end subroutine amrex_fi_write_multifab

     subroutine amrex_fi_build_owner_imultifab (msk, ba, dm, data, geom) bind(c)
       import
       implicit none
       type(c_ptr) :: msk, ba, dm
       type(c_ptr), value :: data, geom
     end subroutine amrex_fi_build_owner_imultifab

     subroutine amrex_fi_multifab_override_sync (mf, geom) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf, geom
     end subroutine amrex_fi_multifab_override_sync

     subroutine amrex_fi_multifab_override_sync_mask (mf, geom, msk) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf, geom, msk
     end subroutine amrex_fi_multifab_override_sync_mask

     subroutine amrex_fi_multifab_sum_boundary (mf, geom, icomp, ncomp) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf, geom
       integer, value :: icomp, ncomp
     end subroutine amrex_fi_multifab_sum_boundary

     subroutine amrex_fi_multifab_average_sync (mf, geom) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf, geom
     end subroutine amrex_fi_multifab_average_sync
  end interface

  interface
     subroutine amrex_fi_new_imultifab (imf,ba,dm,nc,ng,nodal) bind(c)
       import
       implicit none
       type(c_ptr) :: imf, ba, dm
       integer(c_int), value :: nc, ng
       integer(c_int), intent(in) :: nodal(3)
     end subroutine amrex_fi_new_imultifab
     
     subroutine amrex_fi_delete_imultifab (imf) bind(c)
       import
       implicit none
       type(c_ptr), value :: imf
     end subroutine amrex_fi_delete_imultifab

     subroutine amrex_fi_imultifab_dataptr (imf, mfi, dp, lo, hi) bind(c)
       import
       implicit none
       type(c_ptr), value :: imf, mfi
       type(c_ptr) :: dp
       integer(c_int) :: lo(3), hi(3)
     end subroutine amrex_fi_imultifab_dataptr

     subroutine amrex_fi_imultifab_setval (imf, val, ic, nc, ng) bind(c)
       import
       implicit none
       type(c_ptr), value :: imf
       integer(c_int), value :: val
       integer(c_int), value :: ic, nc, ng
     end subroutine amrex_fi_imultifab_setval
  end interface

  interface
     subroutine amrex_fi_new_mfiter_r (mfi, mf, tiling) bind(c)
       import
       implicit none
       type(c_ptr) :: mfi
       type(c_ptr), value :: mf
       integer(c_int), value :: tiling
     end subroutine amrex_fi_new_mfiter_r

     subroutine amrex_fi_new_mfiter_i (mfi, imf, tiling) bind(c)
       import
       implicit none
       type(c_ptr) :: mfi
       type(c_ptr), value :: imf
       integer(c_int), value :: tiling
     end subroutine amrex_fi_new_mfiter_i

     subroutine amrex_fi_delete_mfiter (p) bind(c)
       import
       implicit none
       type(c_ptr), value :: p
     end subroutine amrex_fi_delete_mfiter

     subroutine amrex_fi_increment_mfiter (p, iv) bind(c)
       import
       implicit none
       type(c_ptr), value :: p
       integer(c_int) :: iv
     end subroutine amrex_fi_increment_mfiter

     subroutine amrex_fi_mfiter_is_valid (p, iv) bind(c)
       import
       implicit none
       type(c_ptr), value :: p
       integer(c_int) :: iv
     end subroutine amrex_fi_mfiter_is_valid

     integer(c_int) function amrex_fi_mfiter_grid_index(p) bind(c)
       import
       implicit none
       type(c_ptr), value :: p
     end function amrex_fi_mfiter_grid_index

     subroutine amrex_fi_mfiter_tilebox (p, lo, hi, nodal) bind(c)
       import
       implicit none
       type(c_ptr), value :: p
       integer(c_int) :: lo(3), hi(3), nodal(3)
     end subroutine amrex_fi_mfiter_tilebox

     subroutine amrex_fi_mfiter_nodaltilebox (p, dir, lo, hi, nodal) bind(c)
       import
       implicit none
       type(c_ptr), value :: p
       integer(c_int), value :: dir
       integer(c_int) :: lo(3), hi(3), nodal(3)
     end subroutine amrex_fi_mfiter_nodaltilebox

     subroutine amrex_fi_mfiter_fabbox (p, lo, hi, nodal) bind(c)
       import
       implicit none
       type(c_ptr), value :: p
       integer(c_int) :: lo(3), hi(3), nodal(3)
     end subroutine amrex_fi_mfiter_fabbox
  end interface

contains

  subroutine amrex_multifab_build (mf, ba, dm, nc, ng, nodal)
    type(amrex_multifab), intent(inout) :: mf
    type(amrex_boxarray), intent(in )   :: ba
    type(amrex_distromap),intent(in )   :: dm
    integer, intent(in) :: nc, ng
    logical, intent(in), optional :: nodal(*)
    integer :: inodal(3), dir
    mf%owner = .true.
    mf%nc = nc
    mf%ng = ng
    inodal = 0 
    if (present(nodal)) then
       do dir = 1, ndims
          if (nodal(dir)) inodal(dir) = 1
       end do
    end if
    mf%ba = ba
    mf%dm = dm
    call amrex_fi_new_multifab(mf%p, mf%ba%p, mf%dm%p, mf%nc, mf%ng, inodal)
  end subroutine amrex_multifab_build

  impure elemental subroutine amrex_multifab_destroy (this)
    type(amrex_multifab), intent(inout) :: this
    if (this%owner) then
       if (c_associated(this%p)) then
          call amrex_fi_delete_multifab(this%p)
       end if
    end if
    this%owner = .false.
    this%p = c_null_ptr
    call amrex_boxarray_destroy(this%ba)
    call amrex_distromap_destroy(this%dm)
  end subroutine amrex_multifab_destroy

  subroutine amrex_multifab_assign (dst, src)
    class(amrex_multifab), intent(inout) :: dst
    type (amrex_multifab), intent(in   ) :: src
    call amrex_multifab_destroy(dst)
    dst%owner = .false.
    dst%p     = src%p
    dst%nc    = src%nc
    dst%ng    = src%ng
    dst%ba    = src%ba
    dst%dm    = src%dm
  end subroutine amrex_multifab_assign

  subroutine amrex_multifab_install (this, p)
    class(amrex_multifab), intent(inout) :: this
    type(c_ptr), intent(in) :: p
    this%owner = .false.
    this%p     = p
    this%nc    = amrex_fi_multifab_ncomp(p)
    this%ng    = amrex_fi_multifab_ngrow(p)
    this%ba    = amrex_fi_multifab_boxarray(p)
    this%dm    = amrex_fi_multifab_distromap(p)
  end subroutine amrex_multifab_install

  subroutine amrex_multifab_move (dst, src)
    class(amrex_multifab), intent(inout) :: dst
    type (amrex_multifab), intent(inout) :: src
    call amrex_multifab_destroy(dst)
    dst%owner = src%owner
    dst%p     = src%p
    dst%nc    = src%nc
    dst%ng    = src%ng
    call dst%ba%move(src%ba)
    call dst%dm%move(src%dm)
    src%owner = .false.
    src%p     = c_null_ptr
  end subroutine amrex_multifab_move

  subroutine amrex_multifab_swap(mf1, mf2)
    type(amrex_multifab), intent(inout) :: mf1, mf2
    type(amrex_multifab) :: mftmp
    call mftmp%move(mf1)
    call mf1%move(mf2)
    call mf2%move(mftmp)
    call amrex_multifab_destroy(mftmp)
  end subroutine amrex_multifab_swap

  pure integer function amrex_multifab_ncomp (this)
    class(amrex_multifab), intent(in) :: this
    amrex_multifab_ncomp = this%nc
  end function amrex_multifab_ncomp

  pure integer function amrex_multifab_nghost (this)
    class(amrex_multifab), intent(in) :: this
    amrex_multifab_nghost = this%ng
  end function amrex_multifab_nghost

  function amrex_multifab_dataPtr_iter (this, mfi) result(dp)
    class(amrex_multifab), intent(in) :: this
    type(amrex_mfiter), intent(in) :: mfi
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: dp
    type(c_ptr) :: cp
    real(amrex_real), contiguous, pointer :: fp(:,:,:,:)
    integer(c_int) :: n(4)
    type(amrex_box) :: bx
    call amrex_fi_multifab_dataptr_iter(this%p, mfi%p, cp, bx%lo, bx%hi)
    n(1:3) = bx%hi - bx%lo + 1
    n(4)   = this%ncomp()
    call c_f_pointer(cp, fp, shape=n)
    dp(bx%lo(1):,bx%lo(2):,bx%lo(3):,1:) => fp
  end function amrex_multifab_dataPtr_iter

  function amrex_multifab_dataPtr_int (this, gid) result(dp)
    class(amrex_multifab), intent(in) :: this
    integer, intent(in) :: gid
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: dp
    type(c_ptr) :: cp
    real(amrex_real), contiguous, pointer :: fp(:,:,:,:)
    integer(c_int) :: n(4)
    type(amrex_box) :: bx
    call amrex_fi_multifab_dataptr_int(this%p, gid, cp, bx%lo, bx%hi)
    n(1:3) = bx%hi - bx%lo + 1
    n(4)   = this%ncomp()
    call c_f_pointer(cp, fp, shape=n)
    dp(bx%lo(1):,bx%lo(2):,bx%lo(3):,1:) => fp
  end function amrex_multifab_dataPtr_int

  function amrex_multifab_min (this, comp, nghost) result(r)
    class(amrex_multifab), intent(in) :: this
    integer(c_int), intent(in) :: comp
    integer(c_int), intent(in), optional :: nghost
    real(amrex_real) :: r
    if (present(nghost)) then
       r = amrex_fi_multifab_min(this%p, comp-1, nghost)
    else
       r = amrex_fi_multifab_min(this%p, comp-1, 0)
    end if
  end function amrex_multifab_min

  function amrex_multifab_max (this, comp, nghost) result(r)
    class(amrex_multifab), intent(in) :: this
    integer(c_int), intent(in) :: comp
    integer(c_int), intent(in), optional :: nghost
    real(amrex_real) :: r
    if (present(nghost)) then
       r = amrex_fi_multifab_max(this%p, comp-1, nghost)
    else
       r = amrex_fi_multifab_max(this%p, comp-1, 0)
    end if
  end function amrex_multifab_max

  function amrex_multifab_norm0 (this, comp) result(r)
    class(amrex_multifab), intent(in) :: this
    integer(c_int), intent(in), optional :: comp
    real(amrex_real) :: r
    if (present(comp)) then
       r = amrex_fi_multifab_norm0(this%p, comp-1)
    else
       r = amrex_fi_multifab_norm0(this%p, 0)
    end if
  end function amrex_multifab_norm0

  function amrex_multifab_norm1 (this, comp) result(r)
    class(amrex_multifab), intent(in) :: this
    integer(c_int), intent(in), optional :: comp
    real(amrex_real) :: r
    if (present(comp)) then
       r = amrex_fi_multifab_norm1(this%p, comp-1)
    else
       r = amrex_fi_multifab_norm1(this%p, 0)
    end if
  end function amrex_multifab_norm1

  function amrex_multifab_norm2 (this, comp) result(r)
    class(amrex_multifab), intent(in) :: this
    integer(c_int), intent(in), optional :: comp
    real(amrex_real) :: r
    if (present(comp)) then
       r = amrex_fi_multifab_norm2(this%p, comp-1)
    else
       r = amrex_fi_multifab_norm2(this%p, 0)
    end if
  end function amrex_multifab_norm2

  subroutine amrex_multifab_setval (this, val, icomp, ncomp, nghost)
    class(amrex_multifab), intent(inout) :: this
    real(amrex_real), intent(in) :: val
    integer, intent(in), optional :: icomp, ncomp, nghost
    integer :: ic, nc, ng
    ic = 0;         if (present(icomp))  ic = icomp-1
    nc = this%nc;   if (present(ncomp))  nc = ncomp
    ng = this%ng;   if (present(nghost)) ng = nghost
    call amrex_fi_multifab_setval(this%p, val, ic, nc, ng)
  end subroutine amrex_multifab_setval

  subroutine amrex_multifab_copy (this, srcmf, srccomp, dstcomp, nc, ng)
    class(amrex_multifab) :: this
    type(amrex_multifab), intent(in) :: srcmf
    integer, intent(in) :: srccomp, dstcomp, nc, ng
    call amrex_fi_multifab_copy(this%p, srcmf%p, srccomp-1, dstcomp-1, nc, ng)
  end subroutine amrex_multifab_copy

  subroutine amrex_multifab_parallel_copy (this, srcmf, geom)
    class(amrex_multifab) :: this
    type(amrex_multifab), intent(in) :: srcmf
    type(amrex_geometry), intent(in) :: geom
    call amrex_fi_multifab_parallelcopy(this%p, srcmf%p, 0, 0, this%nc, 0, 0, geom%p)
  end subroutine amrex_multifab_parallel_copy

  subroutine amrex_multifab_parallel_copy_c (this, srcmf, srccomp, dstcomp, nc, geom)
    class(amrex_multifab) :: this
    type(amrex_multifab), intent(in) :: srcmf
    type(amrex_geometry), intent(in) :: geom
    integer, intent(in) :: srccomp, dstcomp, nc
    call amrex_fi_multifab_parallelcopy(this%p, srcmf%p, srccomp-1, dstcomp-1, nc, 0, 0, geom%p)
  end subroutine amrex_multifab_parallel_copy_c

  subroutine amrex_multifab_parallel_copy_cg (this, srcmf, srccomp, dstcomp, nc, srcng, dstng, geom)
    class(amrex_multifab) :: this
    type(amrex_multifab), intent(in) :: srcmf
    type(amrex_geometry), intent(in) :: geom
    integer, intent(in) :: srccomp, dstcomp, nc, srcng, dstng
    call amrex_fi_multifab_parallelcopy(this%p, srcmf%p, srccomp-1, dstcomp-1, nc, srcng, dstng, geom%p)
  end subroutine amrex_multifab_parallel_copy_cg

  subroutine amrex_multifab_fill_boundary (this, geom, cross)
    class(amrex_multifab) :: this
    type(amrex_geometry), intent(in) :: geom
    logical, intent(in), optional :: cross
    call this%amrex_multifab_fill_boundary_c(geom, 1, this%nc, cross)
  end subroutine amrex_multifab_fill_boundary

  subroutine amrex_multifab_fill_boundary_c (this, geom, c, nc, cross)
    class(amrex_multifab) :: this
    type(amrex_geometry), intent(in) :: geom
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
    call amrex_fi_multifab_fill_boundary(this%p, geom%p, c-1, nc, lcross)
  end subroutine amrex_multifab_fill_boundary_c

  subroutine amrex_multifab_write (mf, name)
    type(amrex_multifab), intent(in) :: mf
    character(*), intent(in) :: name
    call amrex_fi_write_multifab(mf%p, amrex_string_f_to_c(name))
  end subroutine amrex_multifab_write

  subroutine amrex_imultifab_build_owner_mask (msk, data, geom)
    type(amrex_imultifab), intent(inout) :: msk
    type(amrex_multifab), intent(in) :: data
    type(amrex_geometry), intent(in) :: geom
    call amrex_imultifab_destroy(msk)
    msk%owner = .true.
    msk%nc = 1
    msk%ng = 0
    call amrex_fi_build_owner_imultifab(msk%p, msk%ba%p, msk%dm%p, data%p, geom%p)
  end subroutine amrex_imultifab_build_owner_mask

  subroutine amrex_multifab_override_sync (this, geom)
    class(amrex_multifab) :: this
    type(amrex_geometry), intent(in) :: geom
    call amrex_fi_multifab_override_sync(this%p, geom%p)
  end subroutine amrex_multifab_override_sync

  subroutine amrex_multifab_override_sync_mask (this, geom, msk)
    class(amrex_multifab) :: this
    type(amrex_geometry), intent(in) :: geom
    type(amrex_imultifab), intent(in) :: msk
    call amrex_fi_multifab_override_sync_mask(this%p, geom%p, msk%p)
  end subroutine amrex_multifab_override_sync_mask

  subroutine amrex_multifab_sum_boundary (this, geom)
    class(amrex_multifab) :: this
    type(amrex_geometry), intent(in) :: geom
    call this%amrex_multifab_sum_boundary_c(geom, 1, this%nc)
  end subroutine amrex_multifab_sum_boundary

  subroutine amrex_multifab_sum_boundary_c (this, geom, c, nc)
    class(amrex_multifab) :: this
    type(amrex_geometry), intent(in) :: geom
    integer, intent(in) :: c, nc
    call amrex_fi_multifab_sum_boundary(this%p, geom%p, c-1, nc)
  end subroutine amrex_multifab_sum_boundary_c

  subroutine amrex_multifab_average_sync (this, geom)
    class(amrex_multifab) :: this
    type(amrex_geometry), intent(in) :: geom
    call amrex_fi_multifab_average_sync(this%p, geom%p)
  end subroutine amrex_multifab_average_sync

!------ imultifab routines ------!

  subroutine amrex_imultifab_build (imf, ba, dm, nc, ng, nodal)
    type(amrex_imultifab) :: imf
    type(amrex_boxarray), intent(in ) :: ba
    type(amrex_distromap), intent(in ) :: dm
    integer, intent(in) :: nc, ng
    logical, intent(in), optional :: nodal(*)
    integer :: inodal(3), dir
    imf%owner = .true.
    imf%nc = nc
    imf%ng = ng
    inodal = 0
    if (present(nodal)) then
       do dir = 1, ndims
          if (nodal(dir)) inodal(dir) = 1
       end do
    end if
    imf%ba = ba
    imf%dm = dm
    call amrex_fi_new_imultifab(imf%p, imf%ba%p, imf%dm%p, imf%nc, imf%ng, inodal)
  end subroutine amrex_imultifab_build

  impure elemental subroutine amrex_imultifab_destroy (this)
    type(amrex_imultifab), intent(inout) :: this
    if (this%owner) then
       if (c_associated(this%p)) then
          call amrex_fi_delete_imultifab(this%p)
       end if
    end if
    this%owner = .false.
    this%p = c_null_ptr
    call amrex_boxarray_destroy(this%ba)
    call amrex_distromap_destroy(this%dm)
  end subroutine amrex_imultifab_destroy

  subroutine amrex_imultifab_assign (dst, src)
    class(amrex_imultifab), intent(inout) :: dst
    type (amrex_imultifab), intent(in   ) :: src
    call amrex_imultifab_destroy(dst)
    dst%owner = .false.
    dst%p     = src%p
    dst%nc    = src%nc
    dst%ng    = src%ng
    dst%ba    = src%ba
    dst%dm    = src%dm
  end subroutine amrex_imultifab_assign

  pure integer function amrex_imultifab_ncomp (this)
    class(amrex_imultifab), intent(in) :: this
    amrex_imultifab_ncomp = this%nc
  end function amrex_imultifab_ncomp

  pure integer function amrex_imultifab_nghost (this)
    class(amrex_imultifab), intent(in) :: this
    amrex_imultifab_nghost = this%ng
  end function amrex_imultifab_nghost

  function amrex_imultifab_dataPtr (this, mfi) result(dp)
    class(amrex_imultifab) :: this
    type(amrex_mfiter), intent(in) :: mfi
    integer, contiguous, pointer, dimension(:,:,:,:) :: dp
    type(c_ptr) :: cp
    integer, contiguous, pointer :: fp(:,:,:,:)
    integer(c_int) :: n(4)
    type(amrex_box) :: bx
    call amrex_fi_imultifab_dataptr(this%p, mfi%p, cp, bx%lo, bx%hi)
    n(1:3) = bx%hi - bx%lo + 1
    n(4)   = this%ncomp()
    call c_f_pointer(cp, fp, shape=n)
    dp(bx%lo(1):,bx%lo(2):,bx%lo(3):,1:) => fp
  end function amrex_imultifab_dataPtr

  subroutine amrex_imultifab_setval (this, val, icomp, ncomp, nghost)
    class(amrex_imultifab), intent(inout) :: this
    integer, intent(in) :: val
    integer, intent(in), optional :: icomp, ncomp, nghost
    integer :: ic, nc, ng
    ic = 0;         if (present(icomp))  ic = icomp-1
    nc = this%nc;   if (present(ncomp))  nc = ncomp
    ng = this%ng;   if (present(nghost)) ng = nghost
    call amrex_fi_imultifab_setval(this%p, val, ic, nc, ng)
  end subroutine amrex_imultifab_setval

!------ MFIter routines ------!

  subroutine amrex_mfiter_build_r (mfi, mf, tiling)
    type(amrex_mfiter) :: mfi
    type(amrex_multifab), intent(in ) :: mf
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
    call amrex_fi_new_mfiter_r(mfi%p, mf%p, t)
  end subroutine amrex_mfiter_build_r

  subroutine amrex_mfiter_build_i (mfi, imf, tiling)
    type(amrex_mfiter) :: mfi
    type(amrex_imultifab), intent(in ) :: imf
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
    call amrex_fi_new_mfiter_i(mfi%p, imf%p, t)
  end subroutine amrex_mfiter_build_i

  subroutine amrex_mfiter_destroy (this)
    type(amrex_mfiter) :: this
    call this%clear()
  end subroutine amrex_mfiter_destroy

  subroutine amrex_mfiter_assign (dst, src)
    class(amrex_mfiter), intent(inout) :: dst
    type (amrex_mfiter), intent(in   ) :: src
    ! No way to disable it at compile time, so ...
    call amrex_abort("amrex_mfiter assignment is disabled")
  end subroutine amrex_mfiter_assign

  subroutine amrex_mfiter_clear (this)
    class(amrex_mfiter) :: this
    this%counter = -1
    if (c_associated(this%p)) then
       call amrex_fi_delete_mfiter(this%p)
    end if
    this%p = c_null_ptr
  end subroutine amrex_mfiter_clear

  logical function amrex_mfiter_next (this)
    class(amrex_mfiter) :: this
    integer(c_int) :: isvalid
    this%counter = this%counter + 1
    if (this%counter == 1) then
       call amrex_fi_mfiter_is_valid(this%p, isvalid)
    else 
       call amrex_fi_increment_mfiter(this%p, isvalid)
    end if
    if (isvalid .eq. 1) then
       amrex_mfiter_next = .true.
    else
       call this%clear()
       amrex_mfiter_next = .false.
    end if
  end function amrex_mfiter_next

  integer function amrex_mfiter_grid_index (this)
    class(amrex_mfiter) :: this
    amrex_mfiter_grid_index = amrex_fi_mfiter_grid_index(this%p)
  end function amrex_mfiter_grid_index

  function amrex_mfiter_tilebox (this) result (bx)
    class(amrex_mfiter), intent(in) :: this
    type(amrex_box) :: bx
    integer :: inodal(3)
    inodal = 0
    call amrex_fi_mfiter_tilebox(this%p, bx%lo, bx%hi, inodal)
    where (inodal .ne. 0) bx%nodal = .true.  ! note default is false
  end function amrex_mfiter_tilebox

  function amrex_mfiter_nodaltilebox (this,dir_in) result (bx)
    class(amrex_mfiter), intent(in) :: this
    type(amrex_box) :: bx
    integer, intent(in), optional :: dir_in
    integer :: dir, inodal(3)
    if (present(dir_in)) then
       dir = dir_in
    else
       dir = -1
    end if
    inodal = 0
    call amrex_fi_mfiter_nodaltilebox(this%p, dir-1, bx%lo, bx%hi, inodal)
    where (inodal .ne. 0) bx%nodal = .true.  ! note default is false
  end function amrex_mfiter_nodaltilebox

  function amrex_mfiter_fabbox (this) result (bx)
    class(amrex_mfiter), intent(in) :: this
    type(amrex_box) :: bx
    integer :: dir, inodal(3)
    inodal = 0
    call amrex_fi_mfiter_fabbox(this%p, bx%lo, bx%hi, inodal)
    where (inodal .ne. 0) bx%nodal = .true.  ! note default is false
  end function amrex_mfiter_fabbox

end module amrex_multifab_module

