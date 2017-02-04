
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
  public :: amrex_mfiter_build, amrex_mfiter_destroy

  type, public   :: amrex_multifab
     logical               :: owner = .false.
     type   (c_ptr)        :: p     =  c_null_ptr
     integer(c_int)        :: nc    =  0
     integer(c_int)        :: ng    =  0
     type(amrex_boxarray)  :: ba
     type(amrex_distromap) :: dm
   contains
     generic   :: assignment(=) => amrex_multifab_assign   ! shallow copy
     procedure :: move          => amrex_multifab_move     ! transfer ownership
     procedure :: ncomp         => amrex_multifab_ncomp
     procedure :: nghost        => amrex_multifab_nghost
     procedure :: dataPtr       => amrex_multifab_dataptr
     procedure :: min           => amrex_multifab_min
     procedure :: max           => amrex_multifab_max
     procedure :: norm0         => amrex_multifab_norm0
     procedure :: norm1         => amrex_multifab_norm1
     procedure :: norm2         => amrex_multifab_norm2
     generic   :: fill_boundary => amrex_multifab_fill_boundary, amrex_multifab_fill_boundary_c
     procedure, private :: amrex_multifab_fill_boundary, amrex_multifab_fill_boundary_c, &
          amrex_multifab_assign
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_multifab_destroy
#endif
  end type amrex_multifab

  type, public :: amrex_mfiter
     type(c_ptr)      :: p       = c_null_ptr
     integer ,private :: counter = -1 
   contains
     generic   :: assignment(=) => amrex_mfiter_assign  ! will abort if called
     procedure :: clear   => amrex_mfiter_clear
     procedure :: next    => amrex_mfiter_next
     procedure :: tilebox => amrex_mfiter_tilebox
     procedure, private :: amrex_mfiter_assign
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_mfiter_destroy
#endif
  end type amrex_mfiter

  interface amrex_multifab_build
     module procedure amrex_multifab_build
  end interface amrex_multifab_build

  interface amrex_mfiter_build
     module procedure amrex_mfiter_build
  end interface amrex_mfiter_build

  ! interfaces to c++ functions

  interface 
     subroutine amrex_fi_new_multifab (mf,ba,dm,nc,ng,nodal) bind(c)
       import
       type(c_ptr) :: mf, ba, dm
       integer(c_int), value :: nc, ng
       integer(c_int), intent(in) :: nodal(3)
     end subroutine amrex_fi_new_multifab
     
     subroutine amrex_fi_delete_multifab (mf) bind(c)
       import
       type(c_ptr), value :: mf
     end subroutine amrex_fi_delete_multifab

     subroutine amrex_fi_multifab_dataptr (mf, mfi, dp, lo, hi) bind(c)
       import
       type(c_ptr), value :: mf, mfi
       type(c_ptr) :: dp
       integer(c_int) :: lo(3), hi(3)
     end subroutine amrex_fi_multifab_dataptr

     function amrex_fi_multifab_min (mf, comp, nghost) bind(c)
       import
       real(amrex_real) :: amrex_fi_multifab_min
       type(c_ptr), value :: mf
       integer(c_int), value :: comp, nghost
     end function amrex_fi_multifab_min

     function amrex_fi_multifab_max (mf, comp, nghost) bind(c)
       import
       real(amrex_real) :: amrex_fi_multifab_max
       type(c_ptr), value :: mf
       integer(c_int), value :: comp, nghost
     end function amrex_fi_multifab_max

     function amrex_fi_multifab_norm0 (mf, comp) bind(c)
       import
       real(amrex_real) :: amrex_fi_multifab_norm0
       type(c_ptr), value :: mf
       integer(c_int), value :: comp
     end function amrex_fi_multifab_norm0

     function amrex_fi_multifab_norm1 (mf, comp) bind(c)
       import
       real(amrex_real) :: amrex_fi_multifab_norm1
       type(c_ptr), value :: mf
       integer(c_int), value :: comp
     end function amrex_fi_multifab_norm1

     function amrex_fi_multifab_norm2 (mf, comp) bind(c)
       import
       real(amrex_real) :: amrex_fi_multifab_norm2
       type(c_ptr), value :: mf
       integer(c_int), value :: comp
     end function amrex_fi_multifab_norm2

     subroutine amrex_fi_multifab_fill_boundary (mf, geom, c, nc, cross) bind(c)
       import
       type(c_ptr), value :: mf, geom
       integer(c_int), value :: c, nc, cross
     end subroutine amrex_fi_multifab_fill_boundary

     subroutine amrex_fi_write_multifab (mf, name) bind(c)
       import
       type(c_ptr), value :: mf
       character(c_char), intent(in) :: name(*)
     end subroutine amrex_fi_write_multifab
  end interface

  interface
     subroutine amrex_fi_new_mfiter (mfi, mf, tiling) bind(c)
       import
       type(c_ptr) :: mfi
       type(c_ptr), value :: mf
       integer(c_int), value :: tiling
     end subroutine amrex_fi_new_mfiter

     subroutine amrex_fi_delete_mfiter (p) bind(c)
       import
       type(c_ptr), value :: p
     end subroutine amrex_fi_delete_mfiter

     subroutine amrex_fi_increment_mfiter (p, iv) bind(c)
       import
       type(c_ptr), value :: p
       integer(c_int) :: iv
     end subroutine amrex_fi_increment_mfiter

     subroutine amrex_fi_mfiter_is_valid (p, iv) bind(c)
       import
       type(c_ptr), value :: p
       integer(c_int) :: iv
     end subroutine amrex_fi_mfiter_is_valid

     subroutine amrex_fi_mfiter_tilebox (p, lo, hi) bind(c)
       import
       type(c_ptr), value :: p
       integer(c_int) :: lo(3), hi(3)
     end subroutine amrex_fi_mfiter_tilebox
  end interface

contains

  subroutine amrex_multifab_build (mf, ba, dm, nc, ng, nodal)
    type(amrex_multifab) :: mf
    type(amrex_boxarray), intent(in ) :: ba
    type(amrex_distromap), intent(in ) :: dm
    integer, intent(in) :: nc, ng
    logical, intent(in), optional :: nodal(:)
    integer :: lnodal(3)
    mf%owner = .true.
    mf%nc = nc
    mf%ng = ng
    lnodal = 0 
    if (present(nodal)) then
       where (nodal .eqv. .true.) lnodal = 1
    end if
    mf%ba = ba
    mf%dm = dm
    call amrex_fi_new_multifab(mf%p, mf%ba%p, mf%dm%p, mf%nc, mf%ng, lnodal)
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

  function amrex_multifab_dataPtr (this, mfi) result(dp)
    class(amrex_multifab) :: this
    type(amrex_mfiter), intent(in) :: mfi
    double precision, contiguous, pointer, dimension(:,:,:,:) :: dp
    type(c_ptr) :: cp
    double precision, contiguous, pointer :: fp(:,:,:,:)
    integer(c_int) :: n(4)
    type(amrex_box) :: bx
    call amrex_fi_multifab_dataptr(this%p, mfi%p, cp, bx%lo, bx%hi)
    n(1:3) = bx%hi - bx%lo + 1
    n(4)   = this%ncomp()
    call c_f_pointer(cp, fp, shape=n)
    dp(bx%lo(1):,bx%lo(2):,bx%lo(3):,1:) => fp
  end function amrex_multifab_dataPtr

  function amrex_multifab_min (this, comp, nghost) result(r)
    class(amrex_multifab), intent(in) :: this
    integer(c_int), intent(in) :: comp
    integer(c_int), intent(in), optional :: nghost
    double precision :: r
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
    double precision :: r
    if (present(nghost)) then
       r = amrex_fi_multifab_max(this%p, comp-1, nghost)
    else
       r = amrex_fi_multifab_max(this%p, comp-1, 0)
    end if
  end function amrex_multifab_max

  function amrex_multifab_norm0 (this, comp) result(r)
    class(amrex_multifab), intent(in) :: this
    integer(c_int), intent(in), optional :: comp
    double precision :: r
    if (present(comp)) then
       r = amrex_fi_multifab_norm0(this%p, comp-1)
    else
       r = amrex_fi_multifab_norm0(this%p, 0)
    end if
  end function amrex_multifab_norm0

  function amrex_multifab_norm1 (this, comp) result(r)
    class(amrex_multifab), intent(in) :: this
    integer(c_int), intent(in), optional :: comp
    double precision :: r
    if (present(comp)) then
       r = amrex_fi_multifab_norm1(this%p, comp-1)
    else
       r = amrex_fi_multifab_norm1(this%p, 0)
    end if
  end function amrex_multifab_norm1

  function amrex_multifab_norm2 (this, comp) result(r)
    class(amrex_multifab), intent(in) :: this
    integer(c_int), intent(in), optional :: comp
    double precision :: r
    if (present(comp)) then
       r = amrex_fi_multifab_norm2(this%p, comp-1)
    else
       r = amrex_fi_multifab_norm2(this%p, 0)
    end if
  end function amrex_multifab_norm2

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
    call amrex_fi_multifab_fill_boundary (this%p, geom%p, c-1, nc, lcross)
  end subroutine amrex_multifab_fill_boundary_c

  subroutine amrex_multifab_write (mf, name)
    type(amrex_multifab), intent(in) :: mf
    character(*), intent(in) :: name
    call amrex_fi_write_multifab(mf%p, amrex_string_f_to_c(name))
  end subroutine amrex_multifab_write

!------ MFIter routines ------!

  subroutine amrex_mfiter_build (mfi, mf, tiling)
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
    call amrex_fi_new_mfiter(mfi%p, mf%p, t)
  end subroutine amrex_mfiter_build

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

  type(amrex_box) function amrex_mfiter_tilebox (this)
    class(amrex_mfiter), intent(in) :: this
    call amrex_fi_mfiter_tilebox(this%p, amrex_mfiter_tilebox%lo, amrex_mfiter_tilebox%hi)
  end function amrex_mfiter_tilebox

end module amrex_multifab_module

