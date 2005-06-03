!! Fortran 90 style and boxes
module fab_module

  use bl_types
  use bl_mem_stat_module
  use box_module

  implicit none

  !! Each Fab type, integer, double, and logical, consists of
  !! DIM:    Dimension
  !  BX:     The Box in index space for which this FAB is defined
  !! IBX:    The index range of the valid data for the FAB
  !! PBX:    The physical box for the FAB
  !! NC:     Number of components
  !! NG:     Number of ghost cells
  
  !! When a FAB is created IBX = BX, unless it is nodal, in which case
  !! IBX = grow(BX, FACE=hi, 1 (in the nodal directions).
  !! PBX = grow(IBX, NG)
  
  !! For parallel systems the BX, IBX, PBX, etc are all defined, but the
  !! underlying pointer will not be allocated.

  !! All FABS are 'Four' dimensional, conventially, (NX,NY,NZ,NC) in size.
  !! NY = 1, NZ = 1, when DIM =1, NZ = 1, when DIM = 2.

  type fab
     integer   :: dim = 0
     type(box) :: bx
     type(box) :: pbx
     type(box) :: ibx
     integer   :: nc = 1
     integer   :: ng = 0
     real(kind = dp_t), pointer, dimension(:,:,:,:) :: p => Null()
  end type fab

  type ifab
     integer   :: dim = 0
     type(box) :: bx
     type(box) :: pbx
     type(box) :: ibx
     integer   :: nc = 1
     integer   :: ng = 0
     integer, pointer, dimension(:,:,:,:) :: p => Null()
  end type ifab

  type lfab
     integer   :: dim = 0
     type(box) :: bx
     type(box) :: pbx
     type(box) :: ibx
     integer   :: nc = 1
     integer   :: ng = 0
     logical, pointer, dimension(:,:,:,:) :: p => Null()
  end type lfab

  !! Returns the Dimension  of the FAB
!  interface get_dim
!     module procedure fab_dim
!     module procedure ifab_dim
!  end interface

  !! Returns whether the FAB has been built, this is different from
  !! whether the underlying pointer has been allocated since on parallel
  !! systems the pointer will not be allocated if the LAYOUT says it
  !! lives on another processor
  interface built_q
     module procedure fab_built_q
     module procedure ifab_built_q
  end interface

  interface conformant_q
     module procedure fab_conformant_q
     module procedure ifab_conformant_q
  end interface

  interface destroy
     module procedure fab_destroy
     module procedure ifab_destroy
     module procedure lfab_destroy
  end interface

  interface max_val
     module procedure fab_max
     module procedure fab_max_bx
     module procedure fab_max_c
     module procedure fab_max_bx_c
     module procedure ifab_max
     module procedure ifab_max_bx
     module procedure ifab_max_c
     module procedure ifab_max_bx_c
  end interface

  interface min_val
     module procedure fab_min
     module procedure fab_min_bx
     module procedure fab_min_c
     module procedure fab_min_bx_c
     module procedure ifab_min
     module procedure ifab_min_bx
     module procedure ifab_min_c
     module procedure ifab_min_bx_c
  end interface

  interface dataptr
     module procedure fab_dataptr
     module procedure fab_dataptr_bx
     module procedure fab_dataptr_c
     module procedure fab_dataptr_bx_c

     module procedure ifab_dataptr
     module procedure ifab_dataptr_bx
     module procedure ifab_dataptr_c
     module procedure ifab_dataptr_bx_c

     module procedure lfab_dataptr
     module procedure lfab_dataptr_bx
     module procedure lfab_dataptr_c
     module procedure lfab_dataptr_bx_c
  end interface

  interface setval
     module procedure fab_setval
     module procedure fab_setval_bx
     module procedure fab_setval_c
     module procedure fab_setval_bx_c

     module procedure ifab_setval
     module procedure ifab_setval_bx
     module procedure ifab_setval_c
     module procedure ifab_setval_bx_c

     module procedure lfab_setval
     module procedure lfab_setval_bx
     module procedure lfab_setval_c
     module procedure lfab_setval_bx_c
  end interface

  interface print
     module procedure fab_print
     module procedure ifab_print
     module procedure lfab_print
  end interface

!  interface lwb
!     module procedure fab_lwb
!     module procedure ifab_lwb
!  end interface

!   interface upb
!      module procedure fab_upb
!      module procedure ifab_upb
!   end interface

  interface get_box
     module procedure fab_get_box
     module procedure ifab_get_box
     module procedure lfab_get_box
  end interface

  interface get_pbox
     module procedure fab_get_pbox
     module procedure ifab_get_pbox
     module procedure lfab_get_pbox
  end interface

  interface get_ibox
     module procedure fab_get_ibox
     module procedure ifab_get_ibox
     module procedure lfab_get_ibox
  end interface

  interface volume
     module procedure fab_volume
     module procedure ifab_volume
     module procedure lfab_volume
  end interface

  interface ncomp
     module procedure fab_ncomp
     module procedure ifab_ncomp
  end interface

  real(kind=dp_t), parameter :: fab_default_init = -Huge(1.0_dp_t)
  integer, parameter :: ifab_default_init = -Huge(1)
  logical, parameter :: lfab_default_init = .FALSE.

  type(mem_stats), private, save ::  fab_ms
  type(mem_stats), private, save :: ifab_ms
  type(mem_stats), private, save :: lfab_ms

contains

  subroutine fab_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    fab_ms = ms
  end subroutine fab_set_mem_stats
  subroutine ifab_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    ifab_ms = ms
  end subroutine ifab_set_mem_stats
  subroutine lfab_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    lfab_ms = ms
  end subroutine lfab_set_mem_stats

  function fab_mem_stats() result(r)
    type(mem_stats) :: r
    r = fab_ms
  end function fab_mem_stats
  function ifab_mem_stats() result(r)
    type(mem_stats) :: r
    r = ifab_ms
  end function ifab_mem_stats
  function lfab_mem_stats() result(r)
    type(mem_stats) :: r
    r = lfab_ms
  end function lfab_mem_stats

  function fab_volume(fb, all) result(r)
    integer(kind=ll_t) :: r
    type(fab), intent(in) :: fb
    logical, optional :: all
    if ( all ) then
       r = volume(get_pbox(fb))
    else
       r = volume(get_box(fb))
    end if
  end function fab_volume
  function ifab_volume(fb, all) result(r)
    integer(kind=ll_t) :: r
    type(ifab), intent(in) :: fb
    logical, optional :: all
    if ( all ) then
       r = volume(get_pbox(fb))
    else
       r = volume(get_box(fb))
    end if
  end function ifab_volume
  function lfab_volume(fb, all) result(r)
    integer(kind=ll_t) :: r
    type(lfab), intent(in) :: fb
    logical, optional :: all
    if ( all ) then
       r = volume(get_pbox(fb))
    else
       r = volume(get_box(fb))
    end if
  end function lfab_volume

  function fab_lwb(fb) result(r)
    type(fab), intent(in) :: fb
    integer :: r(fb%dim)
    r = lwb(fb%bx)
  end function fab_lwb
  function fab_lwb_n(fb,dim) result(r)
    type(fab), intent(in) :: fb
    integer, intent(in) :: dim
    integer :: r
    r = lwb(fb%bx,dim)
  end function fab_lwb_n
  function ifab_lwb(fb) result(r)
    type(ifab), intent(in) :: fb
    integer :: r(fb%dim)
    r = lwb(fb%bx)
  end function ifab_lwb

  function fab_ilwb(fb) result(r)
    type(fab), intent(in) :: fb
    integer :: r(fb%dim)
    r = lwb(fb%ibx)
  end function fab_ilwb
  function fab_ilwb_n(fb,dim) result(r)
    type(fab), intent(in) :: fb
    integer, intent(in) :: dim
    integer :: r
    r = lwb(fb%ibx,dim)
  end function fab_ilwb_n
  function ifab_ilwb(fb) result(r)
    type(ifab), intent(in) :: fb
    integer :: r(fb%dim)
    r = lwb(fb%ibx)
  end function ifab_ilwb

  function fab_upb(fb) result(r)
    type(fab), intent(in) :: fb
    integer :: r(fb%dim)
    r = upb(fb%bx)
  end function fab_upb
  function fab_upb_n(fb,dim) result(r)
    type(fab), intent(in) :: fb
    integer, intent(in) :: dim
    integer :: r
    r = upb(fb%bx, dim)
  end function fab_upb_n
  function ifab_upb(fb) result(r)
    type(ifab), intent(in) :: fb
    integer :: r(fb%dim)
    r = upb(fb%bx)
  end function ifab_upb

  function fab_iupb(fb) result(r)
    type(fab), intent(in) :: fb
    integer :: r(fb%dim)
    r = upb(fb%ibx)
  end function fab_iupb
  function fab_iupb_n(fb,dim) result(r)
    type(fab), intent(in) :: fb
    integer, intent(in) :: dim
    integer :: r
    r = upb(fb%ibx, dim)
  end function fab_iupb_n
  function ifab_iupb(fb) result(r)
    type(ifab), intent(in) :: fb
    integer :: r(fb%dim)
    r = upb(fb%ibx)
  end function ifab_iupb

  function fab_plwb(fb) result(r)
    type(fab), intent(in) :: fb
    integer :: r(fb%dim)
    r = lwb(fb%pbx)
  end function fab_plwb
  function ifab_plwb(fb) result(r)
    type(ifab), intent(in) :: fb
    integer :: r(fb%dim)
    r = lwb(fb%pbx)
  end function ifab_plwb

  function fab_pupb(fb) result(r)
    type(fab), intent(in) :: fb
    integer :: r(fb%dim)
    r = upb(fb%pbx)
  end function fab_pupb
  function ifab_pupb(fb) result(r)
    type(ifab), intent(in) :: fb
    integer :: r(fb%dim)
    r = upb(fb%pbx)
  end function ifab_pupb

  function fab_ncomp(fb) result(r)
    integer :: r
    type(fab), intent(in) :: fb
    r = fb%nc
  end function fab_ncomp
  function ifab_ncomp(fb) result(r)
    integer :: r
    type(ifab), intent(in) :: fb
    r = fb%nc
  end function ifab_ncomp

  function fab_get_box(fb) result(r)
    type(fab), intent(in) :: fb
    type(box) :: r
    r = fb%bx
  end function fab_get_box
  function ifab_get_box(fb) result(r)
    type(ifab), intent(in) :: fb
    type(box) :: r
    r = fb%bx
  end function ifab_get_box
  function lfab_get_box(fb) result(r)
    type(lfab), intent(in) :: fb
    type(box) :: r
    r = fb%bx
  end function lfab_get_box

  function fab_get_pbox(fb) result(r)
    type(fab), intent(in) :: fb
    type(box) :: r
    r = fb%pbx
  end function fab_get_pbox
  function ifab_get_pbox(fb) result(r)
    type(ifab), intent(in) :: fb
    type(box) :: r
    r = fb%pbx
  end function ifab_get_pbox
  function lfab_get_pbox(fb) result(r)
    type(lfab), intent(in) :: fb
    type(box) :: r
    r = fb%pbx
  end function lfab_get_pbox

  function fab_get_ibox(fb) result(r)
    type(fab), intent(in) :: fb
    type(box) :: r
    r = fb%ibx
  end function fab_get_ibox
  function ifab_get_ibox(fb) result(r)
    type(ifab), intent(in) :: fb
    type(box) :: r
    r = fb%ibx
  end function ifab_get_ibox
  function lfab_get_ibox(fb) result(r)
    type(lfab), intent(in) :: fb
    type(box) :: r
    r = fb%ibx
  end function lfab_get_ibox
  
  function fab_built_q(fb) result(r)
    type(fab), intent(in) :: fb
    logical :: r
    r = fb%dim /= 0
  end function fab_built_q
  function ifab_built_q(fb) result(r)
    type(ifab), intent(in) :: fb
    logical :: r
    r = fb%dim /= 0
  end function ifab_built_q

  function fab_conformant_q(a,b) result(r)
    logical :: r
    type(fab), intent(in) :: a, b
    if ( .not. built_q(a) .AND. .not. built_q(b) ) then
       r = .true.
    else if ( built_q(a) .and. built_q(b) ) then
       r = a%dim == b%dim .AND. a%ng == b%ng .AND. &
            a%nc == b%nc .AND. a%ibx == b%ibx
    else
       r = .FALSE.
    end if 
  end function fab_conformant_q
  function ifab_conformant_q(a,b) result(r)
    logical :: r
    type(ifab), intent(in) :: a, b
    if ( .not. built_q(a) .AND. .not. built_q(b) ) then
       r = .true.
    else if ( built_q(a) .and. built_q(b) ) then
       r = a%dim == b%dim .AND. a%ng == b%ng .AND. &
            a%nc == b%nc .AND. a%ibx == b%ibx
    else
       r = .FALSE.
    end if 
  end function ifab_conformant_q

  subroutine fab_set_border_val(fb, val)
    type(fab), intent(inout) :: fb
    real(kind=dp_t), intent(in), optional :: val
    real(kind=dp_t) :: tval
    integer :: i,j,k,n
    integer :: plo(MAX_SPACEDIM), phi(MAX_SPACEDIM)
    integer ::  lo(MAX_SPACEDIM),  hi(MAX_SPACEDIM)
    tval = 0.0_dp_t; if ( present(val) ) tval = val
    plo = 1; plo(1:fb%dim) = lwb(fb%pbx)
    phi = 1; phi(1:fb%dim) = upb(fb%pbx)
    lo = 1;  lo(1:fb%dim) = lwb(fb%ibx)
    hi = 1;  hi(1:fb%dim) = upb(fb%ibx)
    do n = 1, fb%nc
       do k = plo(3), phi(3)
          do j = plo(2), phi(2)
             do i = plo(1), phi(1)
                 if ( k >= lo(3) .and. k <= hi(3) .and. &
                      j >= lo(2) .and. j <= hi(2) .and. &
                      i >= lo(1) .and. i <= hi(1) ) cycle
                 fb%p(i,j,k,n) = tval
             end do
          end do
       end do
    end do
  end subroutine fab_set_border_val
  subroutine ifab_set_border_val(fb, val)
    type(ifab), intent(inout) :: fb
    integer, intent(in), optional :: val
    integer :: tval
    integer :: i,j,k,n
    integer :: plo(MAX_SPACEDIM), phi(MAX_SPACEDIM)
    integer ::  lo(MAX_SPACEDIM),  hi(MAX_SPACEDIM)
    tval = 0; if ( present(val) ) tval = val
    plo = 1; plo(1:fb%dim) = lwb(fb%pbx)
    phi = 1; phi(1:fb%dim) = upb(fb%pbx)
    lo = 1;  lo(1:fb%dim) = lwb(fb%ibx)
    hi = 1;  hi(1:fb%dim) = upb(fb%ibx)
    do n = 1, fb%nc
       do k = plo(3), phi(3)
          do j = plo(2), phi(2)
             do i = plo(1), phi(1)
                 if ( k >= lo(3) .and. k <= hi(3) .and. &
                      j >= lo(2) .and. j <= hi(2) .and. &
                      i >= lo(1) .and. i <= hi(1) ) cycle
                 fb%p(i,j,k,n) = tval
             end do
          end do
       end do
    end do
  end subroutine ifab_set_border_val
  !! lfab_set_border_val:
  !! The default value is .false. since that is sort of mask
  !! like on 'non-valid' data.
  subroutine lfab_set_border_val(fb, val)
    type(lfab), intent(inout) :: fb
    logical, intent(in), optional :: val
    logical :: tval
    integer :: i,j,k,n
    integer :: plo(MAX_SPACEDIM), phi(MAX_SPACEDIM)
    integer ::  lo(MAX_SPACEDIM),  hi(MAX_SPACEDIM)
    tval = .false.; if ( present(val) ) tval = val
    plo = 1; plo(1:fb%dim) = lwb(fb%pbx)
    phi = 1; phi(1:fb%dim) = upb(fb%pbx)
    lo = 1;  lo(1:fb%dim) = lwb(fb%ibx)
    hi = 1;  hi(1:fb%dim) = upb(fb%ibx)
    do n = 1, fb%nc
       do k = plo(3), phi(3)
          do j = plo(2), phi(2)
             do i = plo(1), phi(1)
                 if ( k >= lo(3) .and. k <= hi(3) .and. &
                      j >= lo(2) .and. j <= hi(2) .and. &
                      i >= lo(1) .and. i <= hi(1) ) cycle
                 fb%p(i,j,k,n) = tval
             end do
          end do
       end do
    end do
  end subroutine lfab_set_border_val

  subroutine fab_build(fb, bx, nc, ng, nodal, alloc)
    type(fab), intent(out) :: fb
    type(box), intent(in)  :: bx
    integer, intent(in), optional :: ng, nc
    logical, intent(in), optional :: nodal(:)
    logical, intent(in), optional :: alloc
    integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
    integer :: lnc, lng
    logical :: lal
    lng = 0; if ( present(ng) ) lng = ng
    lnc = 1; if ( present(nc) ) lnc = nc
    lal = .true.; if ( present(alloc) ) lal = alloc
    lo = 1
    hi = 1
    fb%dim = bx%dim
    fb%bx = bx
    fb%ng = lng
    fb%nc = lnc
    fb%ibx = box_nodalize(bx, nodal)
    fb%pbx = grow(fb%ibx, lng)
    lo(1:fb%dim) = fb%pbx%lo(1:fb%dim)
    hi(1:fb%dim) = fb%pbx%hi(1:fb%dim)
    if ( lal ) then
       allocate(fb%p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:lnc))
       call setval(fb, fab_default_init)
    end if
    call mem_stats_alloc(fab_ms, volume(fb, all=.TRUE.))
  end subroutine fab_build

  subroutine ifab_build(fb, bx, nc, ng, nodal, alloc)
    type(ifab), intent(out) :: fb
    type(box), intent(in)  :: bx
    integer, intent(in), optional :: ng, nc
    logical, intent(in), optional :: nodal(:)
    logical, intent(in), optional :: alloc
    integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
    integer :: lnc, lng
    logical :: lal
    lng = 0; if ( present(ng) ) lng = ng
    lnc = 1; if ( present(nc) ) lnc = nc
    lal = .true.; if ( present(alloc) ) lal = alloc
    lo = 1
    hi = 1
    fb%dim = bx%dim
    fb%bx = bx
    fb%ng = lng
    fb%nc = lnc
    fb%ibx = box_nodalize(bx, nodal)
    fb%pbx = grow(fb%ibx, lng)
    lo(1:fb%dim) = fb%pbx%lo(1:fb%dim)
    hi(1:fb%dim) = fb%pbx%hi(1:fb%dim)
    if ( lal ) then
       allocate(fb%p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:lnc))
       call setval(fb, ifab_default_init)
    end if
    call mem_stats_alloc(ifab_ms, volume(fb, all=.TRUE.))
  end subroutine ifab_build

  subroutine lfab_build(fb, bx, nc, ng, nodal, alloc)
    type(lfab), intent(out) :: fb
    type(box), intent(in)  :: bx
    integer, intent(in), optional :: ng, nc
    logical, intent(in), optional :: nodal(:)
    logical, intent(in), optional :: alloc
    integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
    integer :: lnc, lng
    logical :: lal
    lng = 0; if ( present(ng) ) lng = ng
    lnc = 1; if ( present(nc) ) lnc = nc
    lal = .true.; if ( present(alloc) ) lal = alloc
    lo = 1
    hi = 1
    fb%dim = bx%dim
    fb%bx = bx
    fb%ng = lng
    fb%nc = lnc
    fb%ibx = box_nodalize(bx, nodal)
    fb%pbx = grow(fb%ibx, lng)
    lo(1:fb%dim) = fb%pbx%lo(1:fb%dim)
    hi(1:fb%dim) = fb%pbx%hi(1:fb%dim)
    if ( lal ) then
       allocate(fb%p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:lnc))
       call setval(fb, lfab_default_init)
    end if
    call mem_stats_alloc(lfab_ms, volume(fb, all=.TRUE.))
  end subroutine lfab_build

  subroutine fab_destroy(fb)
    type(fab), intent(inout) :: fb
    call mem_stats_dealloc(fab_ms, volume(fb, all=.TRUE.))
    if ( associated(fb%p) ) deallocate(fb%p)
    fb%bx  = nobox(fb%dim)
    fb%dim = 0
    fb%nc  = 0
    fb%ng  = 0
  end subroutine fab_destroy
  subroutine ifab_destroy(fb)
    type(ifab), intent(inout) :: fb
    call mem_stats_dealloc(ifab_ms, volume(fb, all=.TRUE.))
    if ( associated(fb%p) ) deallocate(fb%p)
    fb%bx  = nobox(fb%dim)
    fb%dim = 0
    fb%nc  = 0
    fb%ng  = 0
  end subroutine ifab_destroy
  subroutine lfab_destroy(fb)
    type(lfab), intent(inout) :: fb
    call mem_stats_dealloc(lfab_ms, volume(fb, all=.TRUE.))
    if ( associated(fb%p) ) deallocate(fb%p)
    fb%bx  = nobox(fb%dim)
    fb%dim = 0
    fb%nc  = 0
    fb%ng  = 0
  end subroutine lfab_destroy

  function fab_dim(fb) result(r)
    type(fab), intent(in) :: fb
    integer :: r
    r = box_dim(fb%bx)
  end function fab_dim
  function ifab_dim(fb) result(r)
    type(ifab), intent(in) :: fb
    integer :: r
    r = box_dim(fb%bx)
  end function ifab_dim

  function fab_dataptr(fb) result(r)
    type(fab), intent(in) :: fb
    real(kind=dp_t), pointer :: r(:,:,:,:)
    r => fb%p
  end function fab_dataptr
  function ifab_dataptr(fb) result(r)
    type(ifab), intent(in) :: fb
    integer, pointer :: r(:,:,:,:)
    r => fb%p
  end function ifab_dataptr
  function lfab_dataptr(fb) result(r)
    type(lfab), intent(in) :: fb
    logical, pointer :: r(:,:,:,:)
    r => fb%p
  end function lfab_dataptr

  function fab_dataptr_bx(fb, bx) result(r)
    type(fab), intent(in) :: fb
    type(box), intent(in) :: bx
    real(kind=dp_t), pointer :: r(:,:,:,:)
    if ( .not. contains(fb%pbx, bx) ) call bl_error('fab_dataptr_bx: bx is too large')
    select case (fb%dim)
    case (1)
       r => fb%p(bx%lo(1):bx%hi(1),:,:,:)
    case (2)
       r => fb%p(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),:,:)
    case (3)
       r => fb%p(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),:)
    end select
  end function fab_dataptr_bx
  function ifab_dataptr_bx(fb, bx) result(r)
    type(ifab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, pointer :: r(:,:,:,:)
    if ( .not. contains(fb%pbx, bx) ) call bl_error('ifab_dataptr_bx: bx is too large')
    select case (fb%dim)
    case (1)
       r => fb%p(bx%lo(1):bx%hi(1),:,:,:)
    case (2)
       r => fb%p(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),:,:)
    case (3)
       r => fb%p(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),:)
    end select
  end function ifab_dataptr_bx
  function lfab_dataptr_bx(fb, bx) result(r)
    type(lfab), intent(in) :: fb
    type(box), intent(in) :: bx
    logical, pointer :: r(:,:,:,:)
    if ( .not. contains(fb%pbx, bx) ) call bl_error('lfab_dataptr_bx: bx is too large')
    select case (fb%dim)
    case (1)
       r => fb%p(bx%lo(1):bx%hi(1),:,:,:)
    case (2)
       r => fb%p(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),:,:)
    case (3)
       r => fb%p(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),:)
    end select
  end function lfab_dataptr_bx

  subroutine fab_debug_fill(fb, ii, all, loc)
    type(fab), intent(inout) :: fb
    integer, intent(in) :: ii
    logical, intent(in), optional :: all, loc
    logical :: lall, lloc
    integer :: i, j, k, n
    real(kind=dp_t) :: xx,yy,zz,fc,nn
    real(kind=dp_t) :: sgn
    lall = .false.; if ( present(all) ) lall = all
    lloc = .false.; if ( present(loc) ) lloc = loc
    do n = lbound(fb%p,4), ubound(fb%p,4)
       if ( lall ) then
          sgn = 1
       else
          sgn = -1
       end if
       select case (fb%dim)
       case (3)
          nn = og(n)/100/100/100/100
          do k = lbound(fb%p,3), ubound(fb%p,3)
             zz = og(k)/100/100/100
             do j = lbound(fb%p,2), ubound(fb%p,2)
                yy = og(j)/100/100
                do i = lbound(fb%p,1), ubound(fb%p,1)
                   xx = og(i)/100
                   fc = 0; if ( lloc ) fc = xx + yy + zz + nn
                   fb%p(i,j,k,n) = sgn*(ii+fc)
                end do
             end do
          end do
          do k = fab_ilwb_n(fb,3), fab_iupb_n(fb,3)
             zz = og(k)/100/100/100
             do j = fab_ilwb_n(fb,2), fab_iupb_n(fb,2)
                yy = og(j)/100/100
                do i = fab_ilwb_n(fb,1), fab_iupb_n(fb,1)
                   xx = og(i)/100
                   fc = 0; if ( lloc ) fc = xx + yy + zz + nn
                   fb%p(i,j,k,n) = ii+fc
                end do
             end do
          end do
       case (2)
          nn = og(n)/100/100/100
          do j = lbound(fb%p,2), ubound(fb%p,2)
             yy = og(j)/100/100
             do i = lbound(fb%p,1), ubound(fb%p,1)
                xx = og(i)/100
                fc = 0; if ( lloc ) fc = xx + yy + nn
                fb%p(i,j,1,n) = sgn*(ii+fc)
             end do
          end do
          do j = fab_ilwb_n(fb,2), fab_iupb_n(fb,2)
             yy = og(j)/100/100
             do i = fab_ilwb_n(fb,1), fab_iupb_n(fb,1)
                xx = og(i)/100
                fc = 0; if ( lloc ) fc = xx + yy + nn
                fb%p(i,j,1,n) = ii+fc
             end do
          end do
       case (1)
          nn = og(n)/100/100
          do i = lbound(fb%p,1), ubound(fb%p,1)
             xx = og(i)/100
             fc = 0; if ( lloc ) fc = xx + nn
             fb%p(i,1,1,n) = sgn*(ii+fc)
          end do
          do i = fab_ilwb_n(fb,1), fab_iupb_n(fb,1)
             xx = og(i)/100
             fc = 0; if ( lloc ) fc = xx + nn
             fb%p(i,1,1,n) = ii+fc
          end do
       end select
    end do
  contains
    function og(ii) result(r)
      real(kind=dp_t) :: r
      integer, intent(in) :: ii
      if ( ii >= 0 ) then
         r = ii
      else
         r = 100-abs(ii)
      end if
    end function og
      
  end subroutine fab_debug_fill

  function fab_dataptr_c(fb, c, nc) result(r)
    type(fab), intent(in) :: fb
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    real(kind=dp_t), pointer :: r(:,:,:,:)
    integer :: lnc
    lnc = 1; if ( present(nc) ) lnc = nc
    if ( (c+lnc-1) > fb%nc ) call bl_error('fab_dataptr_c: not enough components')
    r => fb%p(:,:,:,c:c+lnc-1)
  end function fab_dataptr_c
  function ifab_dataptr_c(fb, c, nc) result(r)
    type(ifab), intent(in) :: fb
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    integer, pointer :: r(:,:,:,:)
    integer :: lnc
    lnc = 1; if ( present(nc) ) lnc = nc
    if ( (c+lnc-1) > fb%nc ) call bl_error('ifab_dataptr_c: not enough components')
    r => fb%p(:,:,:,c:c+lnc-1)
  end function ifab_dataptr_c
  function lfab_dataptr_c(fb, c, nc) result(r)
    type(lfab), intent(in) :: fb
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    logical, pointer :: r(:,:,:,:)
    integer :: lnc
    lnc = 1; if ( present(nc) ) lnc = nc
    if ( (c+lnc-1) > fb%nc ) call bl_error('lfab_dataptr_c: not enough components')
    r => fb%p(:,:,:,c:c+lnc-1)
  end function lfab_dataptr_c

  function fab_dataptr_bx_c(fb, bx, c, nc) result(r)
    type(fab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    real(kind=dp_t), pointer :: r(:,:,:,:)
    integer :: lnc
    lnc = 1; if ( present(nc) ) lnc = nc
    if ( (c+lnc-1) > fb%nc ) call bl_error('fab_dataptr_bx_c: not enough components')
    if ( .not. contains(fb%pbx, bx) ) call bl_error('fab_dataptr_bx_c: bx is too large')
    select case (fb%dim)
    case (1)
       r => fb%p(bx%lo(1):bx%hi(1),:,:,c:c+lnc-1)
    case (2)
       r => fb%p(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),:,c:c+lnc-1)
    case (3)
       r => fb%p(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),c:c+lnc-1)
    end select
  end function fab_dataptr_bx_c
  function ifab_dataptr_bx_c(fb, bx, c, nc) result(r)
    type(ifab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    integer, pointer :: r(:,:,:,:)
    integer :: lnc
    lnc = 1; if ( present(nc) ) lnc = nc
    if ( (c+lnc-1) > fb%nc ) call bl_error('ifab_dataptr_bx_c: not enough components')
    if ( .not. contains(fb%pbx, bx) ) call bl_error('ifab_dataptr_bx_c: bx is too large')
    select case (fb%dim)
    case (1)
       r => fb%p(bx%lo(1):bx%hi(1),:,:,c:c+lnc-1)
    case (2)
       r => fb%p(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),:,c:c+lnc-1)
    case (3)
       r => fb%p(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),c:c+lnc-1)
    end select
  end function ifab_dataptr_bx_c
  function lfab_dataptr_bx_c(fb, bx, c, nc) result(r)
    type(lfab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    logical, pointer :: r(:,:,:,:)
    integer :: lnc
    lnc = 1; if ( present(nc) ) lnc = nc
    if ( (c+lnc-1) > fb%nc ) call bl_error('lfab_dataptr_bx_c: not enough components')
    if ( .not. contains(fb%pbx, bx) ) call bl_error('lfab_dataptr_bx_c: bx is too large')
    select case (fb%dim)
    case (1)
       r => fb%p(bx%lo(1):bx%hi(1),:,:,c:c+lnc-1)
    case (2)
       r => fb%p(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),:,c:c+lnc-1)
    case (3)
       r => fb%p(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),c:c+lnc-1)
    end select
  end function lfab_dataptr_bx_c

  subroutine fab_setval(fb, val)
    type(fab), intent(inout) :: fb
    real(kind=dp_t), intent(in) :: val
    if ( .not. associated(fb%p) ) call bl_error("FAB_SETVAL: not associated")
    fb%p = val
  end subroutine fab_setval
  subroutine ifab_setval(fb, val)
    type(ifab), intent(inout) :: fb
    integer, intent(in) :: val
    if ( .not. associated(fb%p) ) call bl_error("FAB_SETVAL: not associated")
    fb%p = val
  end subroutine ifab_setval
  subroutine lfab_setval(fb, val)
    type(lfab), intent(inout) :: fb
    logical, intent(in) :: val
    if ( .not. associated(fb%p) ) call bl_error("FAB_SETVAL: not associated")
    fb%p = val
  end subroutine lfab_setval

  subroutine fab_setval_bx(fb, val, bx)
    type(fab), intent(inout) :: fb
    type(box), intent(in) :: bx
    real(kind=dp_t), intent(in) :: val
    real(kind=dp_t), pointer :: p(:,:,:,:)
    p => fab_dataptr_bx(fb, bx)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine fab_setval_bx
  subroutine ifab_setval_bx(fb, val, bx)
    type(ifab), intent(inout) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: val
    integer, pointer :: p(:,:,:,:)
    p => ifab_dataptr_bx(fb, bx)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine ifab_setval_bx
  subroutine lfab_setval_bx(fb, val, bx)
    type(lfab), intent(inout) :: fb
    type(box), intent(in) :: bx
    logical, intent(in) :: val
    logical, pointer :: p(:,:,:,:)
    p => lfab_dataptr_bx(fb, bx)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine lfab_setval_bx

  subroutine fab_setval_c(fb, val, c, nc)
    type(fab), intent(inout) :: fb
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    real(kind=dp_t), intent(in) :: val
    real(kind=dp_t), pointer :: p(:,:,:,:)
    p => fab_dataptr_c(fb, c, nc)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine fab_setval_c
  subroutine ifab_setval_c(fb, val, c, nc)
    type(ifab), intent(inout) :: fb
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    integer, intent(in) :: val
    integer, pointer :: p(:,:,:,:)
    p => ifab_dataptr_c(fb, c, nc)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine ifab_setval_c
  subroutine lfab_setval_c(fb, val, c, nc)
    type(lfab), intent(inout) :: fb
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    logical, intent(in) :: val
    logical, pointer :: p(:,:,:,:)
    p => lfab_dataptr_c(fb, c, nc)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine lfab_setval_c

  subroutine fab_setval_bx_c(fb, val, bx, c, nc)
    type(fab), intent(inout) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    real(kind=dp_t), intent(in) :: val
    real(kind=dp_t), pointer :: p(:,:,:,:)
    p => fab_dataptr_bx_c(fb, bx, c, nc)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine fab_setval_bx_c
  subroutine ifab_setval_bx_c(fb, val, bx, c, nc)
    type(ifab), intent(inout) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    integer, intent(in) :: val
    integer, pointer :: p(:,:,:,:)
    p => ifab_dataptr_bx_c(fb, bx, c, nc)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine ifab_setval_bx_c
  subroutine lfab_setval_bx_c(fb, val, bx, c, nc)
    type(lfab), intent(inout) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    logical, intent(in) :: val
    logical, pointer :: p(:,:,:,:)
    p => lfab_dataptr_bx_c(fb, bx, c, nc)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine lfab_setval_bx_c

  subroutine fab_print(fb, str, unit, all, data, bx, skip)
    use bl_IO_module
    type(fab), intent(in) :: fb
    character(len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: all, data
    integer, intent(in), optional :: skip
    type(box), intent(in), optional :: bx
    integer :: un
    logical :: lall, ldata
    type(box) :: lbx
    lbx  = box_allbox(fb%dim); if ( present(bx) ) lbx  = bx
    lall = .TRUE.; if ( present(all) ) lall = all
    ldata = .TRUE.; if ( present(data) ) ldata = data
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un, fmt='("FAB")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(": ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" DIM     = ",i2)') fb%dim
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NG      = ",i2)') fb%ng
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NC      = ",i2)') fb%nc
    call unit_skip(un, skip)
    write(unit=un, fmt='(" IBX     = ",i2)', advance = 'no')
    call print(fb%ibx, unit = un)
    call unit_skip(un, skip)
    write(unit=un, fmt='(" PBX     = ",i2)', advance = 'no')
    call print(fb%pbx, unit = un)
    call unit_skip(un, skip)
    write(unit=un, fmt='(" BX      = ",i2)', advance = 'no')
    call print(fb%bx, unit = un)
    if ( .not. associated(fb%p) ) then
       call unit_skip(un, skip)
       write(unit=un) 'NOT ASSOCIATED'
    else
       select case (fb%dim)
       case (1)
          call print_1d(fb%p(:,1,1,:), lbound(fb%p), intersection(fb%ibx,lbx))
       case (2)
          call print_2d(fb%p(:,:,1,:), lbound(fb%p), intersection(fb%ibx,lbx))
       case (3)
          call print_3d(fb%p(:,:,:,:), lbound(fb%p), intersection(fb%ibx,lbx))
       end select
    end if
  contains
    subroutine print_1d(fb, lo, bx)
      integer, intent(in) :: lo(:)
      type(box), intent(in) :: bx
      real(kind=dp_t), intent(in) :: fb(lo(1):,:)
      integer ::n, i
      integer :: nc, hi(1)
      character(len=1) c
      nc = size(fb,dim=2)
      hi(1) = lo(1) + size(fb,dim=1) - 1
      if ( ldata ) then
         do n = 1, nc
            do i = lo(1), hi(1)
               if ( .not. ( lall .or. contains(bx, (/i/)) ) ) cycle
               c = ' '
               if ( .not. contains(bx, (/i/)) ) c = '*'
               call unit_skip(un, skip)
               write(unit=un, fmt='(A1,1X,I3,1(1X,I5),1X,G25.15)') &
                    c, n, i, fb(i,n)
            end do
         end do
      end if
    end subroutine print_1d
    subroutine print_2d(fb, lo, bx)
      integer, intent(in) :: lo(:)
      type(box), intent(in) :: bx
      real(kind=dp_t), intent(in) :: fb(lo(1):,lo(2):,:)
      integer :: n, j, i
      integer :: nc, hi(2)
      character(len=1) c
      nc = size(fb,dim=3)
      do i = 1, 2
         hi(i) = lo(i) + size(fb,dim=i) - 1
      end do
      if ( ldata ) then
         do n = 1, nc
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if ( .not. ( lall .or. contains(bx, (/i,j/)) ) ) cycle
                  c = ' '
                  if ( .not. contains(bx, (/i,j/)) ) c = '*'
                  call unit_skip(un, skip)
                  write(unit=un, fmt='(A1,1X,I3,2(1X,I5),1X,G25.15)') &
                       c, n, i, j, fb(i,j,n)
               end do
            end do
         end do
      end if
    end subroutine print_2d
    subroutine print_3d(fb, lo, bx)
      integer, intent(in) :: lo(:)
      type(box), intent(in) :: bx
      real(kind=dp_t), intent(in) :: fb(lo(1):,lo(2):,lo(3):,:)
      integer :: n, k, j, i
      integer :: nc, hi(3)
      character(len=1) :: c
      nc = size(fb,dim=4)
      do i = 1, 3
         hi(i) = lo(i) + size(fb,dim=i) - 1
      end do
      if ( ldata ) then
         do n = 1, nc
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if ( .not. ( lall .or. contains(bx, (/i,j,k/)) ) ) cycle
                     c = ' '
                     if ( .not. contains(bx, (/i,j,k/)) ) c = '*'
                     call unit_skip(un, skip)
                     write(unit=un, fmt='(A1,1X,I3,3(1X,I5),1X,G25.15)') &
                          c, n, i, j, k, fb(i,j,k,n)
                  end do
               end do
            end do
         end do
      end if
    end subroutine print_3d
  end subroutine fab_print
  subroutine ifab_print(fb, str, unit, all, data, bx, skip)
    use bl_IO_module
    type(ifab), intent(in) :: fb
    character(len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: all, data
    integer, intent(in), optional :: skip
    type(box), intent(in), optional :: bx
    integer :: un
    logical :: lall, ldata
    type(box) :: lbx
    lbx  = box_allbox(fb%dim); if ( present(bx) ) lbx  = bx
    lall = .TRUE.; if ( present(all) ) lall = all
    ldata = .TRUE.; if ( present(data) ) ldata = data
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un, fmt='("IFAB")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(": ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" DIM     = ",i2)') fb%dim
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NG      = ",i2)') fb%ng
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NC      = ",i2)') fb%nc
    call unit_skip(un, skip)
    write(unit=un, fmt='(" IBX     = ",i2)', advance = 'no')
    call print(fb%ibx, unit = un)
    call unit_skip(un, skip)
    write(unit=un, fmt='(" PBX     = ",i2)', advance = 'no')
    call print(fb%pbx, unit = un)
    call unit_skip(un, skip)
    write(unit=un, fmt='(" BX      = ",i2)', advance = 'no')
    call print(fb%bx, unit = un)
    if ( .not. associated(fb%p) ) then
       call unit_skip(un, skip)
       write(unit=un) 'NOT ASSOCIATED'
    else
       select case (fb%dim)
       case (1)
          call print_1d(fb%p(:,1,1,:), lbound(fb%p), intersection(fb%ibx,lbx))
       case (2)
          call print_2d(fb%p(:,:,1,:), lbound(fb%p), intersection(fb%ibx,lbx))
       case (3)
          call print_3d(fb%p(:,:,:,:), lbound(fb%p), intersection(fb%ibx,lbx))
       end select
    end if
  contains
    subroutine print_1d(fb, lo, bx)
      integer, intent(in) :: lo(:)
      type(box), intent(in) :: bx
      integer, intent(in) :: fb(lo(1):,:)
      integer n, i
      integer nc, hi(1)
      character(len=1) c
      nc = size(fb,dim=2)
      hi(1) = lo(1) + size(fb,dim=1) - 1
      if ( ldata ) then
         do n = 1, nc
            do i = lo(1), hi(1)
               if ( .not. ( lall .or. contains(bx, (/i/)) ) ) cycle
               c = ' '
               if ( .not. contains(bx, (/i/)) ) c = '*'
               call unit_skip(un, skip)
               write(unit=un, fmt='(A1,1X,I3,1(1X,I5),1X,G25.15)') &
                    c, n, i, fb(i,n)
            end do
         end do
      end if
    end subroutine print_1d
    subroutine print_2d(fb, lo, bx)
      integer, intent(in) :: lo(:)
      type(box), intent(in) :: bx
      integer, intent(in) :: fb(lo(1):,lo(2):,:)
      integer n, j, i
      integer nc, hi(2)
      character(len=1) c
      nc = size(fb,dim=3)
      do i = 1, 2
         hi(i) = lo(i) + size(fb,dim=i) - 1
      end do
      if ( ldata ) then
         do n = 1, nc
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if ( .not. ( lall .or. contains(bx, (/i,j/)) ) ) cycle
                  c = ' '
                  if ( .not. contains(bx, (/i,j/)) ) c = '*'
                  call unit_skip(un, skip)
                  write(unit=un, fmt='(A1,1X,I3,2(1X,I5),1X,G25.15)') &
                       c, n, i, j, fb(i,j,n)
               end do
            end do
         end do
      end if
    end subroutine print_2d
    subroutine print_3d(fb, lo, bx)
      integer, intent(in) :: lo(:)
      type(box), intent(in) :: bx
      integer, intent(in) :: fb(lo(1):,lo(2):,lo(3):,:)
      integer :: n, k, j, i
      integer :: nc, hi(3)
      character(len=1) :: c
      nc = size(fb,dim=4)
      do i = 1, 3
         hi(i) = lo(i) + size(fb,dim=i) - 1
      end do
      if ( ldata ) then
         do n = 1, nc
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if ( .not. ( lall .or. contains(bx, (/i,j,k/)) ) ) cycle
                     c = ' '
                     if ( .not. contains(bx, (/i,j,k/)) ) c = '*'
                     call unit_skip(un, skip)
                     write(unit=un, fmt='(A1,1X,I3,3(1X,I5),1X,G25.15)') &
                          c, n, i, j, k, fb(i,j,k,n)
                  end do
               end do
            end do
         end do
      end if
    end subroutine print_3d
  end subroutine ifab_print
  subroutine lfab_print(fb, str, unit, all, data, bx, skip)
    use bl_IO_module
    type(lfab), intent(in) :: fb
    character(len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: all, data
    integer, intent(in), optional :: skip
    type(box), intent(in), optional :: bx
    integer :: un
    logical :: lall, ldata
    type(box) :: lbx
    lbx  = box_allbox(fb%dim); if ( present(bx) ) lbx  = bx
    lall = .TRUE.; if ( present(all) ) lall = all
    ldata = .TRUE.; if ( present(data) ) ldata = data
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un, fmt='("LFAB")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(": ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" DIM     = ",i2)') fb%dim
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NG      = ",i2)') fb%ng
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NC      = ",i2)') fb%nc
    call unit_skip(un, skip)
    write(unit=un, fmt='(" IBX     = ",i2)', advance = 'no')
    call print(fb%ibx, unit = un)
    call unit_skip(un, skip)
    write(unit=un, fmt='(" PBX     = ",i2)', advance = 'no')
    call print(fb%pbx, unit = un)
    call unit_skip(un, skip)
    write(unit=un, fmt='(" BX      = ",i2)', advance = 'no')
    call print(fb%bx, unit = un)
    if ( .not. associated(fb%p) ) then
       call unit_skip(un, skip)
       write(unit=un) 'NOT ASSOCIATED'
    else
       select case (fb%dim)
       case (1)
          call print_1d(fb%p(:,1,1,:), lbound(fb%p), intersection(fb%ibx, lbx))
       case (2)
          call print_2d(fb%p(:,:,1,:), lbound(fb%p), intersection(fb%ibx, lbx))
       case (3)
          call print_3d(fb%p(:,:,:,:), lbound(fb%p), intersection(fb%ibx, lbx))
       end select
    end if
  contains
    subroutine print_1d(fb, lo, bx)
      integer, intent(in) :: lo(:)
      type(box), intent(in) :: bx
      logical, intent(in) :: fb(lo(1):,:)
      integer :: n, i
      integer :: nc, hi(1)
      character(len=1) :: c
      nc = size(fb,dim=2)
      hi(1) = lo(1) + size(fb,dim=1) - 1
      if ( ldata ) then
         do n = 1, nc
            do i = lo(1), hi(1)
               if ( .not. ( lall .or. contains(bx, (/i/)) ) ) cycle
               c = ' '
               if ( .not. contains(bx, (/i/)) ) c = '*'
               call unit_skip(un, skip)
               write(unit=un, fmt='(A1,1X,I3,1(1X,I5),1X,G25.15)') &
                    c, n, i, fb(i,n)
            end do
         end do
      end if
    end subroutine print_1d
    subroutine print_2d(fb, lo, bx)
      integer, intent(in) :: lo(:)
      type(box), intent(in) :: bx
      logical, intent(in) :: fb(lo(1):,lo(2):,:)
      integer :: n, j, i
      integer :: nc, hi(2)
      character(len=1) :: c
      nc = size(fb,dim=3)
      do i = 1, 2
         hi(i) = lo(i) + size(fb,dim=i) - 1
      end do
      if ( ldata ) then
         do n = 1, nc
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if ( .not. ( lall .or. contains(bx, (/i,j/)) ) ) cycle
                  c = ' '
                  if ( .not. contains(bx, (/i,j/)) ) c = '*'
                  call unit_skip(un, skip)
                  write(unit=un, fmt='(A1,1X,I3,2(1X,I5),1X,G25.15)') &
                       c, n, i, j, fb(i,j,n)
               end do
            end do
         end do
      end if
    end subroutine print_2d
    subroutine print_3d(fb, lo, bx)
      integer, intent(in) :: lo(:)
      type(box), intent(in) :: bx
      logical, intent(in) :: fb(lo(1):,lo(2):,lo(3):,:)
      integer :: n, k, j, i
      integer :: nc, hi(3)
      character(len=1) :: c
      nc = size(fb,dim=4)
      do i = 1, 3
         hi(i) = lo(i) + size(fb,dim=i) - 1
      end do
      if ( ldata ) then
         do n = 1, nc
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if ( .not. ( lall .or. contains(bx, (/i,j,k/)) ) ) cycle
                     c = ' '
                     if ( .not. contains(bx, (/i,j,k/)) ) c = '*'
                     call unit_skip(un, skip)
                     write(unit=un, fmt='(A1,1X,I3,3(1X,I5),1X,G25.15)') &
                          c, n, i, j, k, fb(i,j,k,n)
                  end do
               end do
            end do
         end do
      end if
    end subroutine print_3d
  end subroutine lfab_print

  function fab_max(fb, all) result(r)
    real(kind=dp_t) :: r
    type(fab), intent(in) :: fb
    logical, intent(in), optional :: all
    logical :: lall
    real(dp_t), pointer :: mp(:,:,:,:)
    lall = .FALSE.; if ( present(all) ) lall = all
    if ( lall ) then
       mp => dataptr(fb, get_pbox(fb))
    else
       mp => dataptr(fb, get_ibox(fb))
    end if
    r = maxval(mp)
  end function fab_max
  function fab_max_c(fb, c, all) result(r)
    real(kind=dp_t) :: r
    type(fab), intent(in) :: fb
    integer, intent(in) :: c
    logical, intent(in), optional :: all
    logical :: lall
    real(dp_t), pointer :: mp(:,:,:,:)
    lall = .FALSE.; if ( present(all) ) lall = all
    if ( lall ) then
       mp => dataptr(fb, get_pbox(fb), c)
    else
       mp => dataptr(fb, get_ibox(fb), c)
    end if
    r = maxval(mp)
  end function fab_max_c
  function fab_max_bx(fb, bx) result(r)
    real(kind=dp_t) :: r
    type(fab), intent(in) :: fb
    type(box), intent(in) :: bx
    real(dp_t), pointer :: mp(:,:,:,:)
    type(box) :: sbx
    sbx = intersection(bx, get_ibox(fb))
    mp => dataptr(fb, sbx)
    r = maxval(mp)
  end function fab_max_bx
  function fab_max_bx_c(fb, bx, c) result(r)
    real(kind=dp_t) :: r
    type(fab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    real(dp_t), pointer :: mp(:,:,:,:)
    type(box) :: sbx
    sbx = intersection(bx, get_ibox(fb))
    mp => dataptr(fb, sbx, c)
    r = maxval(mp)
  end function fab_max_bx_c
  function ifab_max(fb, all) result(r)
    integer :: r
    type(ifab), intent(in) :: fb
    logical, intent(in), optional :: all
    logical :: lall
    integer, pointer :: mp(:,:,:,:)
    lall = .FALSE.; if ( present(all) ) lall = all
    if ( lall ) then
       mp => dataptr(fb, get_pbox(fb))
    else
       mp => dataptr(fb, get_ibox(fb))
    end if
    r = maxval(mp)
  end function ifab_max
  function ifab_max_c(fb, c, all) result(r)
    integer :: r
    type(ifab), intent(in) :: fb
    integer, intent(in) :: c
    logical, intent(in), optional :: all
    logical :: lall
    integer, pointer :: mp(:,:,:,:)
    lall = .FALSE.; if ( present(all) ) lall = all
    if ( lall ) then
       mp => dataptr(fb, get_pbox(fb), c)
    else
       mp => dataptr(fb, get_ibox(fb), c)
    end if
    r = maxval(mp)
  end function ifab_max_c
  function ifab_max_bx(fb, bx) result(r)
    integer :: r
    type(ifab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, pointer :: mp(:,:,:,:)
    type(box) :: sbx
    sbx = intersection(bx, get_ibox(fb))
    mp => dataptr(fb, sbx)
    r = maxval(mp)
  end function ifab_max_bx
  function ifab_max_bx_c(fb, bx, c) result(r)
    integer :: r
    type(ifab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, pointer :: mp(:,:,:,:)
    type(box) :: sbx
    sbx = intersection(bx, get_ibox(fb))
    mp => dataptr(fb, sbx, c)
    r = maxval(mp)
  end function ifab_max_bx_c

  function fab_min(fb, all) result(r)
    real(kind=dp_t) :: r
    type(fab), intent(in) :: fb
    logical, intent(in), optional :: all
    logical :: lall
    real(dp_t), pointer :: mp(:,:,:,:)
    lall = .FALSE.; if ( present(all) ) lall = all
    if ( lall ) then
       mp => dataptr(fb, get_pbox(fb))
    else
       mp => dataptr(fb, get_ibox(fb))
    end if
    r = minval(mp)
  end function fab_min
  function fab_min_c(fb, c, all) result(r)
    real(kind=dp_t) :: r
    type(fab), intent(in) :: fb
    integer, intent(in) :: c
    logical, intent(in), optional :: all
    logical :: lall
    real(dp_t), pointer :: mp(:,:,:,:)
    lall = .FALSE.; if ( present(all) ) lall = all
    if ( lall ) then
       mp => dataptr(fb, get_pbox(fb), c)
    else
       mp => dataptr(fb, get_ibox(fb), c)
    end if
    r = minval(mp)
  end function fab_min_c
  function fab_min_bx(fb, bx) result(r)
    real(kind=dp_t) :: r
    type(fab), intent(in) :: fb
    type(box), intent(in) :: bx
    real(dp_t), pointer :: mp(:,:,:,:)
    type(box) :: sbx
    sbx = intersection(bx, get_ibox(fb))
    mp => dataptr(fb, sbx)
    r = minval(mp)
  end function fab_min_bx
  function fab_min_bx_c(fb, bx, c) result(r)
    real(kind=dp_t) :: r
    type(fab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    real(dp_t), pointer :: mp(:,:,:,:)
    type(box) :: sbx
    sbx = intersection(bx, get_ibox(fb))
    mp => dataptr(fb, sbx, c)
    r = minval(mp)
  end function fab_min_bx_c
  function ifab_min(fb, all) result(r)
    integer :: r
    type(ifab), intent(in) :: fb
    logical, intent(in), optional :: all
    logical :: lall
    integer, pointer :: mp(:,:,:,:)
    lall = .FALSE.; if ( present(all) ) lall = all
    if ( lall ) then
       mp => dataptr(fb, get_pbox(fb))
    else
       mp => dataptr(fb, get_ibox(fb))
    end if
    r = minval(mp)
  end function ifab_min
  function ifab_min_c(fb, c, all) result(r)
    integer :: r
    type(ifab), intent(in) :: fb
    integer, intent(in) :: c
    logical, intent(in), optional :: all
    logical :: lall
    integer, pointer :: mp(:,:,:,:)
    lall = .FALSE.; if ( present(all) ) lall = all
    if ( lall ) then
       mp => dataptr(fb, get_pbox(fb), c)
    else
       mp => dataptr(fb, get_ibox(fb), c)
    end if
    r = minval(mp)
  end function ifab_min_c
  function ifab_min_bx(fb, bx) result(r)
    integer :: r
    type(ifab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, pointer :: mp(:,:,:,:)
    type(box) :: sbx
    sbx = intersection(bx, get_ibox(fb))
    mp => dataptr(fb, sbx)
    r = minval(mp)
  end function ifab_min_bx
  function ifab_min_bx_c(fb, bx, c) result(r)
    integer :: r
    type(ifab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, pointer :: mp(:,:,:,:)
    type(box) :: sbx
    sbx = intersection(bx, get_ibox(fb))
    mp => dataptr(fb, sbx, c)
    r = minval(mp)
  end function ifab_min_bx_c

end module fab_module
