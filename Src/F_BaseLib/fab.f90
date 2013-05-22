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
     private
     integer   :: dim = 0
     type(box) :: bx
     type(box) :: pbx
     type(box) :: ibx
     integer   :: nc = 1
     real(kind = dp_t), pointer, dimension(:,:,:,:) :: p => Null()
  end type fab

  type zfab
     private
     integer   :: dim = 0
     type(box) :: bx
     type(box) :: pbx
     type(box) :: ibx
     integer   :: nc = 1
     complex(kind = dp_t), pointer, dimension(:,:,:,:) :: p => Null()
  end type zfab

  type ifab
     private
     integer   :: dim = 0
     type(box) :: bx
     type(box) :: pbx
     type(box) :: ibx
     integer   :: nc = 1
     integer, pointer, dimension(:,:,:,:) :: p => Null()
  end type ifab

  type lfab
     private
     integer   :: dim = 0
     type(box) :: bx
     type(box) :: pbx
     type(box) :: ibx
     integer   :: nc = 1
     logical, pointer, dimension(:,:,:,:) :: p => Null()
  end type lfab

  !! Returns the Dimension  of the FAB
  interface get_dim
     module procedure fab_dim
     module procedure ifab_dim
     module procedure lfab_dim
     module procedure zfab_dim
  end interface

  interface build
     module procedure fab_build
     module procedure ifab_build
     module procedure lfab_build
     module procedure zfab_build
  end interface

  interface destroy
     module procedure fab_destroy
     module procedure ifab_destroy
     module procedure lfab_destroy
     module procedure zfab_destroy
  end interface

  !! Returns whether the FAB has been built, this is different from
  !! whether the underlying pointer has been allocated since on parallel
  !! systems the pointer will not be allocated if the LAYOUT says it
  !! lives on another processor
  interface built_q
     module procedure fab_built_q
     module procedure zfab_built_q
     module procedure ifab_built_q
     module procedure lfab_built_q
  end interface

  interface max_val
     module procedure fab_max_val
     module procedure fab_max_val_bx
     module procedure fab_max_val_c
     module procedure fab_max_val_bx_c

     module procedure ifab_max_val
     module procedure ifab_max_val_bx
     module procedure ifab_max_val_c
     module procedure ifab_max_val_bx_c
  end interface

  interface min_val
     module procedure fab_min_val
     module procedure fab_min_val_bx
     module procedure fab_min_val_c
     module procedure fab_min_val_bx_c

     module procedure ifab_min_val
     module procedure ifab_min_val_bx
     module procedure ifab_min_val_c
     module procedure ifab_min_val_bx_c
  end interface

  interface dataptr
     module procedure fab_dataptr
     module procedure fab_dataptr_bx
     module procedure fab_dataptr_c
     module procedure fab_dataptr_bx_c

     module procedure zfab_dataptr
     module procedure zfab_dataptr_bx
     module procedure zfab_dataptr_c
     module procedure zfab_dataptr_bx_c

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

     module procedure zfab_setval
     module procedure zfab_setval_bx
     module procedure zfab_setval_c
     module procedure zfab_setval_bx_c

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
     module procedure zfab_print
  end interface

  interface get_box
     module procedure fab_get_box
     module procedure zfab_get_box
     module procedure ifab_get_box
     module procedure lfab_get_box
  end interface

  interface get_pbox
     module procedure fab_get_pbox
     module procedure zfab_get_pbox
     module procedure ifab_get_pbox
     module procedure lfab_get_pbox
  end interface

  interface get_ibox
     module procedure fab_get_ibox
     module procedure zfab_get_pbox
     module procedure ifab_get_ibox
     module procedure lfab_get_ibox
  end interface

  interface volume
     module procedure fab_volume
     module procedure zfab_volume
     module procedure ifab_volume
     module procedure lfab_volume
  end interface

  interface ncomp
     module procedure fab_ncomp
     module procedure zfab_ncomp
     module procedure ifab_ncomp
     module procedure lfab_ncomp
  end interface

  interface contains_nan
     module procedure contains_nan_allc
     module procedure contains_nan_c
     module procedure contains_nan_bx_c
  end interface contains_nan

  interface contains_inf
     module procedure contains_inf_allc
     module procedure contains_inf_c
     module procedure contains_inf_bx_c
  end interface contains_inf

  logical,       private, save ::      do_init_fabs = .false.
  real(dp_t),    private, save ::  fab_default_init = -Huge(1.0_dp_t)
  complex(dp_t), private, save :: zfab_default_init = -Huge(1.0_dp_t)
  integer,       private, save :: ifab_default_init = -Huge(1)
  logical,       private, save :: lfab_default_init = .false.

  type(mem_stats), private, save ::  fab_ms
  type(mem_stats), private, save :: zfab_ms
  type(mem_stats), private, save :: ifab_ms
  type(mem_stats), private, save :: lfab_ms
  !
  ! The high-water mark for number of double precision values allocated in fabs.
  !
  integer(ll_t), private, save :: fab_high_water_mark = 0

  private ::  fab_minval_doit,  fab_maxval_doit

contains
  !
  ! To tell whether or not a value is a IEEE inf.
  !
  ! In time this should be rewritten to use the IEEE_ARITHMETIC module.
  !
  function is_an_inf (val) result(r)

    real(dp_t), intent(in) :: val
    logical                :: r

    integer :: rc

    interface
       subroutine val_is_inf(v, res)
         use bl_types
         real(dp_t), intent(in)  :: v
         integer,    intent(out) :: res
       end subroutine val_is_inf
    end interface

    r = .false.

    call val_is_inf(val,rc)

    if (rc == 1) r = .true.
    
  end function is_an_inf
  !
  ! To tell whether or not a value is a IEEE NaN.
  !
  ! In time this should be rewritten to use the IEEE_ARITHMETIC module.
  !
  function is_a_nan (val) result(r)

    real(dp_t), intent(in) :: val
    logical                :: r

    integer :: rc

    interface
       subroutine val_is_nan(v, res)
         use bl_types
         real(dp_t), intent(in)  :: v
         integer,    intent(out) :: res
       end subroutine val_is_nan
    end interface

    r = .false.

    call val_is_nan(val,rc)

    if (rc == 1) r = .true.
    
  end function is_a_nan
  !
  ! Does a real(dp_t) fab contain a NaN?
  !
  function contains_nan_c(fb,c,nc) result(r)

    use bl_error_module

    logical               :: r
    type(fab), intent(in) :: fb
    integer,   intent(in) :: c, nc

    integer                  :: sz, rc
    real(kind=dp_t), pointer :: pp(:,:,:,:)

    interface
       subroutine fab_contains_nan(dptr, count, res)
         use bl_types
         integer,    intent(in)  :: count
         real(dp_t), intent(in)  :: dptr(count)
         integer,    intent(out) :: res
       end subroutine fab_contains_nan
    end interface

    if ( (c+nc-1) > fb%nc ) call bl_error('contains_nan_c: not enough components')

    r = .false.

    pp => dataptr(fb)

    sz = volume(get_pbox(fb)) * nc

    call fab_contains_nan(pp(:,:,:,c), sz, rc)

    if (rc == 1) r = .true.

  end function contains_nan_c

  function contains_nan_allc(fb) result(r)
    logical               :: r
    type(fab), intent(in) :: fb
    r = contains_nan_c(fb,1,ncomp(fb))
  end function contains_nan_allc

  function contains_nan_bx_c(fb,bx,c,nc) result(r)

    logical               :: r
    type(fab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer,   intent(in) :: c, nc

    integer                      :: sz, rc, i, j, k, n, idx, lo(4), hi(4)
    real(kind=dp_t), allocatable :: d(:)
    real(kind=dp_t), pointer     :: pp(:,:,:,:)

    interface
       subroutine fab_contains_nan(dptr, count, res)
         use bl_types
         integer,    intent(in)  :: count
         real(dp_t), intent(in)  :: dptr(count)
         integer,    intent(out) :: res
       end subroutine fab_contains_nan
    end interface

    r = .false.

    sz = volume(bx) * nc

    allocate(d(sz))

    pp => dataptr(fb,bx,c,nc)

    lo = lbound(pp)
    hi = ubound(pp)

    idx = 1
    do n = lo(4), hi(4)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                d(idx) = pp(i,j,k,n)
                idx = idx + 1
             end do
          end do
       end do
    end do

    call fab_contains_nan(d,sz,rc)

    if (rc == 1) r = .true.

  end function contains_nan_bx_c
  !
  ! Does a real(dp_t) fab contain a Inf?
  !
  function contains_inf_c(fb,c,nc) result(r)

    use bl_error_module

    logical               :: r
    type(fab), intent(in) :: fb
    integer,   intent(in) :: c, nc

    integer                  :: sz, rc
    real(kind=dp_t), pointer :: pp(:,:,:,:)

    interface
       subroutine fab_contains_inf(dptr, count, res)
         use bl_types
         integer,    intent(in)  :: count
         real(dp_t), intent(in)  :: dptr(count)
         integer,    intent(out) :: res
       end subroutine fab_contains_inf
    end interface

    if ( (c+nc-1) > fb%nc ) call bl_error('contains_inf_c: not enough components')

    r = .false.

    pp => dataptr(fb)

    sz = volume(get_pbox(fb)) * nc

    call fab_contains_inf(pp(:,:,:,c), sz, rc)

    if (rc == 1) r = .true.

  end function contains_inf_c

  function contains_inf_allc(fb) result(r)
    logical               :: r
    type(fab), intent(in) :: fb
    r = contains_inf_c(fb,1,ncomp(fb))
  end function contains_inf_allc

  function contains_inf_bx_c(fb,bx,c,nc) result(r)

    logical               :: r
    type(fab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer,   intent(in) :: c, nc

    integer                      :: sz, rc, i, j, k, n, idx, lo(4), hi(4)
    real(kind=dp_t), allocatable :: d(:)
    real(kind=dp_t), pointer     :: pp(:,:,:,:)

    interface
       subroutine fab_contains_inf(dptr, count, res)
         use bl_types
         integer,    intent(in)  :: count
         real(dp_t), intent(in)  :: dptr(count)
         integer,    intent(out) :: res
       end subroutine fab_contains_inf
    end interface

    r = .false.

    sz = volume(bx) * nc

    allocate(d(sz))

    pp => dataptr(fb,bx,c,nc)

    lo = lbound(pp)
    hi = ubound(pp)

    idx = 1
    do n = lo(4), hi(4)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                d(idx) = pp(i,j,k,n)
                idx = idx + 1
             end do
          end do
       end do
    end do

    call fab_contains_inf(d,sz,rc)

    if (rc == 1) r = .true.

  end function contains_inf_bx_c
  !
  ! Toggle whether or not setval() is called immediately after a fab is built.
  !
  subroutine setval_fabs_on_init(yes_or_no)
    logical, intent(in) :: yes_or_no
    do_init_fabs = yes_or_no
  end subroutine setval_fabs_on_init
  !
  ! Set the default value to which fabs are setval()d when built.
  !
  subroutine set_fab_default_init(val)
    real(dp_t), intent(in) :: val
    fab_default_init = val
  end subroutine set_fab_default_init
  subroutine set_zfab_default_init(val)
    complex(dp_t), intent(in) :: val
    zfab_default_init = val
  end subroutine set_zfab_default_init
  subroutine set_ifab_default_init(val)
    integer, intent(in) :: val
    ifab_default_init = val
  end subroutine set_ifab_default_init
  subroutine set_lfab_default_init(val)
    logical, intent(in) :: val
    lfab_default_init = val
  end subroutine set_lfab_default_init

  subroutine print_and_reset_fab_byte_spread()
    use parallel
    integer    :: ioproc, ihi, ilo, iav
    real(dp_t) :: lo, hi, av
    ioproc = parallel_IOProcessorNode()
    call parallel_reduce(lo, real(fab_high_water_mark,dp_t), MPI_MIN, proc = ioproc)
    call parallel_reduce(hi, real(fab_high_water_mark,dp_t), MPI_MAX, proc = ioproc)
    call parallel_reduce(av, real(fab_high_water_mark,dp_t), MPI_SUM, proc = ioproc)
    if ( parallel_IOProcessor() ) then
       !
       ! This assumes sizeof(dp_t) == 8
       !
       ilo = int(8*lo,ll_t)
       ihi = int(8*hi,ll_t)
       iav = int(8*av,ll_t) / parallel_nprocs()
       print*, ''
       write(6,fmt='("FAB byte spread across MPI nodes: [",i11," - ",i11," [avg: ",i11,"]]")') ilo, ihi, iav
       print*, ''
    end if
    fab_high_water_mark = 0
  end subroutine print_and_reset_fab_byte_spread
  !
  ! Returns a pointer to an integer array 0:nprocs-1 of the CPU numbers
  ! [0..nprocs-1], where the CPU numbers are sorted from least amount of
  ! fab memory in use on the CPU to largest amount of fab memory in use
  ! on that CPU.  The pointer must be deallocated by the calling routine.
  !
  function least_used_cpus () result(r)

    use parallel
    use sort_i_module

    integer, pointer     :: r(:)
    integer              :: i
    integer, allocatable :: snd(:), rcv(:), idx(:)

    integer(ll_t) :: val  ! Number of double precision values stored in fabs on this CPU.

    allocate(r(0:parallel_nprocs()-1))

    allocate(snd(1), rcv(parallel_nprocs()), idx(parallel_nprocs()))

    val = fab_ms%num_alloc - fab_ms%num_dealloc

    snd(1) = int(val)

    call parallel_allgather(snd, rcv, 1)

    call sort(rcv, idx)

    r = idx - 1

    if ( .false. .and. parallel_ioprocessor() ) then
       print*, '*** least_used_cpus(): '
       do i = 0, parallel_nprocs()-1
          print*, i, ' : ', r(i)
       end do
    end if

  end function least_used_cpus

  subroutine fab_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    fab_ms = ms
  end subroutine fab_set_mem_stats
  subroutine zfab_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    zfab_ms = ms
  end subroutine zfab_set_mem_stats
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
  function zfab_mem_stats() result(r)
    type(mem_stats) :: r
    r = zfab_ms
  end function zfab_mem_stats
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
    logical, intent(in), optional :: all
    if ( all ) then
       r = volume(get_pbox(fb))
    else
       r = volume(get_box(fb))
    end if
    r = r * fb%nc
  end function fab_volume
  function zfab_volume(fb, all) result(r)
    integer(kind=ll_t) :: r
    type(zfab), intent(in) :: fb
    logical, intent(in), optional :: all
    if ( all ) then
       r = volume(get_pbox(fb))
    else
       r = volume(get_box(fb))
    end if
    r = r * fb%nc
  end function zfab_volume
  function ifab_volume(fb, all) result(r)
    integer(kind=ll_t) :: r
    type(ifab), intent(in) :: fb
    logical, intent(in), optional :: all
    if ( all ) then
       r = volume(get_pbox(fb))
    else
       r = volume(get_box(fb))
    end if
    r = r * fb%nc
  end function ifab_volume
  function lfab_volume(fb, all) result(r)
    integer(kind=ll_t) :: r
    type(lfab), intent(in) :: fb
    logical, intent(in), optional :: all
    if ( all ) then
       r = volume(get_pbox(fb))
    else
       r = volume(get_box(fb))
    end if
    r = r * fb%nc
  end function lfab_volume

  pure function fab_lwb(fb) result(r)
    type(fab), intent(in) :: fb
    integer :: r(fb%dim)
    r = lwb(fb%bx)
  end function fab_lwb
  pure function fab_lwb_n(fb,dim) result(r)
    type(fab), intent(in) :: fb
    integer, intent(in) :: dim
    integer :: r
    r = lwb(fb%bx,dim)
  end function fab_lwb_n
  pure function ifab_lwb(fb) result(r)
    type(ifab), intent(in) :: fb
    integer :: r(fb%dim)
    r = lwb(fb%bx)
  end function ifab_lwb

  pure function zfab_lwb(fb) result(r)
    type(zfab), intent(in) :: fb
    integer :: r(fb%dim)
    r = lwb(fb%bx)
  end function zfab_lwb
  pure function zfab_lwb_n(fb,dim) result(r)
    type(zfab), intent(in) :: fb
    integer, intent(in) :: dim
    integer :: r
    r = lwb(fb%bx,dim)
  end function zfab_lwb_n

  pure function fab_ilwb(fb) result(r)
    type(fab), intent(in) :: fb
    integer :: r(fb%dim)
    r = lwb(fb%ibx)
  end function fab_ilwb
  pure function fab_ilwb_n(fb,dim) result(r)
    type(fab), intent(in) :: fb
    integer, intent(in) :: dim
    integer :: r
    r = lwb(fb%ibx,dim)
  end function fab_ilwb_n
  pure function ifab_ilwb(fb) result(r)
    type(ifab), intent(in) :: fb
    integer :: r(fb%dim)
    r = lwb(fb%ibx)
  end function ifab_ilwb

  pure function fab_upb(fb) result(r)
    type(fab), intent(in) :: fb
    integer :: r(fb%dim)
    r = upb(fb%bx)
  end function fab_upb
  pure function fab_upb_n(fb,dim) result(r)
    type(fab), intent(in) :: fb
    integer, intent(in) :: dim
    integer :: r
    r = upb(fb%bx, dim)
  end function fab_upb_n
  pure function ifab_upb(fb) result(r)
    type(ifab), intent(in) :: fb
    integer :: r(fb%dim)
    r = upb(fb%bx)
  end function ifab_upb

  pure function fab_iupb(fb) result(r)
    type(fab), intent(in) :: fb
    integer :: r(fb%dim)
    r = upb(fb%ibx)
  end function fab_iupb
  pure function fab_iupb_n(fb,dim) result(r)
    type(fab), intent(in) :: fb
    integer, intent(in) :: dim
    integer :: r
    r = upb(fb%ibx, dim)
  end function fab_iupb_n
  pure function ifab_iupb(fb) result(r)
    type(ifab), intent(in) :: fb
    integer :: r(fb%dim)
    r = upb(fb%ibx)
  end function ifab_iupb

  pure function fab_plwb(fb) result(r)
    type(fab), intent(in) :: fb
    integer :: r(fb%dim)
    r = lwb(fb%pbx)
  end function fab_plwb
  pure function ifab_plwb(fb) result(r)
    type(ifab), intent(in) :: fb
    integer :: r(fb%dim)
    r = lwb(fb%pbx)
  end function ifab_plwb

  pure function fab_pupb(fb) result(r)
    type(fab), intent(in) :: fb
    integer :: r(fb%dim)
    r = upb(fb%pbx)
  end function fab_pupb
  pure function ifab_pupb(fb) result(r)
    type(ifab), intent(in) :: fb
    integer :: r(fb%dim)
    r = upb(fb%pbx)
  end function ifab_pupb

  pure function fab_ncomp(fb) result(r)
    integer :: r
    type(fab), intent(in) :: fb
    r = fb%nc
  end function fab_ncomp
  pure function zfab_ncomp(fb) result(r)
    integer :: r
    type(zfab), intent(in) :: fb
    r = fb%nc
  end function zfab_ncomp
  pure function ifab_ncomp(fb) result(r)
    integer :: r
    type(ifab), intent(in) :: fb
    r = fb%nc
  end function ifab_ncomp
  pure function lfab_ncomp(fb) result(r)
    integer :: r
    type(lfab), intent(in) :: fb
    r = fb%nc
  end function lfab_ncomp

  pure function fab_get_box(fb) result(r)
    type(fab), intent(in) :: fb
    type(box) :: r
    r = fb%bx
  end function fab_get_box
  pure function zfab_get_box(fb) result(r)
    type(zfab), intent(in) :: fb
    type(box) :: r
    r = fb%bx
  end function zfab_get_box
  pure function ifab_get_box(fb) result(r)
    type(ifab), intent(in) :: fb
    type(box) :: r
    r = fb%bx
  end function ifab_get_box
  pure function lfab_get_box(fb) result(r)
    type(lfab), intent(in) :: fb
    type(box) :: r
    r = fb%bx
  end function lfab_get_box

  pure function fab_get_pbox(fb) result(r)
    type(fab), intent(in) :: fb
    type(box) :: r
    r = fb%pbx
  end function fab_get_pbox
  pure function zfab_get_pbox(fb) result(r)
    type(zfab), intent(in) :: fb
    type(box) :: r
    r = fb%pbx
  end function zfab_get_pbox
  pure function ifab_get_pbox(fb) result(r)
    type(ifab), intent(in) :: fb
    type(box) :: r
    r = fb%pbx
  end function ifab_get_pbox
  pure function lfab_get_pbox(fb) result(r)
    type(lfab), intent(in) :: fb
    type(box) :: r
    r = fb%pbx
  end function lfab_get_pbox

  pure function fab_get_ibox(fb) result(r)
    type(fab), intent(in) :: fb
    type(box) :: r
    r = fb%ibx
  end function fab_get_ibox
  pure function zfab_get_ibox(fb) result(r)
    type(zfab), intent(in) :: fb
    type(box) :: r
    r = fb%ibx
  end function zfab_get_ibox
  pure function ifab_get_ibox(fb) result(r)
    type(ifab), intent(in) :: fb
    type(box) :: r
    r = fb%ibx
  end function ifab_get_ibox
  pure function lfab_get_ibox(fb) result(r)
    type(lfab), intent(in) :: fb
    type(box) :: r
    r = fb%ibx
  end function lfab_get_ibox
  
  pure function fab_built_q(fb) result(r)
    type(fab), intent(in) :: fb
    logical :: r
    r = fb%dim /= 0
  end function fab_built_q
  pure function zfab_built_q(fb) result(r)
    type(zfab), intent(in) :: fb
    logical :: r
    r = fb%dim /= 0
  end function zfab_built_q
  pure function ifab_built_q(fb) result(r)
    type(ifab), intent(in) :: fb
    logical :: r
    r = fb%dim /= 0
  end function ifab_built_q
  pure function lfab_built_q(fb) result(r)
    type(lfab), intent(in) :: fb
    logical :: r
    r = fb%dim /= 0
  end function lfab_built_q

  subroutine fab_build(fb, bx, nc, ng, nodal, alloc, stencil)
    type(fab), intent(out) :: fb
    type(box), intent(in)  :: bx
    integer, intent(in), optional :: ng, nc
    logical, intent(in), optional :: nodal(:)
    logical, intent(in), optional :: alloc
    logical, intent(in), optional :: stencil
    integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
    integer :: lnc, lng
    logical :: lal, lst
    lng = 0; if ( present(ng) ) lng = ng
    lnc = 1; if ( present(nc) ) lnc = nc
    lal = .true. ; if ( present(alloc)   ) lal = alloc
    lst = .false.; if ( present(stencil) ) lst = stencil
    lo = 1
    hi = 1
    fb%dim = bx%dim
    fb%bx = bx
    fb%nc = lnc
    fb%ibx = box_nodalize(bx, nodal)
    fb%pbx = grow(fb%ibx, lng)
    lo(1:fb%dim) = fb%pbx%lo(1:fb%dim)
    hi(1:fb%dim) = fb%pbx%hi(1:fb%dim)
    if ( lal ) then
       if ( lst) then
          allocate(fb%p(1:lnc,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
       else
          allocate(fb%p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:lnc))
       end if
       if ( do_init_fabs ) call setval(fb, fab_default_init)
       call mem_stats_alloc(fab_ms, volume(fb, all=.TRUE.))
       if ( (fab_ms%num_alloc-fab_ms%num_dealloc) > fab_high_water_mark ) then
          fab_high_water_mark = (fab_ms%num_alloc-fab_ms%num_dealloc)
       end if
    end if
  end subroutine fab_build

  subroutine zfab_build(fb, bx, nc, ng, nodal, alloc)
    type(zfab), intent(out) :: fb
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
    fb%nc = lnc
    fb%ibx = box_nodalize(bx, nodal)
    fb%pbx = grow(fb%ibx, lng)
    lo(1:fb%dim) = fb%pbx%lo(1:fb%dim)
    hi(1:fb%dim) = fb%pbx%hi(1:fb%dim)
    if ( lal ) then
       allocate(fb%p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:lnc))
       if ( do_init_fabs ) call setval(fb, zfab_default_init)
       call mem_stats_alloc(zfab_ms, volume(fb, all=.TRUE.))
    end if
  end subroutine zfab_build

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
    fb%nc = lnc
    fb%ibx = box_nodalize(bx, nodal)
    fb%pbx = grow(fb%ibx, lng)
    lo(1:fb%dim) = fb%pbx%lo(1:fb%dim)
    hi(1:fb%dim) = fb%pbx%hi(1:fb%dim)
    if ( lal ) then
       allocate(fb%p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:lnc))
       if ( do_init_fabs ) call setval(fb, ifab_default_init)
       call mem_stats_alloc(ifab_ms, volume(fb, all=.TRUE.))
    end if
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
    fb%nc = lnc
    fb%ibx = box_nodalize(bx, nodal)
    fb%pbx = grow(fb%ibx, lng)
    lo(1:fb%dim) = fb%pbx%lo(1:fb%dim)
    hi(1:fb%dim) = fb%pbx%hi(1:fb%dim)
    if ( lal ) then
       allocate(fb%p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:lnc))
       if ( do_init_fabs ) call setval(fb, lfab_default_init)
       call mem_stats_alloc(lfab_ms, volume(fb, all=.TRUE.))
    end if
  end subroutine lfab_build

  subroutine fab_destroy(fb)
    type(fab), intent(inout) :: fb
    if ( associated(fb%p) ) then
       deallocate(fb%p)
       call mem_stats_dealloc(fab_ms, volume(fb, all=.TRUE.))
    end if
    fb%bx  = nobox(fb%dim)
    fb%dim = 0
    fb%nc  = 0
  end subroutine fab_destroy
  subroutine zfab_destroy(fb)
    type(zfab), intent(inout) :: fb
    if ( associated(fb%p) ) then
       deallocate(fb%p)
       call mem_stats_dealloc(zfab_ms, volume(fb, all=.TRUE.))
    end if
    fb%bx  = nobox(fb%dim)
    fb%dim = 0
    fb%nc  = 0
  end subroutine zfab_destroy
  subroutine ifab_destroy(fb)
    type(ifab), intent(inout) :: fb
    if ( associated(fb%p) ) then
       deallocate(fb%p)
       call mem_stats_dealloc(ifab_ms, volume(fb, all=.TRUE.))
    end if
    fb%bx  = nobox(fb%dim)
    fb%dim = 0
    fb%nc  = 0
  end subroutine ifab_destroy
  subroutine lfab_destroy(fb)
    type(lfab), intent(inout) :: fb
    if ( associated(fb%p) ) then
       deallocate(fb%p)
       call mem_stats_dealloc(lfab_ms, volume(fb, all=.TRUE.))
    end if
    fb%bx  = nobox(fb%dim)
    fb%dim = 0
    fb%nc  = 0
  end subroutine lfab_destroy

  pure function fab_dim(fb) result(r)
    type(fab), intent(in) :: fb
    integer :: r
    r = fb%dim
  end function fab_dim
  pure function ifab_dim(fb) result(r)
    type(ifab), intent(in) :: fb
    integer :: r
    r = fb%dim
  end function ifab_dim
  pure function lfab_dim(fb) result(r)
    type(lfab), intent(in) :: fb
    integer :: r
    r = fb%dim
  end function lfab_dim
  pure function zfab_dim(fb) result(r)
    type(zfab), intent(in) :: fb
    integer :: r
    r = fb%dim
  end function zfab_dim

  function fab_dataptr(fb) result(r)
    type(fab), intent(in) :: fb
    real(kind=dp_t), pointer :: r(:,:,:,:)
    r => fb%p
  end function fab_dataptr
  function zfab_dataptr(fb) result(r)
    type(zfab), intent(in) :: fb
    complex(kind=dp_t), pointer :: r(:,:,:,:)
    r => fb%p
  end function zfab_dataptr
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
    use bl_error_module
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
  function zfab_dataptr_bx(fb, bx) result(r)
    use bl_error_module
    type(zfab), intent(in) :: fb
    type(box), intent(in) :: bx
    complex(kind=dp_t), pointer :: r(:,:,:,:)
    if ( .not. contains(fb%pbx, bx) ) call bl_error('fab_dataptr_bx: bx is too large')
    select case (fb%dim)
    case (1)
       r => fb%p(bx%lo(1):bx%hi(1),:,:,:)
    case (2)
       r => fb%p(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),:,:)
    case (3)
       r => fb%p(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),:)
    end select
  end function zfab_dataptr_bx
  function ifab_dataptr_bx(fb, bx) result(r)
    use bl_error_module
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
    use bl_error_module
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

  function fab_dataptr_c(fb, c, nc) result(r)
    use bl_error_module
    type(fab), intent(in) :: fb
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    real(kind=dp_t), pointer :: r(:,:,:,:)
    integer :: lnc
    lnc = 1; if ( present(nc) ) lnc = nc
    if ( (c+lnc-1) > fb%nc ) call bl_error('fab_dataptr_c: not enough components')
    r => fb%p(:,:,:,c:c+lnc-1)
  end function fab_dataptr_c
  function zfab_dataptr_c(fb, c, nc) result(r)
    use bl_error_module
    type(zfab), intent(in) :: fb
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    complex(kind=dp_t), pointer :: r(:,:,:,:)
    integer :: lnc
    lnc = 1; if ( present(nc) ) lnc = nc
    if ( (c+lnc-1) > fb%nc ) call bl_error('fab_dataptr_c: not enough components')
    r => fb%p(:,:,:,c:c+lnc-1)
  end function zfab_dataptr_c
  function ifab_dataptr_c(fb, c, nc) result(r)
    use bl_error_module
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
    use bl_error_module
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
    use bl_error_module
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
  function zfab_dataptr_bx_c(fb, bx, c, nc) result(r)
    use bl_error_module
    type(zfab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    complex(kind=dp_t), pointer :: r(:,:,:,:)
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
  end function zfab_dataptr_bx_c
  function ifab_dataptr_bx_c(fb, bx, c, nc) result(r)
    use bl_error_module
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
    use bl_error_module
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
    use bl_error_module
    type(fab), intent(inout) :: fb
    real(kind=dp_t), intent(in) :: val
    if ( .not. associated(fb%p) ) call bl_error("FAB_SETVAL: not associated")
    fb%p = val
  end subroutine fab_setval
  subroutine zfab_setval(fb, val)
    use bl_error_module
    type(zfab), intent(inout) :: fb
    complex(kind=dp_t), intent(in) :: val
    if ( .not. associated(fb%p) ) call bl_error("FAB_SETVAL: not associated")
    fb%p = val
  end subroutine zfab_setval
  subroutine ifab_setval(fb, val)
    use bl_error_module
    type(ifab), intent(inout) :: fb
    integer, intent(in) :: val
    if ( .not. associated(fb%p) ) call bl_error("FAB_SETVAL: not associated")
    fb%p = val
  end subroutine ifab_setval
  subroutine lfab_setval(fb, val)
    use bl_error_module
    type(lfab), intent(inout) :: fb
    logical, intent(in) :: val
    if ( .not. associated(fb%p) ) call bl_error("FAB_SETVAL: not associated")
    fb%p = val
  end subroutine lfab_setval

  subroutine fab_setval_bx(fb, val, bx)
    use bl_error_module
    type(fab), intent(inout) :: fb
    type(box), intent(in) :: bx
    real(kind=dp_t), intent(in) :: val
    real(kind=dp_t), pointer :: p(:,:,:,:)
    p => fab_dataptr_bx(fb, bx)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine fab_setval_bx
  subroutine zfab_setval_bx(fb, val, bx)
    use bl_error_module
    type(zfab), intent(inout) :: fb
    type(box), intent(in) :: bx
    complex(kind=dp_t), intent(in) :: val
    complex(kind=dp_t), pointer :: p(:,:,:,:)
    p => zfab_dataptr_bx(fb, bx)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine zfab_setval_bx
  subroutine ifab_setval_bx(fb, val, bx)
    use bl_error_module
    type(ifab), intent(inout) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: val
    integer, pointer :: p(:,:,:,:)
    p => ifab_dataptr_bx(fb, bx)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine ifab_setval_bx
  subroutine lfab_setval_bx(fb, val, bx)
    use bl_error_module
    type(lfab), intent(inout) :: fb
    type(box), intent(in) :: bx
    logical, intent(in) :: val
    logical, pointer :: p(:,:,:,:)
    p => lfab_dataptr_bx(fb, bx)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine lfab_setval_bx

  subroutine fab_setval_c(fb, val, c, nc)
    use bl_error_module
    type(fab), intent(inout) :: fb
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    real(kind=dp_t), intent(in) :: val
    real(kind=dp_t), pointer :: p(:,:,:,:)
    p => fab_dataptr_c(fb, c, nc)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine fab_setval_c
  subroutine zfab_setval_c(fb, val, c, nc)
    use bl_error_module
    type(zfab), intent(inout) :: fb
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    complex(kind=dp_t), intent(in) :: val
    complex(kind=dp_t), pointer :: p(:,:,:,:)
    p => zfab_dataptr_c(fb, c, nc)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine zfab_setval_c
  subroutine ifab_setval_c(fb, val, c, nc)
    use bl_error_module
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
    use bl_error_module
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
    use bl_error_module
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
  subroutine zfab_setval_bx_c(fb, val, bx, c, nc)
    use bl_error_module
    type(zfab), intent(inout) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    complex(kind=dp_t), intent(in) :: val
    complex(kind=dp_t), pointer :: p(:,:,:,:)
    p => zfab_dataptr_bx_c(fb, bx, c, nc)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine zfab_setval_bx_c
  subroutine ifab_setval_bx_c(fb, val, bx, c, nc)
    use bl_error_module
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
    use bl_error_module
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

  subroutine fab_print(fb, comp, str, unit, all, data, bx, skip)
    use bl_IO_module
    type(fab), intent(in) :: fb
    integer, intent(in), optional :: comp
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
       write(unit=un, fmt=*) 'NOT ASSOCIATED'
    else
       if ( present(comp) ) then
        select case (fb%dim)
        case (1)
          call print_1d(fb%p(:,1,1,comp:comp), lbound(fb%p), intersection(fb%ibx,lbx))
        case (2)
          call print_2d(fb%p(:,:,1,comp:comp), lbound(fb%p), intersection(fb%ibx,lbx))
        case (3)
          call print_3d(fb%p(:,:,:,comp:comp), lbound(fb%p), intersection(fb%ibx,lbx))
        end select
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
       write(unit=un, fmt=*) 'NOT ASSOCIATED'
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
       write(unit=un, fmt=*) 'NOT ASSOCIATED'
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
  subroutine zfab_print(fb, str, unit, all, data, bx, skip)
    use bl_IO_module
    type(zfab), intent(in) :: fb
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
    write(unit=un, fmt='("ZFAB")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(": ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" DIM     = ",i2)') fb%dim
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
       write(unit=un, fmt=*) 'NOT ASSOCIATED'
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
      complex(kind=dp_t), intent(in) :: fb(lo(1):,:)
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
               write(unit=un, fmt='(A1,1X,I3,1(1X,I5),1X,2G25.15)') &
                    c, n, i, fb(i,n)
            end do
         end do
      end if
    end subroutine print_1d
    subroutine print_2d(fb, lo, bx)
      integer, intent(in) :: lo(:)
      type(box), intent(in) :: bx
      complex(kind=dp_t), intent(in) :: fb(lo(1):,lo(2):,:)
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
                  write(unit=un, fmt='(A1,1X,I3,2(1X,I5),1X,2G25.15)') &
                       c, n, i, j, fb(i,j,n)
               end do
            end do
         end do
      end if
    end subroutine print_2d
    subroutine print_3d(fb, lo, bx)
      integer, intent(in) :: lo(:)
      type(box), intent(in) :: bx
      complex(kind=dp_t), intent(in) :: fb(lo(1):,lo(2):,lo(3):,:)
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
                     write(unit=un, fmt='(A1,1X,I3,3(1X,I5),1X,2G25.15)') &
                          c, n, i, j, k, fb(i,j,k,n)
                  end do
               end do
            end do
         end do
      end if
    end subroutine print_3d
  end subroutine zfab_print

  function fab_minval_doit(ap) result(r)

    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t)          :: r, r1

    integer :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! minval(ap)

    r1 = Huge(r)

    !$OMP PARALLEL PRIVATE(i,j,k,n) REDUCTION(MIN : r1) IF((hi(3)-lo(3)).ge.7)
    do n = lo(4), hi(4)
       !$OMP DO
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                r1 = min(r1,ap(i,j,k,n))
             end do
          end do
       end do
       !$OMP END DO NOWAIT
    end do
    !$OMP END PARALLEL

    r = r1

  end function fab_minval_doit

  function fab_maxval_doit(ap) result(r)

    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t)          :: r, r1

    integer :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! maxval(ap)

    r1 = -Huge(r)

    !$OMP PARALLEL PRIVATE(i,j,k,n) REDUCTION(MAX : r1) IF((hi(3)-lo(3)).ge.7)
    do n = lo(4), hi(4)
       !$OMP DO
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                r1 = max(r1,ap(i,j,k,n))
             end do
          end do
       end do
       !$OMP END DO NOWAIT
    end do
    !$OMP END PARALLEL

    r = r1

  end function fab_maxval_doit

  function fab_max_val(fb, all) result(r)
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
    r = fab_maxval_doit(mp)
  end function fab_max_val
  function fab_max_val_c(fb, c, nc, all) result(r)
    real(kind=dp_t) :: r
    type(fab), intent(in) :: fb
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    logical :: lall
    real(dp_t), pointer :: mp(:,:,:,:)
    lall = .FALSE.; if ( present(all) ) lall = all
    if ( lall ) then
       mp => dataptr(fb, get_pbox(fb), c, nc)
    else
       mp => dataptr(fb, get_ibox(fb), c, nc)
    end if
    r = fab_maxval_doit(mp)
  end function fab_max_val_c
  function fab_max_val_bx(fb, bx) result(r)
    real(kind=dp_t) :: r
    type(fab), intent(in) :: fb
    type(box), intent(in) :: bx
    real(dp_t), pointer :: mp(:,:,:,:)
    type(box) :: sbx
    sbx = intersection(bx, get_ibox(fb))
    mp => dataptr(fb, sbx)
    r = fab_maxval_doit(mp)
  end function fab_max_val_bx
  function fab_max_val_bx_c(fb, bx, c, nc) result(r)
    real(kind=dp_t) :: r
    type(fab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    real(dp_t), pointer :: mp(:,:,:,:)
    type(box) :: sbx
    sbx = intersection(bx, get_ibox(fb))
    mp => dataptr(fb, sbx, c, nc)
    r = fab_maxval_doit(mp)
  end function fab_max_val_bx_c
  function ifab_max_val(fb, all) result(r)
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
  end function ifab_max_val
  function ifab_max_val_c(fb, c, nc, all) result(r)
    integer :: r
    type(ifab), intent(in) :: fb
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    logical :: lall
    integer, pointer :: mp(:,:,:,:)
    lall = .FALSE.; if ( present(all) ) lall = all
    if ( lall ) then
       mp => dataptr(fb, get_pbox(fb), c, nc)
    else
       mp => dataptr(fb, get_ibox(fb), c, nc)
    end if
    r = maxval(mp)
  end function ifab_max_val_c
  function ifab_max_val_bx(fb, bx) result(r)
    integer :: r
    type(ifab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, pointer :: mp(:,:,:,:)
    type(box) :: sbx
    sbx = intersection(bx, get_ibox(fb))
    mp => dataptr(fb, sbx)
    r = maxval(mp)
  end function ifab_max_val_bx
  function ifab_max_val_bx_c(fb, bx, c, nc) result(r)
    integer :: r
    type(ifab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    integer, pointer :: mp(:,:,:,:)
    type(box) :: sbx
    sbx = intersection(bx, get_ibox(fb))
    mp => dataptr(fb, sbx, c, nc)
    r = maxval(mp)
  end function ifab_max_val_bx_c

  function fab_min_val(fb, all) result(r)
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
    r = fab_minval_doit(mp)
  end function fab_min_val
  function fab_min_val_c(fb, c, nc, all) result(r)
    real(kind=dp_t) :: r
    type(fab), intent(in) :: fb
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    logical :: lall
    real(dp_t), pointer :: mp(:,:,:,:)
    lall = .FALSE.; if ( present(all) ) lall = all
    if ( lall ) then
       mp => dataptr(fb, get_pbox(fb), c, nc)
    else
       mp => dataptr(fb, get_ibox(fb), c, nc)
    end if
    r = fab_minval_doit(mp)
  end function fab_min_val_c
  function fab_min_val_bx(fb, bx) result(r)
    real(kind=dp_t) :: r
    type(fab), intent(in) :: fb
    type(box), intent(in) :: bx
    real(dp_t), pointer :: mp(:,:,:,:)
    type(box) :: sbx
    sbx = intersection(bx, get_ibox(fb))
    mp => dataptr(fb, sbx)
    r = fab_minval_doit(mp)
  end function fab_min_val_bx
  function fab_min_val_bx_c(fb, bx, c, nc) result(r)
    real(kind=dp_t) :: r
    type(fab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    real(dp_t), pointer :: mp(:,:,:,:)
    type(box) :: sbx
    sbx = intersection(bx, get_ibox(fb))
    mp => dataptr(fb, sbx, c, nc)
    r = fab_minval_doit(mp)
  end function fab_min_val_bx_c
  function ifab_min_val(fb, all) result(r)
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
  end function ifab_min_val
  function ifab_min_val_c(fb, c, nc, all) result(r)
    integer :: r
    type(ifab), intent(in) :: fb
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    logical :: lall
    integer, pointer :: mp(:,:,:,:)
    lall = .FALSE.; if ( present(all) ) lall = all
    if ( lall ) then
       mp => dataptr(fb, get_pbox(fb), c, nc)
    else
       mp => dataptr(fb, get_ibox(fb), c, nc)
    end if
    r = minval(mp)
  end function ifab_min_val_c
  function ifab_min_val_bx(fb, bx) result(r)
    integer :: r
    type(ifab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, pointer :: mp(:,:,:,:)
    type(box) :: sbx
    sbx = intersection(bx, get_ibox(fb))
    mp => dataptr(fb, sbx)
    r = minval(mp)
  end function ifab_min_val_bx
  function ifab_min_val_bx_c(fb, bx, c, nc) result(r)
    integer :: r
    type(ifab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    integer, pointer :: mp(:,:,:,:)
    type(box) :: sbx
    sbx = intersection(bx, get_ibox(fb))
    mp => dataptr(fb, sbx, c, nc)
    r = minval(mp)
  end function ifab_min_val_bx_c

  function lfab_count(fb, all) result(r)
    integer :: r
    type(lfab), intent(in) :: fb
    logical, intent(in), optional :: all
    logical, pointer :: lp(:,:,:,:)
    logical :: lall
    lall = .false. ; if ( present(all) ) lall = all
    if ( lall ) then
       lp => dataptr(fb, get_pbox(fb))
    else
       lp => dataptr(fb, get_ibox(fb))
    end if
    r = count(lp)
  end function lfab_count

end module fab_module
