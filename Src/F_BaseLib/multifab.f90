module multifab_module

  use layout_module
  use fab_module
  use bl_mem_stat_module

  implicit none

  type multifab
     logical :: bound = .false.
     integer :: dim = 0
     integer :: nboxes = 0
     integer :: nc = 1
     integer :: ng = 0
     logical, pointer :: nodal(:) => Null()
     type(layout) :: la
     type(fab), pointer :: fbs(:) => Null()
  end type multifab

  type imultifab
     logical :: bound = .false.
     integer :: dim = 0
     integer :: nboxes = 0
     integer :: nc = 1
     integer :: ng = 0
     logical, pointer :: nodal(:) => Null()
     type(layout) :: la
     type(ifab), pointer :: fbs(:) => Null()
  end type imultifab

  type lmultifab
     logical :: bound = .false.
     integer :: dim = 0
     integer :: nboxes = 0
     integer :: nc = 1
     integer :: ng = 0
     logical, pointer :: nodal(:) => Null()
     type(layout) :: la
     type(lfab), pointer :: fbs(:) => Null()
  end type lmultifab

  interface cell_centered_q
     module procedure multifab_cell_centered_q
     module procedure imultifab_cell_centered_q
     module procedure lmultifab_cell_centered_q
  end interface

  interface nodal_q
     module procedure multifab_nodal_q
     module procedure imultifab_nodal_q
     module procedure lmultifab_nodal_q
  end interface

  interface built_q
     module procedure multifab_built_q
     module procedure imultifab_built_q
     module procedure lmultifab_built_q
  end interface

  interface build
     module procedure multifab_build
     module procedure multifab_build_copy

     module procedure imultifab_build
     module procedure imultifab_build_copy

     module procedure lmultifab_build
  end interface

  interface destroy
     module procedure multifab_destroy
     module procedure imultifab_destroy
     module procedure lmultifab_destroy
  end interface

  interface copy
     module procedure multifab_copy_c
     module procedure multifab_copy
     module procedure imultifab_copy_c
     module procedure imultifab_copy
  end interface

  interface conformant_q
     module procedure multifab_conformant_q
  end interface

  interface dataptr
     module procedure multifab_dataptr
     module procedure multifab_dataptr_bx
     module procedure multifab_dataptr_bx_c
     module procedure multifab_dataptr_c

     module procedure imultifab_dataptr
     module procedure imultifab_dataptr_bx
     module procedure imultifab_dataptr_bx_c
     module procedure imultifab_dataptr_c

     module procedure lmultifab_dataptr
     module procedure lmultifab_dataptr_bx
     module procedure lmultifab_dataptr_bx_c
     module procedure lmultifab_dataptr_c
  end interface

  interface print
     module procedure multifab_print
     module procedure imultifab_print
  end interface

  interface setval
     module procedure multifab_setval
     module procedure multifab_setval_bx
     module procedure multifab_setval_c
     module procedure multifab_setval_bx_c

     module procedure imultifab_setval
     module procedure imultifab_setval_bx
     module procedure imultifab_setval_c
     module procedure imultifab_setval_bx_c

     module procedure lmultifab_setval
     module procedure lmultifab_setval_bx
     module procedure lmultifab_setval_c
     module procedure lmultifab_setval_bx_c

     module procedure lmultifab_setval_ba
  end interface

  interface get_boxarray
     module procedure multifab_get_boxarray
     module procedure imultifab_get_boxarray
     module procedure lmultifab_get_boxarray
  end interface

  interface get_pbox
     module procedure multifab_get_pbox
     module procedure imultifab_get_pbox
     module procedure lmultifab_get_pbox
  end interface

  interface get_ibox
     module procedure multifab_get_ibox
     module procedure imultifab_get_ibox
     module procedure lmultifab_get_ibox
  end interface
  
  interface get_box
     module procedure multifab_get_box
     module procedure imultifab_get_box
     module procedure lmultifab_get_box
  end interface

  interface nboxes
     module procedure multifab_nboxes
     module procedure imultifab_nboxes
  end interface

  interface volume
     module procedure multifab_volume
     module procedure imultifab_volume
     module procedure lmultifab_volume
  end interface

  interface local
     module procedure multifab_local
     module procedure imultifab_local
     module procedure lmultifab_local
  end interface

  interface remote
     module procedure multifab_remote
     module procedure imultifab_remote
     module procedure lmultifab_remote
  end interface

  interface get_layout
     module procedure multifab_get_layout
     module procedure imultifab_get_layout
     module procedure lmultifab_get_layout
  end interface

  interface rescale
     module procedure multifab_rescale_c
     module procedure multifab_rescale
  end interface

  interface saxpy
     module procedure multifab_saxpy_3_c
     module procedure multifab_saxpy_3
     module procedure multifab_saxpy_4
     module procedure multifab_saxpy_5
  end interface

  interface fill_boundary
     module procedure multifab_fill_boundary
     module procedure multifab_fill_boundary_c
     module procedure imultifab_fill_boundary
     module procedure lmultifab_fill_boundary
  end interface

  interface internal_sync
     module procedure multifab_internal_sync
     module procedure lmultifab_internal_sync
  end interface

  interface ncomp
     module procedure multifab_ncomp
     module procedure imultifab_ncomp
     module procedure lmultifab_ncomp
  end interface

  interface norm_l2
     module procedure multifab_norm_l2_c
     module procedure multifab_norm_l2
  end interface

  interface norm_inf
     module procedure multifab_norm_inf_c
     module procedure multifab_norm_inf
  end interface

  interface inorm_inf
     module procedure imultifab_norm_inf_c
     module procedure imultifab_norm_inf
  end interface

  interface norm_l1
     module procedure multifab_norm_l1_c
     module procedure multifab_norm_l1
  end interface

  interface dot
     module procedure multifab_dot
     module procedure multifab_dot_c
     module procedure multifab_dot_cc
  end interface

  interface plus_plus
     module procedure multifab_plus_plus
     module procedure multifab_plus_plus_s
     module procedure multifab_plus_plus_c
     module procedure multifab_plus_plus_s_c
  end interface

  interface div_div
     module procedure multifab_div_div
     module procedure multifab_div_div_s
     module procedure multifab_div_div_c
     module procedure multifab_div_div_s_c
  end interface

  interface mult_mult
     module procedure multifab_mult_mult
     module procedure multifab_mult_mult_s
     module procedure multifab_mult_mult_c
     module procedure multifab_mult_mult_s_c
  end interface

  interface sub_sub
     module procedure multifab_sub_sub
     module procedure multifab_sub_sub_s
     module procedure multifab_sub_sub_c
     module procedure multifab_sub_sub_s_c
  end interface

  type(mem_stats), private, save ::  multifab_ms
  type(mem_stats), private, save :: imultifab_ms
  type(mem_stats), private, save :: lmultifab_ms

! interface assignment(=)
!    module procedure multifab_die_die
! end interface

  private :: build_nodal_dot_mask

  ! buffers used for parallel fill_boundary
  real(kind=dp_t), allocatable, save, private  :: g_snd_d(:), g_rcv_d(:)
  integer, allocatable, save, private  :: g_snd_i(:), g_rcv_i(:)
  logical, allocatable, save, private  :: g_snd_l(:), g_rcv_l(:)

  logical, private :: d_fb_fancy = .true.
  logical, private :: i_fb_fancy = .true.
  logical, private :: l_fb_fancy = .true.

  logical, private :: d_fb_async = .false. ! Do both recv's and send's asynchronously?
  logical, private :: i_fb_async = .false. ! Do both recv's and send's asynchronously?
  logical, private :: l_fb_async = .false. ! Do both recv's and send's asynchronously?

contains

  subroutine multifab_get_behavior(fb_async, fb_fancy)
    logical, intent(out), optional :: fb_async
    logical, intent(out), optional :: fb_fancy
    if ( present(fb_async) ) fb_async = d_fb_async
    if ( present(fb_fancy) ) fb_fancy = d_fb_fancy
  end subroutine multifab_get_behavior
  subroutine imultifab_get_behavior(fb_async, fb_fancy)
    logical, intent(out), optional :: fb_async
    logical, intent(out), optional :: fb_fancy
    if ( present(fb_async) ) fb_async = i_fb_async
    if ( present(fb_fancy) ) fb_fancy = i_fb_fancy
  end subroutine imultifab_get_behavior
  subroutine lmultifab_get_behavior(fb_async, fb_fancy)
    logical, intent(out), optional :: fb_async
    logical, intent(out), optional :: fb_fancy
    if ( present(fb_async) ) fb_async = l_fb_async
    if ( present(fb_fancy) ) fb_fancy = l_fb_fancy
  end subroutine lmultifab_get_behavior

  subroutine multifab_set_behavior(fb_async, fb_fancy)
    logical, intent(in), optional :: fb_async
    logical, intent(in), optional :: fb_fancy
    if ( present(fb_async) ) d_fb_async = fb_async
    if ( present(fb_fancy) ) d_fb_fancy = fb_fancy
  end subroutine multifab_set_behavior
  subroutine imultifab_set_behavior(fb_async, fb_fancy)
    logical, intent(in), optional :: fb_async
    logical, intent(in), optional :: fb_fancy
    if ( present(fb_async) ) i_fb_async = fb_async
    if ( present(fb_fancy) ) i_fb_fancy = fb_fancy
  end subroutine imultifab_set_behavior
  subroutine lmultifab_set_behavior(fb_async, fb_fancy)
    logical, intent(in), optional :: fb_async
    logical, intent(in), optional :: fb_fancy
    if ( present(fb_async) ) l_fb_async = fb_async
    if ( present(fb_fancy) ) l_fb_fancy = fb_fancy
  end subroutine lmultifab_set_behavior

! subroutine multifab_die_die(mf,mf1)
!   type(multifab), intent(in) :: mf1
!   type(multifab), intent(inout) :: mf
!   call bl_error("DIE")
! end subroutine multifab_die_die

  function multifab_ncomp(mf) result(r)
    integer :: r
    type(multifab), intent(in) :: mf
    r = mf%nc
  end function multifab_ncomp
  function imultifab_ncomp(mf) result(r)
    integer :: r
    type(imultifab), intent(in) :: mf
    r = mf%nc
  end function imultifab_ncomp
  function lmultifab_ncomp(mf) result(r)
    integer :: r
    type(lmultifab), intent(in) :: mf
    r = mf%nc
  end function lmultifab_ncomp

  subroutine multifab_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    multifab_ms = ms
  end subroutine multifab_set_mem_stats
  subroutine imultifab_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    imultifab_ms = ms
  end subroutine imultifab_set_mem_stats
  subroutine lmultifab_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    lmultifab_ms = ms
  end subroutine lmultifab_set_mem_stats

  function multifab_mem_stats() result(r)
    type(mem_stats) :: r
    r = multifab_ms
  end function multifab_mem_stats
  function imultifab_mem_stats() result(r)
    type(mem_stats) :: r
    r = imultifab_ms
  end function imultifab_mem_stats
  function lmultifab_mem_stats() result(r)
    type(mem_stats) :: r
    r = lmultifab_ms
  end function lmultifab_mem_stats

  subroutine check_conformance(a,b)
    type(multifab), intent(in) :: a, b
    if ( not_equal(a%la, b%la) ) call bl_error("MULTIFAB::check_conformance:: not conformant")
  end subroutine check_conformance

  function multifab_cell_centered_q(mf) result(r)
    logical :: r
    type(multifab), intent(in) :: mf
    r = .not. any(mf%nodal)
  end function multifab_cell_centered_q
  function imultifab_cell_centered_q(mf) result(r)
    logical :: r
    type(imultifab), intent(in) :: mf
    r = .not. any(mf%nodal)
  end function imultifab_cell_centered_q
  function lmultifab_cell_centered_q(mf) result(r)
    logical :: r
    type(lmultifab), intent(in) :: mf
    r = .not. any(mf%nodal)
  end function lmultifab_cell_centered_q
  
  function multifab_nodal_q(mf) result(r)
    logical :: r
    type(multifab), intent(in) :: mf
    r = all(mf%nodal)
  end function multifab_nodal_q
  function imultifab_nodal_q(mf) result(r)
    logical :: r
    type(imultifab), intent(in) :: mf
    r = all(mf%nodal)
  end function imultifab_nodal_q
  function lmultifab_nodal_q(mf) result(r)
    logical :: r
    type(lmultifab), intent(in) :: mf
    r = mf%dim /= 0
  end function lmultifab_nodal_q
  
  function multifab_built_q(mf) result(r)
    logical :: r
    type(multifab), intent(in) :: mf
    r = mf%dim /= 0
  end function multifab_built_q
  function imultifab_built_q(mf) result(r)
    logical :: r
    type(imultifab), intent(in) :: mf
    r = mf%dim /= 0
  end function imultifab_built_q
  function lmultifab_built_q(mf) result(r)
    logical :: r
    type(lmultifab), intent(in) :: mf
    r = mf%dim /= 0
  end function lmultifab_built_q
  
  function multifab_conformant_q(a,b) result(r)
    logical :: r
    type(multifab), intent(in) :: a, b
    r = equal(a%la, b%la) .and. a%dim == b%dim .and. a%ng == b%ng .and. a%nc == b%nc
  end function multifab_conformant_q
  function imultifab_conformant_q(a,b) result(r)
    logical :: r
    type(imultifab), intent(in) :: a, b
    r = equal(a%la, b%la) .and. a%dim == b%dim .and. a%ng == b%ng .and. a%nc == b%nc
  end function imultifab_conformant_q

  function multifab_get_layout(mf) result(r)
    type(layout) :: r
    type(multifab), intent(in) :: mf
    r = mf%la
  end function multifab_get_layout
  function imultifab_get_layout(mf) result(r)
    type(layout) :: r
    type(imultifab), intent(in) :: mf
    r = mf%la
  end function imultifab_get_layout
  function lmultifab_get_layout(mf) result(r)
    type(layout) :: r
    type(lmultifab), intent(in) :: mf
    r = mf%la
  end function lmultifab_get_layout

  subroutine multifab_build(mf, la, nc, ng, nodal)
    type(multifab), intent(out) :: mf
    type(layout), intent(in) :: la
    integer, intent(in), optional :: nc, ng
    logical, intent(in), optional :: nodal(:)
    integer :: i
    integer :: lnc, lng
    if ( built_q(mf) ) call bl_error("MULTIFAB_BUILD: already built")
    lng = 0; if ( present(ng) ) lng = ng
    lnc = 0; if ( present(nc) ) lnc = nc
    mf%dim = layout_dim(la)
    mf%la = la
    mf%nboxes = layout_nboxes(la)
    mf%nc = lnc
    mf%ng = lng
    allocate(mf%nodal(mf%dim))
    mf%nodal = .False.; if ( present(nodal) ) mf%nodal = nodal
    allocate(mf%fbs(mf%nboxes))
    do i = 1, mf%nboxes
      call fab_build( &
           mf%fbs(i), get_box(mf%la, i), &
           mf%nc, mf%ng, mf%nodal,  &
           alloc = multifab_local(mf, i))
    end do
    call mem_stats_alloc(multifab_ms, volume(mf, all = .TRUE.))
  end subroutine multifab_build

  subroutine imultifab_build(mf, la, nc, ng, nodal)
    type(imultifab), intent(out) :: mf
    type(layout), intent(in) :: la
    integer, intent(in), optional :: nc, ng
    logical, intent(in), optional :: nodal(:)
    integer :: i
    integer :: lnc, lng
    if ( built_q(mf) ) call bl_error("MULTIFAB_BUILD: already built")
    lng = 0; if ( present(ng) ) lng = ng
    lnc = 0; if ( present(nc) ) lnc = nc
    mf%dim = layout_dim(la)
    mf%la = la
    mf%nboxes = layout_nboxes(la)
    mf%nc = lnc
    mf%ng = lng
    allocate(mf%nodal(mf%dim))
    mf%nodal = .False.; if ( present(nodal) ) mf%nodal = nodal
    allocate(mf%fbs(mf%nboxes))
    do i = 1, mf%nboxes
       call ifab_build(mf%fbs(i), get_box(mf%la, i), &
            mf%nc, mf%ng, mf%nodal, &
            alloc = imultifab_local(mf, i))
    end do
    call mem_stats_alloc(imultifab_ms, volume(mf, all = .TRUE.))
  end subroutine imultifab_build

  subroutine lmultifab_build(mf, la, nc, ng, nodal)
    type(lmultifab), intent(out) :: mf
    type(layout), intent(in) :: la
    integer, intent(in),optional :: nc, ng
    logical, intent(in),optional :: nodal(:)
    integer :: i
    integer :: lnc, lng

    if ( built_q(mf) ) call bl_error("MULTIFAB_BUILD: already built")
    lng = 0; if ( present(ng) ) lng = ng
    lnc = 0; if ( present(nc) ) lnc = nc
    mf%dim = layout_dim(la)
    mf%la = la
    mf%nboxes = layout_nboxes(la)
    mf%nc = lnc
    mf%ng = lng
    allocate(mf%nodal(mf%dim))
    mf%nodal = .False.; if ( present(nodal) ) mf%nodal = nodal
    allocate(mf%fbs(mf%nboxes))
    do i = 1, mf%nboxes
       call lfab_build(mf%fbs(i), get_box(mf%la, i), &
            mf%nc, mf%ng, mf%nodal, &
            alloc = lmultifab_local(mf, i))
    end do
    call mem_stats_alloc(lmultifab_ms, volume(mf, all = .TRUE.))
  end subroutine lmultifab_build

  subroutine multifab_build_copy(m1, m2)
    type(multifab), intent(inout) :: m1
    type(multifab), intent(in) :: m2
    real(kind=dp_t), pointer :: m1p(:,:,:,:)
    real(kind=dp_t), pointer :: m2p(:,:,:,:)
    integer :: i
    if ( built_q(m1) ) call bl_error("MULTIFAB_BUILD_COPY: already built")
    if ( built_q(m1) ) call multifab_destroy(m1)
    m1%dim = m2%dim
    m1%la = m2%la
    m1%nboxes = m2%nboxes
    m1%nc = m2%nc
    m1%ng = m2%nc
    allocate(m1%nodal(m1%dim))
    m1%nodal = m2%nodal
    allocate(m1%fbs(m1%nboxes))
    do i = 1, m1%nboxes
       if ( multifab_remote(m1, i) ) cycle
       call fab_build(m1%fbs(i), get_box(m1%la, i), m1%nc, m1%ng, m1%nodal)
       m1p => dataptr(m1,i)
       m2p => dataptr(m2,i)
       m1p = m2p
    end do
    call mem_stats_alloc(multifab_ms, volume(m1, all = .TRUE.))
  end subroutine multifab_build_copy
  subroutine imultifab_build_copy(m1, m2)
    type(imultifab), intent(inout) :: m1
    type(imultifab), intent(in) :: m2
    integer, pointer :: m1p(:,:,:,:)
    integer, pointer :: m2p(:,:,:,:)
    integer :: i
    if ( built_q(m1) ) call bl_error("IMULTIFAB_BUILD_COPY: already built")
    if ( built_q(m1) ) call imultifab_destroy(m1)
    m1%dim = m2%dim
    m1%la = m2%la
    m1%nboxes = m2%nboxes
    m1%nc = m2%nc
    m1%ng = m2%nc
    allocate(m1%nodal(m1%dim))
    m1%nodal = m2%nodal
    allocate(m1%fbs(m1%nboxes))
    do i = 1, m1%nboxes
       if ( imultifab_remote(m1, i) ) cycle
       call ifab_build(m1%fbs(i), get_box(m1%la, i), m1%nc, m1%ng, m1%nodal)
       m1p => dataptr(m1,i)
       m2p => dataptr(m2,i)
       m1p = m2p
    end do
    call mem_stats_alloc(imultifab_ms, volume(m1, all = .TRUE.))
  end subroutine imultifab_build_copy

  subroutine multifab_destroy(mf)
    type(multifab), intent(inout) :: mf
    integer :: i
    if ( .not. mf%bound ) then
       call mem_stats_dealloc(multifab_ms, volume(mf, all = .TRUE.))
       do i = 1, mf%nboxes
          call fab_destroy(mf%fbs(i))
       end do
    end if
    deallocate(mf%fbs)
    deallocate(mf%nodal)
    mf%dim = 0
    mf%nc  = 0
    mf%ng  = 0
    mf%nboxes = 0
  end subroutine multifab_destroy
  subroutine imultifab_destroy(mf)
    type(imultifab), intent(inout) :: mf
    integer :: i
    if ( .not. mf%bound ) then
       call mem_stats_dealloc(imultifab_ms, volume(mf, all = .TRUE.))
       do i = 1, mf%nboxes
          call ifab_destroy(mf%fbs(i))
       end do
    end if
    deallocate(mf%fbs)
    deallocate(mf%nodal)
    mf%dim = 0
    mf%nc  = 0
    mf%ng  = 0
    mf%nboxes = 0
  end subroutine imultifab_destroy
  subroutine lmultifab_destroy(mf)
    type(lmultifab), intent(inout) :: mf
    integer :: i
    if ( .not. mf%bound ) then
       call mem_stats_dealloc(lmultifab_ms, volume(mf, all = .TRUE.))
       do i = 1, mf%nboxes
          call lfab_destroy(mf%fbs(i))
       end do
    end if
    deallocate(mf%fbs)
    deallocate(mf%nodal)
    mf%dim = 0
    mf%nc  = 0
    mf%ng  = 0
    mf%nboxes = 0
  end subroutine lmultifab_destroy

  function multifab_volume(mf, all) result(r)
    integer(kind=ll_t) :: r
    type(multifab), intent(in) :: mf
    logical, optional :: all
    integer :: i
    if ( all ) then
       r = 0_ll_t
       do i = 1, mf%nboxes
          r = r + volume(get_pbox(mf, i))
       end do
    else
       r = volume(get_boxarray(mf))
    end if
  end function multifab_volume
  function imultifab_volume(mf, all) result(r)
    integer(kind=ll_t) :: r
    type(imultifab), intent(in) :: mf
    logical, optional :: all
    integer :: i
    if ( all ) then
       r = 0_ll_t
       do i = 1, mf%nboxes
          r = r + volume(get_pbox(mf, i))
       end do
    else
       r = volume(get_boxarray(mf))
    end if
  end function imultifab_volume
  function lmultifab_volume(mf, all) result(r)
    integer(kind=ll_t) :: r
    type(lmultifab), intent(in) :: mf
    logical, optional :: all
    integer :: i
    if ( all ) then
       r = 0_ll_t
       do i = 1, mf%nboxes
          r = r + volume(get_pbox(mf, i))
       end do
    else
       r = volume(get_boxarray(mf))
    end if
  end function lmultifab_volume

  function multifab_nboxes(mf) result(r)
    type(multifab), intent(in) :: mf
    integer :: r
    r = mf%nboxes
  end function multifab_nboxes
  function imultifab_nboxes(mf) result(r)
    type(imultifab), intent(in) :: mf
    integer :: r
    r = mf%nboxes
  end function imultifab_nboxes

  function multifab_remote(mf, i) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i
    logical :: r
    r = layout_remote(mf%la, i)
  end function multifab_remote
  function imultifab_remote(mf, i) result(r)
    type(imultifab), intent(in) :: mf
    integer, intent(in) :: i
    logical :: r
    r = layout_remote(mf%la, i)
  end function imultifab_remote
  function lmultifab_remote(mf, i) result(r)
    type(lmultifab), intent(in) :: mf
    integer, intent(in) :: i
    logical :: r
    r = layout_remote(mf%la, i)
  end function lmultifab_remote

  function multifab_local(mf, i) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i
    logical :: r
    r = layout_local(mf%la, i)
  end function multifab_local
  function imultifab_local(mf, i) result(r)
    type(imultifab), intent(in) :: mf
    integer, intent(in) :: i
    logical :: r
    r = layout_local(mf%la, i)
  end function imultifab_local
  function lmultifab_local(mf, i) result(r)
    type(lmultifab), intent(in) :: mf
    integer, intent(in) :: i
    logical :: r
    r = layout_local(mf%la, i)
  end function lmultifab_local

  function multifab_get_boxarray(mf) result(r)
    type(boxarray) :: r
    type(multifab), intent(in) :: mf
    r = get_boxarray(mf%la)
  end function multifab_get_boxarray
  function imultifab_get_boxarray(mf) result(r)
    type(boxarray) :: r
    type(imultifab), intent(in) :: mf
    r = get_boxarray(mf%la)
  end function imultifab_get_boxarray
  function lmultifab_get_boxarray(mf) result(r)
    type(boxarray) :: r
    type(lmultifab), intent(in) :: mf
    r = get_boxarray(mf%la)
  end function lmultifab_get_boxarray

  function multifab_get_box(mf, i) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = get_box(mf%la, i)
  end function multifab_get_box
  function imultifab_get_box(mf, i) result(r)
    type(imultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = get_box(mf%la, i)
  end function imultifab_get_box
  function lmultifab_get_box(mf, i) result(r)
    type(lmultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = get_box(mf%la, i)
  end function lmultifab_get_box

  function multifab_get_ibox(mf, i) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = mf%fbs(i)%ibx
  end function multifab_get_ibox
  function imultifab_get_ibox(mf, i) result(r)
    type(imultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = mf%fbs(i)%ibx
  end function imultifab_get_ibox
  function lmultifab_get_ibox(mf, i) result(r)
    type(lmultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = mf%fbs(i)%ibx
  end function lmultifab_get_ibox

  function multifab_get_pbox(mf, i) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = mf%fbs(i)%pbx
  end function multifab_get_pbox
  function imultifab_get_pbox(mf, i) result(r)
    type(imultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = mf%fbs(i)%pbx
  end function imultifab_get_pbox
  function lmultifab_get_pbox(mf, i) result(r)
    type(lmultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = mf%fbs(i)%pbx
  end function lmultifab_get_pbox

  function multifab_dataptr(mf, i) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i
    real(dp_t), pointer :: r(:,:,:,:)
    r => mf%fbs(i)%p
  end function multifab_dataptr
  function imultifab_dataptr(mf, i) result(r)
    type(imultifab), intent(in) :: mf
    integer, intent(in) :: i
    integer, pointer :: r(:,:,:,:)
    r => mf%fbs(i)%p
  end function imultifab_dataptr
  function lmultifab_dataptr(mf, i) result(r)
    type(lmultifab), intent(in) :: mf
    integer, intent(in) :: i
    logical, pointer :: r(:,:,:,:)
    r => mf%fbs(i)%p
  end function lmultifab_dataptr

  function multifab_dataptr_c(mf, i, c, nc) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i, c
    integer, intent(in), optional :: nc
    real(dp_t), pointer :: r(:,:,:,:)
    r => fab_dataptr_c(mf%fbs(i), c, nc)
  end function multifab_dataptr_c
  function imultifab_dataptr_c(mf, i, c, nc) result(r)
    type(imultifab), intent(in) :: mf
    integer, intent(in) :: i, c
    integer, intent(in), optional :: nc
    integer, pointer :: r(:,:,:,:)
    r => ifab_dataptr_c(mf%fbs(i), c, nc)
  end function imultifab_dataptr_c
  function lmultifab_dataptr_c(mf, i, c, nc) result(r)
    type(lmultifab), intent(in) :: mf
    integer, intent(in) :: i, c
    integer, intent(in), optional :: nc
    logical, pointer :: r(:,:,:,:)
    r => lfab_dataptr_c(mf%fbs(i), c, nc)
  end function lmultifab_dataptr_c

  function multifab_dataptr_bx_c(mf, i, bx, c, nc) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i, c
    integer, intent(in), optional :: nc
    type(box), intent(in) :: bx
    real(dp_t), pointer :: r(:,:,:,:)
    r => fab_dataptr_bx_c(mf%fbs(i), bx, c, nc)
  end function multifab_dataptr_bx_c
  function imultifab_dataptr_bx_c(mf, i, bx, c, nc) result(r)
    type(imultifab), intent(in) :: mf
    integer, intent(in) :: i, c
    integer, intent(in), optional :: nc
    type(box), intent(in) :: bx
    integer, pointer :: r(:,:,:,:)
    r => ifab_dataptr_bx_c(mf%fbs(i), bx, c, nc)
  end function imultifab_dataptr_bx_c
  function lmultifab_dataptr_bx_c(mf, i, bx, c, nc) result(r)
    type(lmultifab), intent(in) :: mf
    integer, intent(in) :: i, c
    integer, intent(in), optional :: nc
    type(box), intent(in) :: bx
    logical, pointer :: r(:,:,:,:)
    r => lfab_dataptr_bx_c(mf%fbs(i), bx, c, nc)
  end function lmultifab_dataptr_bx_c

  function multifab_dataptr_bx(mf, i, bx) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box), intent(in) :: bx
    real(dp_t), pointer :: r(:,:,:,:)
    r => fab_dataptr_bx(mf%fbs(i), bx)
  end function multifab_dataptr_bx
  function imultifab_dataptr_bx(mf, i, bx) result(r)
    type(imultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box), intent(in) :: bx
    integer, pointer :: r(:,:,:,:)
    r => ifab_dataptr_bx(mf%fbs(i), bx)
  end function imultifab_dataptr_bx
  function lmultifab_dataptr_bx(mf, i, bx) result(r)
    type(lmultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box), intent(in) :: bx
    logical, pointer :: r(:,:,:,:)
    r => lfab_dataptr_bx(mf%fbs(i), bx)
  end function lmultifab_dataptr_bx

  subroutine multifab_setval(mf, val, all)
    type(multifab), intent(inout) :: mf
    real(kind=dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    type(box) :: bx
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,bx) SHARED(val)
    do i = 1, mf%nboxes
       if ( multifab_remote(mf, i) ) cycle
       if ( lall ) then
          call setval(mf%fbs(i), val)
       else
          bx = get_ibox(mf, i)
          call setval(mf%fbs(i), val, bx)
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_setval
  subroutine imultifab_setval(mf, val, all)
    type(imultifab), intent(inout) :: mf
    integer, intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    type(box) :: bx
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,bx) SHARED(val)
    do i = 1, mf%nboxes
       if ( imultifab_remote(mf, i) ) cycle
       if ( lall ) then
          call setval(mf%fbs(i), val)
       else
          bx = get_ibox(mf, i)
          call setval(mf%fbs(i), val, bx)
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine imultifab_setval
  subroutine lmultifab_setval(mf, val, all)
    type(lmultifab), intent(inout) :: mf
    logical, intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    type(box) :: bx
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,bx) SHARED(val)
    do i = 1, mf%nboxes
       if ( lmultifab_remote(mf, i) ) cycle
       if ( lall ) then
          call setval(mf%fbs(i), val)
       else
          bx = get_ibox(mf, i)
          call setval(mf%fbs(i), val, bx)
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine lmultifab_setval

  subroutine multifab_setval_bx(mf, val, bx, all)
    type(multifab), intent(inout) :: mf
    type(box), intent(in) :: bx
    real(kind=dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    type(box) :: bx1
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,bx1) SHARED(val)
    do i = 1, mf%nboxes
       if ( multifab_remote(mf, i) ) cycle
       if ( lall ) then
          bx1 = intersection(bx, get_pbox(mf, i))
          if ( .not. empty(bx1) ) &
               call setval(mf%fbs(i), val, bx1)
       else
          bx1 = intersection(bx, get_ibox(mf, i))
          if ( .not. empty(bx1) ) &
               call setval(mf%fbs(i), val, bx1)
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_setval_bx
  subroutine imultifab_setval_bx(mf, val, bx, all)
    type(imultifab), intent(inout) :: mf
    type(box), intent(in) :: bx
    integer, intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    type(box) :: bx1
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,bx1) SHARED(val)
    do i = 1, mf%nboxes
       if ( imultifab_remote(mf, i) ) cycle
       if ( lall ) then
          bx1 = intersection(bx, get_pbox(mf, i))
          if ( .not. empty(bx1) ) &
               call setval(mf%fbs(i), val, bx1)
       else
          bx1 = intersection(bx, get_ibox(mf, i))
          if ( .not. empty(bx1) ) &
               call setval(mf%fbs(i), val, bx1)
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine imultifab_setval_bx
  subroutine lmultifab_setval_bx(mf, val, bx, all)
    type(lmultifab), intent(inout) :: mf
    type(box), intent(in) :: bx
    logical, intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    type(box) :: bx1
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,bx1) SHARED(val)
    do i = 1, mf%nboxes
       if ( lmultifab_remote(mf, i) ) cycle
       if ( lall ) then
          bx1 = intersection(bx, get_pbox(mf, i))
          if ( .not. empty(bx1) ) &
               call setval(mf%fbs(i), val, bx1)
       else
          bx1 = intersection(bx, get_ibox(mf, i))
          if ( .not. empty(bx1) ) &
               call setval(mf%fbs(i), val, bx1)
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine lmultifab_setval_bx

  subroutine lmultifab_setval_ba(mf, val, ba, all)
    type(lmultifab), intent(inout) :: mf
    type(boxarray), intent(in) :: ba
    logical, intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i, n
    type(box) :: bx1
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,bx1) SHARED(val)
    do n = 1, ba%nboxes
       do i = 1, mf%nboxes
          if ( lmultifab_remote(mf, i) ) cycle
          if ( lall ) then
             bx1 = intersection(get_box(ba,n), get_pbox(mf, i))
             if ( .not. empty(bx1) ) &
                  call setval(mf%fbs(i), val, bx1)
          else
             bx1 = intersection(get_box(ba,n), get_ibox(mf, i))
             if ( .not. empty(bx1) ) &
                  call setval(mf%fbs(i), val, bx1)
          end if
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine lmultifab_setval_ba

  subroutine multifab_setval_bx_c(mf, val, bx, c, nc, all)
    type(multifab), intent(inout) :: mf
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    real(kind=dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    type(box) :: bx1
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,bx1) SHARED(val,c)
    do i = 1, mf%nboxes
       if ( multifab_remote(mf, i) ) cycle
       if ( lall ) then
          bx1 = intersection(bx, get_pbox(mf, i))
          if ( .not. empty(bx1) ) &
               call setval(mf%fbs(i), val, bx1, c, nc)
       else
          bx1 = intersection(bx, get_ibox(mf, i))
          if ( .not. empty(bx1) ) &
               call setval(mf%fbs(i), val, bx1, c, nc)
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_setval_bx_c
  subroutine imultifab_setval_bx_c(mf, val, bx, c, nc, all)
    type(imultifab), intent(inout) :: mf
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    integer, intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    type(box) :: bx1
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,bx1) SHARED(val,c)
    do i = 1, mf%nboxes
       if ( imultifab_remote(mf, i) ) cycle
       if ( lall ) then
          bx1 = intersection(bx, get_pbox(mf, i))
          if ( .not. empty(bx1) ) &
               call setval(mf%fbs(i), val, bx1, c, nc)
       else
          bx1 = intersection(bx, get_ibox(mf, i))
          if ( .not. empty(bx1) ) &
               call setval(mf%fbs(i), val, bx1, c, nc)
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine imultifab_setval_bx_c
  subroutine lmultifab_setval_bx_c(mf, val, bx, c, nc, all)
    type(lmultifab), intent(inout) :: mf
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    logical, intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    type(box) :: bx1
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,bx1) SHARED(val,c)
    do i = 1, mf%nboxes
       if ( lmultifab_remote(mf, i) ) cycle
       if ( lall ) then
          bx1 = intersection(bx, get_pbox(mf, i))
          if ( .not. empty(bx1) ) &
               call setval(mf%fbs(i), val, bx1, c, nc)
       else
          bx1 = intersection(bx, get_ibox(mf, i))
          if ( .not. empty(bx1) ) &
               call setval(mf%fbs(i), val, bx1, c, nc)
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine lmultifab_setval_bx_c

  subroutine multifab_setval_c(mf, val, c, nc, all)
    type(multifab), intent(inout) :: mf
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    real(kind=dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    type(box) :: bx
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,bx) SHARED(val,c)
    do i = 1, mf%nboxes
       if ( multifab_remote(mf, i) ) cycle
       if ( lall ) then
          call setval(mf%fbs(i), val, c, nc)
       else
          bx = get_ibox(mf, i)
          call setval(mf%fbs(i), val, bx, c, nc)
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_setval_c
  subroutine imultifab_setval_c(mf, val, c, nc, all)
    type(imultifab), intent(inout) :: mf
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    integer, intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    type(box) :: bx
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,bx) SHARED(val,c)
    do i = 1, mf%nboxes
       if ( imultifab_remote(mf, i) ) cycle
       if ( lall ) then
          call setval(mf%fbs(i), val, c, nc)
       else
          bx = get_ibox(mf, i)
          call setval(mf%fbs(i), val, bx, c, nc)
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine imultifab_setval_c
  subroutine lmultifab_setval_c(mf, val, c, nc, all)
    type(lmultifab), intent(inout) :: mf
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    logical, intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    type(box) :: bx
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,bx) SHARED(val,c)
    do i = 1, mf%nboxes
       if ( lmultifab_remote(mf, i) ) cycle
       if ( lall ) then
          call setval(mf%fbs(i), val, c, nc)
       else
          bx = get_ibox(mf, i)
          call setval(mf%fbs(i), val, bx, c, nc)
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine lmultifab_setval_c

  subroutine multifab_debug_fill(mf, all, loc)
    type(multifab), intent(inout) :: mf
    logical, intent(in), optional :: all
    logical, intent(in), optional :: loc
    integer :: i
    do i = 1, mf%nboxes
       if ( multifab_remote(mf, i) ) cycle
       call fab_debug_fill(mf%fbs(i), i, all, loc)
    end do
  end subroutine multifab_debug_fill
  subroutine imultifab_debug_fill(mf, all)
    type(imultifab), intent(inout) :: mf
    logical, intent(in), optional :: all
    integer :: i
    type(box) :: bx
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    do i = 1, mf%nboxes
       if ( imultifab_remote(mf, i) ) cycle
       if ( lall ) then
          bx = get_pbox(mf, i)
       else
          call setval(mf%fbs(i), -i)
          bx = get_ibox(mf, i)
       end if
       call setval(mf%fbs(i), i, bx)
    end do
  end subroutine imultifab_debug_fill

  !! MULTIFAB_SET_BOUNDARY_VAL
  !! BUG: sets ghost cell values of the fabs, but
  !! will corrupt f/f ghost cell values
  subroutine multifab_set_border_val(mf, val)
    type(multifab), intent(inout) :: mf
    real(kind=dp_t), intent(in), optional :: val
    integer :: i
    do i = 1, mf%nboxes; if ( multifab_remote(mf, i) ) cycle
       call fab_set_border_val(mf%fbs(i), val)
    end do
  end subroutine multifab_set_border_val

  subroutine multifab_fill_boundary(mf, ng, nocomm, cross)
    type(multifab), intent(inout) :: mf
    integer, intent(in), optional :: ng
    logical, intent(in), optional :: nocomm, cross

    integer :: lng
    logical :: lnocomm, lcross

    lcross  = .false.; if ( present(cross) )  lcross  = cross
    lnocomm = .false.; if ( present(nocomm) ) lnocomm = nocomm
    lng     = mf%ng;   if ( present(ng) )     lng     = ng

    if ( lng > mf%ng ) call bl_error("MULTIFAB_FILL_BOUNDARY: ng too large", ng)

    if ( lng < 1 ) return

    if ( d_fb_fancy ) then
       call fancy()
    else
       call easy()
    end if

  contains

    subroutine easy()
      real(dp_t), dimension(:,:,:,:), pointer :: pdst, psrc
      type(boxarray)                          :: bxai
      type(box)                               :: abx
      integer                                 :: i, j, ii, proc
      integer                                 :: shft(2*3**mf%dim,mf%dim)
      integer, parameter                      :: tag = 1101

      do i = 1, mf%nboxes
         call boxarray_bndry_periodic(bxai, mf%la%lap%pd, get_box(mf,i), mf%nodal, mf%la%lap%pmask, lng, shft)
         do j = 1, mf%nboxes
            if ( remote(mf,i) .and. remote(mf,j) ) cycle
            do ii = 1, bxai%nboxes
               abx = intersection(get_ibox(mf,j), bxai%bxs(ii))
               if ( .not. empty(abx) ) then
                  if ( local(mf,i) .and. local(mf,j) ) then
                     psrc => dataptr(mf, j, abx)
                     pdst => dataptr(mf, i, shift(abx,-shft(ii,:)))
                     pdst = psrc
                  else if ( .not. lnocomm ) then
                     if ( local(mf,j) ) then ! must send
                        psrc => dataptr(mf, j, abx)
                        proc = get_proc(mf%la, i)
                        call parallel_send(psrc, proc, tag)
                     else if ( local(mf,i) ) then  ! must recv
                        pdst => dataptr(mf, i, shift(abx,-shft(ii,:)))
                        proc = get_proc(mf%la,j)
                        call parallel_recv(pdst, proc, tag)
                     end if
                  end if
               end if
            end do
         end do
         call destroy(bxai)
      end do
    end subroutine easy

    subroutine fancy()
      real(kind=dp_t), dimension(:,:,:,:), pointer :: p1, p2
      real(kind=dp_t), dimension(:,:,:,:), pointer :: p
      integer, allocatable :: rst(:)
      integer, allocatable :: sst(:)
      integer :: i, ii, jj
      type(box) ::sbx, dbx
      integer, parameter :: tag = 1102
      integer :: nc, sh(MAX_SPACEDIM+1)
      type(boxassoc) :: bxasc

      bxasc = layout_boxassoc(mf%la, lng, mf%nodal, lcross)
!call boxassoc_print(bxasc, unit = 1)

      !$OMP PARALLEL DO PRIVATE(i,ii,jj,sbx,dbx,p1,p2)
      do i = 1, bxasc%l_con%ncpy
         ii = bxasc%l_con%cpy(i)%nd
         jj = bxasc%l_con%cpy(i)%ns
         sbx = bxasc%l_con%cpy(i)%sbx
         dbx = bxasc%l_con%cpy(i)%dbx
         p1 => dataptr(mf%fbs(ii), dbx)
         p2 => dataptr(mf%fbs(jj), sbx)
!print *, i
!print *, ii, dbx
!print *, jj, sbx
         p1 = p2
      end do
      !$OMP END PARALLEL DO

      if ( lnocomm ) return
      
      if ( bxasc%r_con%svol > 0 ) then
         if ( allocated(g_snd_d) ) then
            if ( size(g_snd_d) < bxasc%r_con%svol ) then
               deallocate(g_snd_d)
               allocate(g_snd_d(bxasc%r_con%svol))
            end if
         else
            allocate(g_snd_d(bxasc%r_con%svol))
         end if
      end if
      if ( bxasc%r_con%rvol > 0 ) then
         if ( allocated(g_rcv_d) ) then
            if ( size(g_rcv_d) < bxasc%r_con%rvol ) then
               deallocate(g_rcv_d)
               allocate(g_rcv_d(bxasc%r_con%rvol))
            end if
         else
            allocate(g_rcv_d(bxasc%r_con%rvol))
         end if
      end if

      ! the boxassoc contains size data in terms of numpts, must use
      ! nc to get the actual volume

      nc = mf%nc
      do i = 1, bxasc%r_con%nsnd
         p => dataptr(mf, bxasc%r_con%snd(i)%ns, bxasc%r_con%snd(i)%sbx)
         g_snd_d(1 + nc*bxasc%r_con%snd(i)%pv:nc*bxasc%r_con%snd(i)%av) = &
              reshape(p, nc*bxasc%r_con%snd(i)%s1)
      end do
      allocate(rst(bxasc%r_con%nrp), sst(bxasc%r_con%nsp))
      !
      ! Always do recv's asynchronously.
      !
      do i = 1, bxasc%r_con%nrp
         rst(i) = parallel_irecv_dv(g_rcv_d(1+nc*bxasc%r_con%rtr(i)%pv:), &
              nc*bxasc%r_con%rtr(i)%sz, bxasc%r_con%rtr(i)%pr, tag)
      end do
      if ( d_fb_async ) then
         do i = 1, bxasc%r_con%nsp
            sst(i) = parallel_isend_dv(g_snd_d(1+nc*bxasc%r_con%str(i)%pv:), &
                 nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
         end do
      else
         do i = 1, bxasc%r_con%nsp
            call parallel_send_dv(g_snd_d(1+nc*bxasc%r_con%str(i)%pv), &
                 nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
         end do
      end if
      call parallel_wait(rst)
      do i = 1, bxasc%r_con%nrcv
         sh = bxasc%r_con%rcv(i)%sh
         sh(4) = nc
         p => dataptr(mf, bxasc%r_con%rcv(i)%nd, bxasc%r_con%rcv(i)%dbx)
         p =  reshape( &
              g_rcv_d(1 + nc*bxasc%r_con%rcv(i)%pv:nc*bxasc%r_con%rcv(i)%av), &
              sh)
      end do
      if ( d_fb_async) call parallel_wait(sst)
    end subroutine fancy
  end subroutine multifab_fill_boundary

  subroutine multifab_fill_boundary_c(mf, c, nc, ng, nocomm, cross)
    type(multifab), intent(inout) :: mf
    integer, intent(in)           :: c, nc
    integer, intent(in), optional :: ng
    logical, intent(in), optional :: nocomm, cross

    integer :: lng
    logical :: lnocomm, lcross

    lcross  = .false.; if ( present(cross) )  lcross  = cross
    lnocomm = .false.; if ( present(nocomm) ) lnocomm = nocomm
    lng     = mf%ng;   if ( present(ng) )     lng     = ng

    if ( lng > mf%ng ) call bl_error("MULTIFAB_FILL_BOUNDARY: ng too large", ng)

    if ( lng < 1 ) return

    if ( d_fb_fancy ) then
       call fancy()
    else
       call easy()
    end if

  contains

    subroutine easy()
      real(dp_t), dimension(:,:,:,:), pointer :: pdst, psrc
      type(boxarray)                          :: bxai
      type(box)                               :: abx
      integer                                 :: i, j, ii, proc
      integer                                 :: shft(2*3**mf%dim,mf%dim)
      integer, parameter                      :: tag = 1101

      do i = 1, mf%nboxes
         call boxarray_bndry_periodic(bxai, mf%la%lap%pd, get_box(mf,i), mf%nodal, mf%la%lap%pmask, lng, shft)
         do j = 1, mf%nboxes
            if ( remote(mf,i) .and. remote(mf,j) ) cycle
            do ii = 1, bxai%nboxes
               abx = intersection(get_ibox(mf,j), bxai%bxs(ii))
               if ( .not. empty(abx) ) then
                  if ( local(mf,i) .and. local(mf,j) ) then
                     psrc => dataptr(mf, j, abx, c, nc)
                     pdst => dataptr(mf, i, shift(abx,-shft(ii,:)), c, nc)
                     pdst = psrc
                  else if ( .not. lnocomm ) then
                     if ( local(mf,j) ) then ! must send
                        psrc => dataptr(mf, j, abx, c, nc)
                        proc = get_proc(mf%la, i)
                        call parallel_send(psrc, proc, tag)
                     else if ( local(mf,i) ) then  ! must recv
                        pdst => dataptr(mf, i, shift(abx,-shft(ii,:)), c, nc)
                        proc = get_proc(mf%la,j)
                        call parallel_recv(pdst, proc, tag)
                     end if
                  end if
               end if
            end do
         end do
         call destroy(bxai)
      end do
    end subroutine easy

    subroutine fancy()
      real(kind=dp_t), dimension(:,:,:,:), pointer :: p1, p2
      real(kind=dp_t), dimension(:,:,:,:), pointer :: p
      integer, allocatable :: rst(:)
      integer, allocatable :: sst(:)
      integer :: i, ii, jj
      type(box) ::sbx, dbx
      integer, parameter :: tag = 1102
      integer :: sh(MAX_SPACEDIM+1)
      type(boxassoc) :: bxasc

      bxasc = layout_boxassoc(mf%la, lng, mf%nodal, lcross)
!call boxassoc_print(bxasc, unit = 1)

      !$OMP PARALLEL DO PRIVATE(i,ii,jj,sbx,dbx,p1,p2)
      do i = 1, bxasc%l_con%ncpy
         ii = bxasc%l_con%cpy(i)%nd
         jj = bxasc%l_con%cpy(i)%ns
         sbx = bxasc%l_con%cpy(i)%sbx
         dbx = bxasc%l_con%cpy(i)%dbx
         p1 => dataptr(mf%fbs(ii), dbx, c, nc)
         p2 => dataptr(mf%fbs(jj), sbx, c, nc)
!print *, i
!print *, ii, dbx
!print *, jj, sbx
         p1 = p2
      end do
      !$OMP END PARALLEL DO

      if ( lnocomm ) return
      
      if ( bxasc%r_con%svol > 0 ) then
         if ( allocated(g_snd_d) ) then
            if ( size(g_snd_d) < bxasc%r_con%svol ) then
               deallocate(g_snd_d)
               allocate(g_snd_d(bxasc%r_con%svol))
            end if
         else
            allocate(g_snd_d(bxasc%r_con%svol))
         end if
      end if
      if ( bxasc%r_con%rvol > 0 ) then
         if ( allocated(g_rcv_d) ) then
            if ( size(g_rcv_d) < bxasc%r_con%rvol ) then
               deallocate(g_rcv_d)
               allocate(g_rcv_d(bxasc%r_con%rvol))
            end if
         else
            allocate(g_rcv_d(bxasc%r_con%rvol))
         end if
      end if

      ! the boxassoc contains size data in terms of numpts, must use
      ! nc to get the actual volume

      do i = 1, bxasc%r_con%nsnd
         p => dataptr(mf, bxasc%r_con%snd(i)%ns, bxasc%r_con%snd(i)%sbx, c, nc)
         g_snd_d(1 + nc*bxasc%r_con%snd(i)%pv:nc*bxasc%r_con%snd(i)%av) = &
              reshape(p, nc*bxasc%r_con%snd(i)%s1)
      end do
      allocate(rst(bxasc%r_con%nrp), sst(bxasc%r_con%nsp))
      !
      ! Always do recv's asynchronously.
      !
      do i = 1, bxasc%r_con%nrp
         rst(i) = parallel_irecv_dv(g_rcv_d(1+nc*bxasc%r_con%rtr(i)%pv:), &
              nc*bxasc%r_con%rtr(i)%sz, bxasc%r_con%rtr(i)%pr, tag)
      end do
      if ( d_fb_async ) then
         do i = 1, bxasc%r_con%nsp
            sst(i) = parallel_isend_dv(g_snd_d(1+nc*bxasc%r_con%str(i)%pv:), &
                 nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
         end do
      else
         do i = 1, bxasc%r_con%nsp
            call parallel_send_dv(g_snd_d(1+nc*bxasc%r_con%str(i)%pv), &
                 nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
         end do
      end if
      call parallel_wait(rst)
      do i = 1, bxasc%r_con%nrcv
         sh = bxasc%r_con%rcv(i)%sh
         sh(4) = nc
         p => dataptr(mf, bxasc%r_con%rcv(i)%nd, bxasc%r_con%rcv(i)%dbx, c, nc)
         p =  reshape( &
              g_rcv_d(1 + nc*bxasc%r_con%rcv(i)%pv:nc*bxasc%r_con%rcv(i)%av), &
              sh)
      end do
      if ( d_fb_async) call parallel_wait(sst)
    end subroutine fancy
  end subroutine multifab_fill_boundary_c

  subroutine imultifab_fill_boundary(mf, ng, nocomm, cross)
    type(imultifab), intent(inout) :: mf
    integer, intent(in), optional  :: ng
    logical, intent(in), optional  :: nocomm, cross

    integer :: lng
    logical :: lnocomm, lcross

    lcross  = .false.; if ( present(cross) )  lcross  = cross
    lnocomm = .false.; if ( present(nocomm) ) lnocomm = nocomm
    lng     = mf%ng;   if ( present(ng) )     lng     = ng

    if ( lng > mf%ng ) &
         call bl_error("IMULTIFAB_FILL_BOUNDARY: ng too large", ng)
    if ( lng < 1 ) return

    if ( i_fb_fancy ) then
       call fancy()
    else
       call easy()
    end if

  contains

    subroutine easy()
      integer, dimension(:,:,:,:), pointer :: pdst, psrc
      type(boxarray)                       :: bxai
      type(box)                            :: abx
      integer                              :: i, j, ii, proc
      integer                              :: shft(2*3**mf%dim,mf%dim)
      integer, parameter                   :: tag = 1101

      do i = 1, mf%nboxes
         call boxarray_bndry_periodic(bxai, mf%la%lap%pd, get_box(mf,i), mf%nodal, mf%la%lap%pmask, lng, shft)
         do j = 1, mf%nboxes
            if ( remote(mf,i) .and. remote(mf,j) ) cycle
            do ii = 1, bxai%nboxes
               abx = intersection(get_ibox(mf,j), bxai%bxs(ii))
               if ( .not. empty(abx) ) then
                  if ( local(mf,i) .and. local(mf,j) ) then
                     psrc => dataptr(mf, j, abx)
                     pdst => dataptr(mf, i, shift(abx,-shft(ii,:)))
                     pdst = psrc
                  else if ( .not. lnocomm ) then
                     if ( local(mf,j) ) then ! must send
                        psrc => dataptr(mf, j, abx)
                        proc = get_proc(mf%la, i)
                        call parallel_send(psrc, proc, tag)
                     else if ( local(mf,i) ) then  ! must recv
                        pdst => dataptr(mf, i, shift(abx,-shft(ii,:)))
                        proc = get_proc(mf%la,j)
                        call parallel_recv(pdst, proc, tag)
                     end if
                  end if
               end if
            end do
         end do
         call destroy(bxai)
      end do
    end subroutine easy

    subroutine fancy()
      integer, dimension(:,:,:,:), pointer :: p1, p2
      integer, dimension(:,:,:,:), pointer :: p
      integer, allocatable :: rst(:)
      integer, allocatable :: sst(:)
      integer :: i, ii, jj
      type(box) :: sbx, dbx
      integer, parameter :: tag = 1102
      integer :: nc, sh(MAX_SPACEDIM+1)
      type(boxassoc) :: bxasc

      bxasc = layout_boxassoc(mf%la, lng, mf%nodal, lcross)

      do i = 1, bxasc%l_con%ncpy
         ii = bxasc%l_con%cpy(i)%nd
         jj = bxasc%l_con%cpy(i)%ns
         sbx = bxasc%l_con%cpy(i)%sbx
         dbx = bxasc%l_con%cpy(i)%dbx
         p1 => dataptr(mf%fbs(ii), dbx)
         p2 => dataptr(mf%fbs(jj), sbx)
         p1 = p2
      end do

      if ( lnocomm ) return
      
      if ( bxasc%r_con%svol > 0 ) then
         if ( allocated(g_snd_i) ) then
            if ( size(g_snd_i) < bxasc%r_con%svol ) then
               deallocate(g_snd_i)
               allocate(g_snd_i(bxasc%r_con%svol))
            end if
         else
            allocate(g_snd_i(bxasc%r_con%svol))
         end if
      end if
      if ( bxasc%r_con%rvol > 0 ) then
         if ( allocated(g_rcv_i) ) then
            if ( size(g_rcv_i) < bxasc%r_con%rvol ) then
               deallocate(g_rcv_i)
               allocate(g_rcv_i(bxasc%r_con%rvol))
            end if
         else
            allocate(g_rcv_i(bxasc%r_con%rvol))
         end if
      end if

      ! the boxassoc contains size data in terms of numpts, must use
      ! nc to get the actual volume

      nc = mf%nc
      do i = 1, bxasc%r_con%nsnd
         p => dataptr(mf, bxasc%r_con%snd(i)%ns, bxasc%r_con%snd(i)%sbx)
         g_snd_i(1 + nc*bxasc%r_con%snd(i)%pv:nc*bxasc%r_con%snd(i)%av) = &
              reshape(p, nc*bxasc%r_con%snd(i)%s1)
      end do
      allocate(rst(bxasc%r_con%nrp), sst(bxasc%r_con%nsp))
      !
      ! Always do recv's asynchronously.
      !
      do i = 1, bxasc%r_con%nrp
         rst(i) = parallel_irecv_iv(g_rcv_i(1+nc*bxasc%r_con%rtr(i)%pv:), &
              nc*bxasc%r_con%rtr(i)%sz, bxasc%r_con%rtr(i)%pr, tag)
      end do
      if ( i_fb_async ) then
         do i = 1, bxasc%r_con%nsp
            sst(i) = parallel_isend_iv(g_snd_i(1+nc*bxasc%r_con%str(i)%pv:), &
                 nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
         end do
      else
         do i = 1, bxasc%r_con%nsp
            call parallel_send_iv(g_snd_i(1+nc*bxasc%r_con%str(i)%pv:), &
                 nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
         end do
      end if
      call parallel_wait(rst)
      do i = 1, bxasc%r_con%nrcv
         sh = bxasc%r_con%rcv(i)%sh
         sh(4) = nc
         p => dataptr(mf, bxasc%r_con%rcv(i)%nd, bxasc%r_con%rcv(i)%dbx)
         p =  reshape( &
              g_rcv_i(1 + nc*bxasc%r_con%rcv(i)%pv:nc*bxasc%r_con%rcv(i)%av), &
              sh)
      end do
      if ( i_fb_async) call parallel_wait(sst)
    end subroutine fancy
  end subroutine imultifab_fill_boundary

  subroutine lmultifab_fill_boundary(mf, ng, nocomm, cross)
    type(lmultifab), intent(inout) :: mf
    integer, intent(in), optional  :: ng
    logical, intent(in), optional  :: nocomm, cross

    integer :: lng
    logical :: lnocomm, lcross

    lcross  = .false.; if ( present(cross) )  lcross  = cross
    lnocomm = .false.; if ( present(nocomm) ) lnocomm = nocomm
    lng     = mf%ng;   if ( present(ng) )     lng     = ng

    if ( lng > mf%ng ) &
         call bl_error("LMULTIFAB_FILL_BOUNDARY: ng too large", ng)
    if ( lng < 1 ) return

    if ( l_fb_fancy ) then
       call fancy()
    else
       call easy()
    end if

  contains

    subroutine easy()
      logical, dimension(:,:,:,:), pointer :: pdst, psrc
      type(boxarray)                       :: bxai
      type(box)                            :: abx
      integer                              :: i, j, ii, proc
      integer                              :: shft(2*3**mf%dim,mf%dim)
      integer, parameter                   :: tag = 1101

      do i = 1, mf%nboxes
         call boxarray_bndry_periodic(bxai, mf%la%lap%pd, get_box(mf,i), mf%nodal, mf%la%lap%pmask, lng, shft)
         do j = 1, mf%nboxes
            if ( remote(mf,i) .and. remote(mf,j) ) cycle
            do ii = 1, bxai%nboxes
               abx = intersection(get_ibox(mf,j), bxai%bxs(ii))
               if ( .not. empty(abx) ) then
                  if ( local(mf,i) .and. local(mf,j) ) then
                     psrc => dataptr(mf, j, abx)
                     pdst => dataptr(mf, i, shift(abx,-shft(ii,:)))
                     pdst = psrc
                  else if ( .not. lnocomm ) then
                     if ( local(mf,j) ) then ! must send
                        psrc => dataptr(mf, j, abx)
                        proc = get_proc(mf%la, i)
                        call parallel_send(psrc, proc, tag)
                     else if ( local(mf,i) ) then  ! must recv
                        pdst => dataptr(mf, i, shift(abx,-shft(ii,:)))
                        proc = get_proc(mf%la,j)
                        call parallel_recv(pdst, proc, tag)
                     end if
                  end if
               end if
            end do
         end do
         call destroy(bxai)
      end do
    end subroutine easy

    subroutine fancy()
      logical, dimension(:,:,:,:), pointer :: p1, p2
      logical, dimension(:,:,:,:), pointer :: p
      integer, allocatable :: rst(:)
      integer, allocatable :: sst(:)
      integer :: i, ii, jj
      type(box) :: sbx, dbx
      integer, parameter :: tag = 1102
      integer :: nc, sh(MAX_SPACEDIM+1)
      type(boxassoc) :: bxasc

      bxasc = layout_boxassoc(mf%la, lng, mf%nodal, lcross)

      do i = 1, bxasc%l_con%ncpy
         ii = bxasc%l_con%cpy(i)%nd
         jj = bxasc%l_con%cpy(i)%ns
         sbx = bxasc%l_con%cpy(i)%sbx
         dbx = bxasc%l_con%cpy(i)%dbx
         p1 => dataptr(mf%fbs(ii), dbx)
         p2 => dataptr(mf%fbs(jj), sbx)
         p1 = p2
      end do

      if ( lnocomm ) return
      
      if ( bxasc%r_con%svol > 0 ) then
         if ( allocated(g_snd_l) ) then
            if ( size(g_snd_l) < bxasc%r_con%svol ) then
               deallocate(g_snd_l)
               allocate(g_snd_l(bxasc%r_con%svol))
            end if
         else
            allocate(g_snd_l(bxasc%r_con%svol))
         end if
      end if
      if ( bxasc%r_con%rvol > 0 ) then
         if ( allocated(g_rcv_l) ) then
            if ( size(g_rcv_l) < bxasc%r_con%rvol ) then
               deallocate(g_rcv_l)
               allocate(g_rcv_l(bxasc%r_con%rvol))
            end if
         else
            allocate(g_rcv_l(bxasc%r_con%rvol))
         end if
      end if

      ! the boxassoc contains size data in terms of numpts, must use
      ! nc to get the actual volume

      nc = mf%nc
      do i = 1, bxasc%r_con%nsnd
         p => dataptr(mf, bxasc%r_con%snd(i)%ns, bxasc%r_con%snd(i)%sbx)
         g_snd_l(1 + nc*bxasc%r_con%snd(i)%pv:nc*bxasc%r_con%snd(i)%av) = &
              reshape(p, nc*bxasc%r_con%snd(i)%s1)
      end do
      allocate(rst(bxasc%r_con%nrp), sst(bxasc%r_con%nsp))
      !
      ! Always do recv's asynchronously.
      !
      do i = 1, bxasc%r_con%nrp
         rst(i) = parallel_irecv_lv(g_rcv_l(1+nc*bxasc%r_con%rtr(i)%pv:), &
              nc*bxasc%r_con%rtr(i)%sz, bxasc%r_con%rtr(i)%pr, tag)
      end do
      if ( l_fb_async ) then
         do i = 1, bxasc%r_con%nsp
            sst(i) = parallel_isend_lv(g_snd_l(1+nc*bxasc%r_con%str(i)%pv:), &
                 nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
         end do
      else
         do i = 1, bxasc%r_con%nsp
            call parallel_send_lv(g_snd_l(1+nc*bxasc%r_con%str(i)%pv:), &
                 nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
         end do
      end if
      call parallel_wait(rst)
      do i = 1, bxasc%r_con%nrcv
         sh = bxasc%r_con%rcv(i)%sh
         sh(4) = nc
         p => dataptr(mf, bxasc%r_con%rcv(i)%nd, bxasc%r_con%rcv(i)%dbx)
         p =  reshape( &
              g_rcv_l(1 + nc*bxasc%r_con%rcv(i)%pv:nc*bxasc%r_con%rcv(i)%av), &
              sh)
      end do
      if ( l_fb_async) call parallel_wait(sst)
    end subroutine fancy

  end subroutine lmultifab_fill_boundary

  subroutine multifab_internal_sync_shift(dmn,bx,pmask,nodal,shft,cnt)
    type(box), intent(in)  :: dmn,bx
    integer,   intent(out) :: shft(:,:),cnt
    logical,   intent(in)  :: pmask(:),nodal(:)

    type(box) :: dom,src
    integer   :: nbeg(3),nend(3),ldom(3),r(3),ri,rj,rk,l(3)
    !
    ! First a zero shift to represent the original box.
    !
    cnt = 1
    shft(cnt,:) = 0

    if (all(pmask .eqv. .false.)) return

    dom = box_nodalize(dmn,nodal)

    if (any(nodal .eqv. .true.)) then
       if (box_contains_strict(dom,bx)) return
    else
       if (box_contains(dom,bx)) return
    end if

    l(:) = 0; where(nodal) l = 1
    
    nbeg           = 0
    nend           = 0
    nbeg(1:bx%dim) = -1
    nend(1:bx%dim) = +1
    src            = bx
    ldom           = (/ extent(dom,1)-l(1),extent(dom,2)-l(2),extent(dom,3)-l(2) /)

    do ri = nbeg(1), nend(1)
       if (ri /= 0 .and. (.not. is_periodic(1))) cycle
       if (ri /= 0 .and. is_periodic(1)) src = shift(src,ri*ldom(1),1)

       do rj = nbeg(2), nend(2)
          if (rj /= 0 .and. (.not. is_periodic(2))) cycle
          if (rj /= 0 .and. is_periodic(2)) src = shift(src,rj*ldom(2),2)

          do rk = nbeg(3), nend(3)
             if (rk /= 0 .and. (.not. is_periodic(3))) cycle
             if (rk /= 0 .and. is_periodic(3)) src = shift(src,rk*ldom(3),3)

             if (ri == 0 .and. rj == 0 .and. rk == 0) cycle

             if (intersects(dom,src)) then
                cnt = cnt + 1
                r = (/ri*ldom(1),rj*ldom(2),rk*ldom(3)/)
                shft(cnt,1:bx%dim) = r(1:bx%dim)
             end if

             if (rk /= 0 .and. is_periodic(3)) src = shift(src,-rk*ldom(3),3)
          end do

          if (rj /= 0 .and. is_periodic(2)) src = shift(src,-rj*ldom(2),2)
       end do

       if (ri /= 0 .and. is_periodic(1)) src = shift(src,-ri*ldom(1),1)
    end do

    contains

      function is_periodic(i) result(r)
        integer, intent(in) :: i
        logical             :: r
        r = .false.
        if (i >= 1 .and. i <= bx%dim) r = pmask(i) .eqv. .true.
      end function is_periodic

  end subroutine multifab_internal_sync_shift
  !!
  !! Internal Sync makes sure that any overlapped values are reconciled
  !! by copying values from the lower index number fabs to the higher index
  !! numbered boxes.  Works cell centered and node centered.  Though in a typical
  !! cell-centered multifab, there are no overlaps to reconcile.
  !! if ALL is true then even ghost cell data is 'reconciled'
  !!
  subroutine multifab_internal_sync(mf, all, filter)
    type(multifab), intent(inout)               :: mf
    logical, intent(in), optional               :: all
    type(box)                                   :: ibx, jbx, abx, sbx
    real(dp_t), dimension(:,:,:,:), pointer     :: pdst, psrc
    real(dp_t), dimension(:,:,:,:), allocatable :: pt
    integer                                     :: i, j, jj, proc, cnt
    integer                                     :: shft(3**mf%dim,mf%dim)
    integer, parameter                          :: tag = 1104
    logical                                     :: lall

    interface
       subroutine filter(out, in)
         use bl_types
         real(dp_t), intent(inout) :: out(:,:,:,:)
         real(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface

    optional filter

    lall = .false. ; if ( present(all) ) lall = all

    do j = 1, mf%nboxes
       if ( lall ) then
          jbx = get_pbox(mf,j)
       else
          jbx = get_ibox(mf,j)
       end if
       call multifab_internal_sync_shift(mf%la%lap%pd, jbx, mf%la%lap%pmask, mf%nodal, shft, cnt)
       do jj = 1, cnt
          do i = j, mf%nboxes
             !
             ! Do not overwrite ourselves.
             !
             if (i == j .and. .not. any(shft(jj,:) /= 0)) cycle
             if ( lall ) then
                ibx = get_pbox(mf,i)
             else
                ibx = get_ibox(mf,i)
             end if
             abx = intersection(ibx, shift(jbx,shft(jj,:)))
             if ( .not. empty(abx) ) then
                if ( local(mf, i) .and. local(mf, j) ) then
                   pdst => dataptr(mf, i, abx)
                   psrc => dataptr(mf, j, shift(abx,-shft(jj,:)))
                   if ( present(filter) ) then
                      call filter(pdst, psrc)
                   else
                      pdst = psrc
                   end if
                else if ( local(mf, j) ) then ! must send
                   proc = get_proc(mf%la, i)
                   psrc => dataptr(mf, j, shift(abx,-shft(jj,:)))
                   call parallel_send(psrc, proc, tag)
                else if ( local(mf, i) ) then  ! must recv
                   proc = get_proc(mf%la,j)
                   pdst => dataptr(mf, i, abx)
                   if ( present(filter) ) then
                      allocate(pt(size(pdst,1),size(pdst,2),size(pdst,3),size(pdst,4)))
                      call parallel_recv(pt, proc, tag)
                      call filter(pdst, pt)
                      deallocate(pt)
                   else
                      call parallel_recv(pdst, proc, tag)
                   end if
                end if
             end if
          end do
       end do
    end do
  end subroutine multifab_internal_sync
  !!
  !! Internal Sync makes sure that any overlapped values are reconciled
  !! by coping values from the lower index number fabs to the higher index
  !! numbered boxes.  Works cell centered and node centered.  Though in a typical
  !! cell-centered multifab, there are no overlaps to reconcile.
  !!
  subroutine lmultifab_internal_sync(mf, all, filter)
    type(lmultifab), intent(inout)           :: mf
    logical, intent(in), optional            :: all
    type(box)                                :: ibx, jbx, abx, sbx
    logical, dimension(:,:,:,:), pointer     :: pdst, psrc
    logical, dimension(:,:,:,:), allocatable :: pt
    integer                                  :: i, j, jj, proc, cnt
    integer                                  :: shft(3**mf%dim,mf%dim)
    integer, parameter                       :: tag = 1104
    logical                                  :: lall

    interface
       subroutine filter(out, in)
         logical, intent(inout) :: out(:,:,:,:)
         logical, intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface

    optional filter

    lall = .false. ; if ( present(all) ) lall = all

    do j = 1, mf%nboxes
       if ( lall ) then
          jbx = get_pbox(mf,j)
       else
          jbx = get_ibox(mf,j)
       end if
       call multifab_internal_sync_shift(mf%la%lap%pd, jbx, mf%la%lap%pmask, mf%nodal, shft, cnt)
       do jj = 1, cnt
          do i = j, mf%nboxes
             !
             ! Do not overwrite ourselves.
             !
             if (i == j .and. .not. any(shft(jj,:) /= 0)) cycle
             if ( lall ) then
                ibx = get_pbox(mf,i)
             else
                ibx = get_ibox(mf,i)
             end if
             abx = intersection(ibx, shift(jbx,shft(jj,:)))
             if ( .not. empty(abx) ) then
                if ( local(mf, i) .and. local(mf, j) ) then
                   pdst => dataptr(mf, i, abx)
                   psrc => dataptr(mf, j, shift(abx,-shft(jj,:)))
                   if ( present(filter) ) then
                      call filter(pdst, psrc)
                   else
                      pdst = psrc
                   end if
                else if ( local(mf, j) ) then ! must send
                   proc = get_proc(mf%la, i)
                   psrc => dataptr(mf, j, shift(abx,-shft(jj,:)))
                   call parallel_send(psrc, proc, tag)
                else if ( local(mf, i) ) then  ! must recv
                   proc = get_proc(mf%la,j)
                   pdst => dataptr(mf, i, abx)
                   if ( present(filter) ) then
                      allocate(pt(size(pdst,1),size(pdst,2),size(pdst,3),size(pdst,4)))
                      call parallel_recv(pt, proc, tag)
                      call filter(pdst, pt)
                      deallocate(pt)
                   else
                      call parallel_recv(pdst, proc, tag)
                   end if
                end if
             end if
          end do
       end do
    end do
  end subroutine lmultifab_internal_sync

  subroutine multifab_print(mf, str, unit, all, data, skip)
    use bl_IO_module
    type(multifab), intent(in) :: mf
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: all, data
    integer, intent(in), optional :: skip
    integer :: i, ii
    integer :: un
    character(len=5) :: fn
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un, fmt='("MULTIFAB")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(": ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" DIM     = ",i2)') mf%dim
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NC      = ",i2)') mf%nc
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NG      = ",i2)') mf%ng
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NODAL   = ",3(L2,1X))') mf%nodal
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NBOXES  = ",i2)') mf%nboxes
    do ii = 0, parallel_nprocs()
       if ( ii == parallel_myproc() ) then
          do i = 1, mf%nboxes; if ( remote(mf,i) ) cycle
             write(unit=fn, fmt='(i5)') i
             call print(mf%fbs(i), str = fn, unit = unit, all = all, data = data, &
                  skip = unit_get_skip(skip) + 2)
          end do
       end if
       call parallel_barrier()
    end do
  end subroutine multifab_print
  subroutine imultifab_print(mf, str, unit, all, data, skip)
    use bl_IO_module
    type(imultifab), intent(in) :: mf
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: all, data
    integer, intent(in), optional :: skip
    integer :: i, ii
    integer :: un
    character(len=5) :: fn
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un, fmt='("IMULTIFAB")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(": ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" DIM     = ",i2)') mf%dim
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NC      = ",i2)') mf%nc
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NG      = ",i2)') mf%ng
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NODAL   = ",3(L2,1X))') mf%nodal
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NBOXES  = ",i2)') mf%nboxes
    do ii = 0, parallel_nprocs()
       if ( ii == parallel_myproc() ) then
          do i = 1, mf%nboxes; if ( remote(mf,i) ) cycle
             write(unit=fn, fmt='(i5)') i
             call print(mf%fbs(i), str = fn, unit = unit, all = all, data = data, &
                  skip = unit_get_skip(skip) + 2)
          end do
       end if
       call parallel_barrier()
    end do
  end subroutine imultifab_print
  subroutine lmultifab_print(mf, str, unit, all, data, skip)
    use bl_IO_module
    type(lmultifab), intent(in) :: mf
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: all, data
    integer, intent(in), optional :: skip
    integer :: i, ii
    integer :: un
    character(len=5) :: fn
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un, fmt='("LMULTIFAB")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(": ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" DIM     = ",i2)') mf%dim
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NC      = ",i2)') mf%nc
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NG      = ",i2)') mf%ng
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NODAL   = ",3(L2,1X))') mf%nodal
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NBOXES  = ",i2)') mf%nboxes
    do ii = 0, parallel_nprocs()
       if ( ii == parallel_myproc() ) then
          do i = 1, mf%nboxes; if ( remote(mf,i) ) cycle
             write(unit=fn, fmt='(i5)') i
             call print(mf%fbs(i), str = fn, unit = unit, all = all, data = data, &
                  skip = unit_get_skip(skip) + 2)
          end do
       end if
       call parallel_barrier()
    end do
  end subroutine lmultifab_print

  subroutine multifab_copy_c(mf1, targ, mf2, src, nc, all)
    type(multifab), intent(inout) :: mf1
    type(multifab), intent(in)  :: mf2
    logical, intent(in), optional :: all
    integer, intent(in) :: targ, src
    integer, intent(in), optional :: nc
    real(dp_t), pointer :: mp1(:,:,:,:)
    real(dp_t), pointer :: mp2(:,:,:,:)
    logical :: lall
    integer :: i, n, lnc

    lnc = 1;     if ( present(nc) ) lnc = nc
    lall = .false.; if ( present(all) ) lall = all
    if ( mf1%la == mf2%la ) then
       !$OMP PARALLEL DO PRIVATE(i,mp1,mp2) SHARED(lall)
       do n = 0, lnc-1
          do i = 1, mf1%nboxes
             if ( multifab_remote(mf1,i) ) cycle
             if ( lall ) then
                mp1 => dataptr(mf1, i, get_pbox(mf1, i), targ+n)
                mp2 => dataptr(mf2, i, get_pbox(mf2, i), src +n)
             else
                mp1 => dataptr(mf1, i, get_ibox(mf1, i), targ+n)
                mp2 => dataptr(mf2, i, get_ibox(mf2, i), src +n)
             end if
             mp1 = mp2
          end do
       end do
       !$OMP END PARALLEL DO
    else
       call easy()
    end if

  contains

    subroutine easy()
      type(box) :: bx, abx
      real(dp_t), pointer :: p1(:,:,:,:)
      real(dp_t), pointer :: p2(:,:,:,:)
      integer, parameter :: tag = 1102
      integer :: i, j, proc
      do n = 0, lnc - 1
         do i = 1, mf1%nboxes
            if ( lall ) then
               bx = get_pbox(mf1, i)
            else
               bx = get_ibox(mf1, i)
            end if
            do j = 1, mf2%nboxes
               if ( remote(mf1,i) .and. remote(mf2,j) ) cycle
               abx = intersection(bx, get_ibox(mf2, j))
               if ( .not. empty(abx) ) then
                  if ( local(mf1,i) .and. local(mf2,j) ) then
                     p1 => dataptr(mf1, i, abx, targ+n)
                     p2 => dataptr(mf2, j, abx, src +n)
                     p1 = p2
                  else if ( local(mf2,j) ) then ! must send
                     p2 => dataptr(mf2, j, abx, src+n)
                     proc = get_proc(mf2%la, i)
                     call parallel_send(p2, proc, tag)
                  else if ( local(mf1,i) ) then ! must recv
                     p1 => dataptr(mf1, i, abx, targ+n)
                     proc = get_proc(mf1%la, j)
                     call parallel_recv(p1, proc, tag)
                  end if
               end if
            end do
         end do
      end do

    end subroutine easy

  end subroutine multifab_copy_c
  subroutine imultifab_copy_c(mf1, targ, mf2, src, nc, all)
    type(imultifab), intent(inout) :: mf1
    type(imultifab), intent(in)  :: mf2
    logical, intent(in), optional :: all
    integer, intent(in) :: targ, src
    integer, intent(in), optional :: nc
    integer, pointer :: mp1(:,:,:,:)
    integer, pointer :: mp2(:,:,:,:)
    logical :: lall
    integer :: i, n, lnc

    lnc = 1;     if ( present(nc) ) lnc = nc
    lall = .false.; if ( present(all) ) lall = all
    if ( mf1%la == mf2%la ) then
       !$OMP PARALLEL DO PRIVATE(i,mp1,mp2)
       do n = 0, lnc-1
          do i = 1, mf1%nboxes
             if ( remote(mf1,i) ) cycle
             if ( lall ) then
                mp1 => dataptr(mf1, i, get_pbox(mf1, i), targ + n)
                mp2 => dataptr(mf2, i, get_pbox(mf2, i), src  + n)
             else
                mp1 => dataptr(mf1, i, get_ibox(mf1, i), targ + n)
                mp2 => dataptr(mf2, i, get_ibox(mf2, i), src  + n)
             end if
             mp1 = mp2
          end do
       end do
       !$OMP END PARALLEL DO
    else
       call easy()
    end if

  contains

    subroutine easy()
      type(box) :: bx, abx
      integer, pointer :: p1(:,:,:,:)
      integer, pointer :: p2(:,:,:,:)
      integer, parameter :: tag = 1102
      integer :: i, j, proc
      do n = 0, lnc-1
         do i = 1, mf1%nboxes
            if ( lall ) then
               bx = get_pbox(mf1, i)
            else
               bx = get_ibox(mf1, i)
            end if
            do j = 1, mf2%nboxes
               abx = intersection(bx, get_ibox(mf2, j))
               if ( .not. empty(abx) ) then
                  if ( local(mf1,i) .and. local(mf2,j) ) then
                     p1 => dataptr(mf1, i, abx, targ + n)
                     p2 => dataptr(mf2, j, abx, src  + n)
                     p1 = p2
                  else if ( local(mf2,j) ) then ! must send
                     p2 => dataptr(mf2, j, abx, src + n)
                     proc = get_proc(mf2%la, i)
                     call parallel_send(p2, proc, tag)
                  else if ( local(mf1,i) ) then ! must recv
                     p1 => dataptr(mf1, i, abx, targ + n)
                     proc = get_proc(mf1%la, j)
                     call parallel_recv(p1, proc, tag)
                  end if
               end if
            end do
         end do
      end do
    end subroutine easy

  end subroutine imultifab_copy_c
  subroutine multifab_copy(mf1, mf2, all)
    type(multifab), intent(inout) :: mf1
    type(multifab), intent(in)  :: mf2
    logical, intent(in), optional :: all
    real(dp_t), pointer :: mp1(:,:,:,:)
    real(dp_t), pointer :: mp2(:,:,:,:)
    logical :: lall
    integer :: i

    lall = .false.; if ( present(all) ) lall = all
    if ( mf1%la == mf2%la ) then
       !$OMP PARALLEL DO PRIVATE(i,mp1,mp2) SHARED(lall)
       do i = 1, mf1%nboxes
          if ( multifab_remote(mf1,i) ) cycle
          if ( lall ) then
             mp1 => dataptr(mf1, i, get_pbox(mf1, i))
             mp2 => dataptr(mf2, i, get_pbox(mf2, i))
          else
             mp1 => dataptr(mf1, i, get_ibox(mf1, i))
             mp2 => dataptr(mf2, i, get_ibox(mf2, i))
          end if
          mp1 = mp2
       end do
       !$OMP END PARALLEL DO
    else
       call easy()
    end if

  contains

    subroutine easy()
      type(box) :: bx, abx
      real(dp_t), pointer :: p1(:,:,:,:)
      real(dp_t), pointer :: p2(:,:,:,:)
      integer, parameter :: tag = 1102
      integer i, j, proc
      do i = 1, mf1%nboxes
         if ( lall ) then
            bx = get_pbox(mf1, i)
         else
            bx = get_ibox(mf1, i)
         end if
         do j = 1, mf2%nboxes
            if ( remote(mf1,i) .and. remote(mf2,j) ) cycle
            abx = intersection(bx, get_ibox(mf2, j))
            if ( .not. empty(abx) ) then
               if ( local(mf1,i) .and. local(mf2,j) ) then
                  p1 => dataptr(mf1, i, abx)
                  p2 => dataptr(mf2, j, abx)
                  p1 = p2
               else if ( local(mf2,j) ) then ! must send
                  p2 => dataptr(mf2, j, abx)
                  proc = get_proc(mf2%la, i)
                  call parallel_send(p2, proc, tag)
               else if ( local(mf1,i) ) then ! must recv
                  p1 => dataptr(mf1, i, abx)
                  proc = get_proc(mf1%la, j)
                  call parallel_recv(p1, proc, tag)
               end if
            end if
         end do
      end do

    end subroutine easy

  end subroutine multifab_copy
  subroutine imultifab_copy(mf1, mf2, all)
    type(imultifab), intent(inout) :: mf1
    type(imultifab), intent(in)  :: mf2
    logical, intent(in), optional :: all
    integer, pointer :: mp1(:,:,:,:)
    integer, pointer :: mp2(:,:,:,:)
    logical :: lall
    integer :: i
    lall = .false.; if ( present(all) ) lall = all
    if ( mf1%la == mf2%la ) then
       !$OMP PARALLEL DO PRIVATE(i,mp1,mp2)
       do i = 1, mf1%nboxes
          if ( remote(mf1,i) ) cycle
          if ( lall ) then
             mp1 => dataptr(mf1, i, get_pbox(mf1, i))
             mp2 => dataptr(mf2, i, get_pbox(mf2, i))
          else
             mp1 => dataptr(mf1, i, get_ibox(mf1, i))
             mp2 => dataptr(mf2, i, get_ibox(mf2, i))
          end if
          mp1 = mp2
       end do
       !$OMP END PARALLEL DO
    else
       call easy()
    end if

  contains

    subroutine easy()
      type(box) :: bx, abx
      integer, pointer :: p1(:,:,:,:)
      integer, pointer :: p2(:,:,:,:)
      integer, parameter :: tag = 1102
      integer :: i, j, proc
      do i = 1, mf1%nboxes
         if ( lall ) then
            bx = get_pbox(mf1, i)
         else
            bx = get_ibox(mf1, i)
         end if
         do j = 1, mf2%nboxes
            abx = intersection(bx, get_ibox(mf2, j))
            if ( .not. empty(abx) ) then
               if ( local(mf1,i) .and. local(mf2,j) ) then
                  p1 => dataptr(mf1, i, abx)
                  p2 => dataptr(mf2, j, abx)
                  p1 = p2
               else if ( local(mf2,j) ) then ! must send
                  p2 => dataptr(mf2, j, abx)
                  proc = get_proc(mf2%la, i)
                  call parallel_send(p2, proc, tag)
               else if ( local(mf1,i) ) then ! must recv
                  p1 => dataptr(mf1, i, abx)
                  proc = get_proc(mf1%la, j)
                  call parallel_recv(p1, proc, tag)
               end if
            end if
         end do
      end do

    end subroutine easy

  end subroutine imultifab_copy

  subroutine build_nodal_dot_mask(mask, mf)
    type(multifab), intent(in)  :: mf
    type(multifab), intent(out) :: mask
    integer :: i, d
    type(box) :: full_box, shrunk_box, inner_box

    call build(mask, mf%la, 1, 0, mf%nodal)
    do i = 1, mf%nboxes
       if ( multifab_remote(mf,i) ) cycle
       full_box = get_ibox(mf,i)

       select case ( full_box%dim )
       case (1)
          call setval(mask%fbs(i), 0.5_dp_t, full_box)
       case (2)
          call setval(mask%fbs(i), 0.25_dp_t, full_box)
          do d = 1,2
             shrunk_box = grow(grow(full_box,-1,d,-1),-1,d,+1)
             call setval(mask%fbs(i), 0.5_dp_t, shrunk_box)
          end do
       case (3)
          call setval(mask%fbs(i), 0.125_dp_t, full_box)

          do d = 1,3
             shrunk_box = grow(grow(full_box,-1,d,-1),-1,d,+1)
             call setval(mask%fbs(i), 0.25_dp_t, shrunk_box)
          end do

          d = 1
          shrunk_box = grow(grow(full_box,-1,d,-1),-1,d,+1)
          d = 2
          shrunk_box = grow(grow(shrunk_box,-1,d,-1),-1,d,+1)
          call setval(mask%fbs(i), 0.5_dp_t, shrunk_box)

          d = 1
          shrunk_box = grow(grow(full_box,-1,d,-1),-1,d,+1)
          d = 3
          shrunk_box = grow(grow(shrunk_box,-1,d,-1),-1,d,+1)
          call setval(mask%fbs(i), 0.5_dp_t, shrunk_box)

          d = 2
          shrunk_box = grow(grow(full_box,-1,d,-1),-1,d,+1)
          d = 3
          shrunk_box = grow(grow(shrunk_box,-1,d,-1),-1,d,+1)
          call setval(mask%fbs(i), 0.5_dp_t, shrunk_box)
       end select
       inner_box = grow(full_box,-1)
       call setval(mask%fbs(i),1.0_dp_t,inner_box)
    end do
  end subroutine build_nodal_dot_mask

  function multifab_dot_cc(mf, comp, mf1, comp1, all, mask) result(r)
    real(dp_t) :: r
    logical, intent(in), optional :: all
    type(multifab), intent(in) :: mf
    type(multifab), intent(in)  :: mf1
    type(multifab) :: tmask
    type(lmultifab), optional, intent(in) :: mask
    integer, intent(in) :: comp, comp1
    real(dp_t), pointer :: mp(:,:,:,:)
    real(dp_t), pointer :: mp1(:,:,:,:)
    real(dp_t), pointer :: ma(:,:,:,:)
    logical, pointer :: lmp(:,:,:,:)
    real(dp_t) :: r1
    integer :: i
    logical :: lall

    r1 = 0_dp_t
    if ( present(mask) ) then
       if ( present(all) ) call bl_error('DONT SAY ALL IN MASKED FUNCTION DOT')
       if ( ncomp(mask) /= 1 ) call bl_error('Mask array is multicomponent')
       if ( cell_centered_q(mf) ) then
          !$OMP PARALLEL DO PRIVATE(i,mp,mp1) REDUCTION(+:r1)
          do i = 1, mf%nboxes
             if ( multifab_remote(mf,i) ) cycle
             mp => dataptr(mf, i, get_ibox(mf, i), comp)
             mp1 => dataptr(mf1, i, get_ibox(mf1, i), comp1)
             lmp => dataptr(mask, i, get_ibox(mask, i), 1)
             r1 = r1 + sum(mp*mp1,mask = lmp)
          end do
          !$OMP END PARALLEL DO
       else if ( nodal_q(mf) ) then
          call build_nodal_dot_mask(tmask, mf)
          do i = 1, mf%nboxes
             if ( multifab_remote(mf,i) ) cycle

             mp => dataptr(mf, i, get_ibox(mf, i), comp)
             mp1 => dataptr(mf1, i, get_ibox(mf1, i), comp1)
             ma => dataptr(tmask, i, get_ibox(tmask, i))

             r1 = r1 + sum(ma*mp*mp1)
          end do
          call destroy(tmask)
       else
          call bl_error("MULTIFAB_DOT_CC, failes when not nodal or cell-centered, can be fixed")
       end if
    else
       lall = .FALSE.; if ( present(all) ) lall = all
       if ( cell_centered_q(mf) ) then
          !$OMP PARALLEL DO PRIVATE(i,mp,mp1) REDUCTION(+:r1)
          do i = 1, mf%nboxes
             if ( multifab_remote(mf,i) ) cycle
             if ( lall ) then
                mp => dataptr(mf, i, get_pbox(mf, i), comp)
                mp1 => dataptr(mf1, i, get_pbox(mf1, i), comp1)
             else
                mp => dataptr(mf, i, get_ibox(mf, i), comp)
                mp1 => dataptr(mf1, i, get_ibox(mf1, i), comp1)
             end if
             r1 = r1 + sum(mp*mp1)
          end do
          !$OMP END PARALLEL DO
       else if ( nodal_q(mf) ) then
          if ( lall ) call bl_error('DONT SAY ALL IN NODAL FUNCTION DOT')
          call build_nodal_dot_mask(tmask, mf)
          do i = 1, mf%nboxes
             if ( multifab_remote(mf,i) ) cycle

             mp => dataptr(mf, i, get_ibox(mf, i), comp)
             mp1 => dataptr(mf1, i, get_ibox(mf1, i), comp1)
             ma => dataptr(tmask, i, get_ibox(tmask, i))

             r1 = r1 + sum(ma*mp*mp1)
          end do
          call destroy(tmask)
       else
          call bl_error("MULTIFAB_DOT_CC, failes when not nodal or cell-centered, can be fixed")
       end if
    end if
    call parallel_reduce(r,r1,MPI_SUM)
  end function multifab_dot_cc

  function multifab_dot_c(mf, mf1, comp, all) result(r)
    real(dp_t) :: r
    logical, intent(in), optional :: all
    type(multifab), intent(in) :: mf
    type(multifab), intent(in)  :: mf1
    type(multifab) :: mask
    integer, intent(in) :: comp
    real(dp_t), pointer :: mp(:,:,:,:)
    real(dp_t), pointer :: mp1(:,:,:,:)
    real(dp_t), pointer :: ma(:,:,:,:)
    real(dp_t) :: r1
    integer :: i
    logical :: lall
    lall = .FALSE.; if ( present(all) ) lall = all
    r1 = 0_dp_t
    if ( cell_centered_q(mf) ) then
       !$OMP PARALLEL DO PRIVATE(i,mp,mp1) REDUCTION(+:r1)
       do i = 1, mf%nboxes
          if ( multifab_remote(mf,i) ) cycle
          if ( lall ) then
             mp => dataptr(mf, i, get_pbox(mf, i), comp)
             mp1 => dataptr(mf1, i, get_pbox(mf1, i), comp)
          else
             mp => dataptr(mf, i, get_ibox(mf, i), comp)
             mp1 => dataptr(mf1, i, get_ibox(mf1, i), comp)
          end if
          r1 = r1 + sum(mp*mp1)
       end do
       !$OMP END PARALLEL DO
    else if (nodal_q(mf) ) then
       if ( lall ) call bl_error('DONT SAY ALL IN NODAL FUNCTION DOT')
       call build_nodal_dot_mask(mask, mf)
       do i = 1, mf%nboxes
          if ( multifab_remote(mf,i) ) cycle

          mp => dataptr(mf, i, get_ibox(mf, i), comp)
          mp1 => dataptr(mf1, i, get_ibox(mf1, i), comp)
          ma => dataptr(mask, i, get_ibox(mask, i))

          r1 = r1 + sum(ma*mp*mp1)
       end do
       call destroy(mask)
    else
       call bl_error("MULTIFAB_DOT_C fails when not nodal or cell centered, can be fixed")
    end if
    call parallel_reduce(r,r1,MPI_SUM)
  end function multifab_dot_c

  function multifab_dot(mf, mf1, all) result(r)
    real(dp_t) :: r
    logical, intent(in), optional :: all
    type(multifab), intent(in) :: mf
    type(multifab), intent(in)  :: mf1
    type(multifab) :: mask
    real(dp_t), pointer :: mp(:,:,:,:)
    real(dp_t), pointer :: mp1(:,:,:,:)
    real(dp_t), pointer :: ma(:,:,:,:)
    real(dp_t) :: r1
    integer :: i
    logical :: lall
    lall = .FALSE.; if ( present(all) ) lall = all
    r1 = 0_dp_t
    if ( cell_centered_q(mf) ) then
       !$OMP PARALLEL DO PRIVATE(i,mp,mp1) REDUCTION(+:r1)
       do i = 1, mf%nboxes
          if ( multifab_remote(mf,i) ) cycle
          if ( lall ) then
             mp => dataptr(mf, i, get_pbox(mf, i))
             mp1 => dataptr(mf1, i, get_pbox(mf1, i))
          else
             mp => dataptr(mf, i, get_ibox(mf, i))
             mp1 => dataptr(mf1, i, get_ibox(mf1, i))
          end if
          r1 = r1 + sum(mp*mp1)
       end do
       !$OMP END PARALLEL DO
    else if ( nodal_q(mf) ) then
       if ( lall ) call bl_error('DONT SAY ALL IN NODAL FUNCTION DOT')
       call build_nodal_dot_mask(mask, mf)
       do i = 1, mf%nboxes
          if ( multifab_remote(mf,i) ) cycle

          mp => dataptr(mf, i, get_ibox(mf, i))
          mp1 => dataptr(mf1, i, get_ibox(mf1, i))
          ma => dataptr(mask, i, get_ibox(mask, i))

          r1 = r1 + sum(ma*mp*mp1)
       end do
       call destroy(mask)
    else
       call bl_error("MULTIFAB_DOT fails when not nodal or cell centered, can be fixed")
    end if
    call parallel_reduce(r,r1,MPI_SUM)
  end function multifab_dot

  subroutine multifab_rescale_c(mf, c, val, off)
    real(dp_t), intent(in) :: val
    integer, intent(in) :: c
    real(dp_t), intent(in), optional :: off
    type(multifab), intent(inout) :: mf
    real(dp_t), pointer :: mp(:,:,:,:)
    integer :: i
    !$OMP PARALLEL DO PRIVATE(i,mp)
    do i = 1, mf%nboxes
       if ( multifab_remote(mf,i) ) cycle
       mp => dataptr(mf, i, get_ibox(mf, i), c)
       if ( present(off) ) then
          mp = mp*val + off
       else
          mp = mp*val
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_rescale_c
  subroutine multifab_rescale(mf, val, off)
    real(dp_t), intent(in) :: val
    real(dp_t), intent(in), optional :: off
    type(multifab), intent(inout) :: mf
    real(dp_t), pointer :: mp(:,:,:,:)
    integer :: i
    !$OMP PARALLEL DO PRIVATE(i,mp)
    do i = 1, mf%nboxes
       if ( multifab_remote(mf,i) ) cycle
       mp => dataptr(mf, i, get_ibox(mf, i))
       if ( present(off) ) then
          mp = mp*val + off
       else
          mp = mp*val
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_rescale

  subroutine multifab_saxpy_5(a, b1, b, c1, c)
    real(dp_t), intent(in) :: b1, c1
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b, c
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    real(dp_t), pointer :: cp(:,:,:,:)
    integer :: i
    !$OMP PARALLEL DO PRIVATE(i,ap,bp,cp)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       ap => dataptr(a, i, get_ibox(a, i))
       bp => dataptr(b, i, get_ibox(b, i))
       cp => dataptr(c, i, get_ibox(c, i))
       ap = b1*bp + c1*cp
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_saxpy_5

  subroutine multifab_saxpy_4(a, b, c1, c)
    real(dp_t), intent(in) :: c1
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b,c
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    real(dp_t), pointer :: cp(:,:,:,:)
    integer :: i
    !$OMP PARALLEL DO PRIVATE(i,ap,bp,cp)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       ap => dataptr(a, i, get_ibox(a, i))
       bp => dataptr(b, i, get_ibox(b, i))
       cp => dataptr(c, i, get_ibox(c, i))
       ap = bp + c1*cp
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_saxpy_4

  subroutine multifab_saxpy_3(a, b1, b, all)
    real(dp_t), intent(in) :: b1
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    logical, intent(in), optional :: all
    integer :: i
    logical :: lall

    lall = .false.; if ( present(all) ) lall = all

    if (lall) then
       !$OMP PARALLEL DO PRIVATE(i,ap,bp)
       do i = 1, a%nboxes
          if ( multifab_remote(a,i) ) cycle
          ap => dataptr(a,i)
          bp => dataptr(b,i)
          ap = ap + b1*bp
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,ap,bp)
       do i = 1, a%nboxes
          if ( multifab_remote(a,i) ) cycle
          ap => dataptr(a, i, get_ibox(a, i))
          bp => dataptr(b, i, get_ibox(b, i))
          ap = ap + b1*bp
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine multifab_saxpy_3

  subroutine multifab_saxpy_3_c(a, ia, b1, b, all)
    real(dp_t), intent(in) :: b1
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    integer, intent(in) :: ia
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    logical, intent(in), optional :: all
    integer :: i
    logical :: lall

    lall = .false.; if ( present(all) ) lall = all

    if (lall) then
       !$OMP PARALLEL DO PRIVATE(i,ap,bp)
       do i = 1, a%nboxes; if ( multifab_remote(a,i) ) cycle
          ap => dataptr(a,i,ia)
          bp => dataptr(b,i)
          ap = ap + b1*bp
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,ap,bp)
       do i = 1, a%nboxes; if ( multifab_remote(a,i) ) cycle
          ap => dataptr(a, i, get_ibox(a, i), ia)
          bp => dataptr(b, i, get_ibox(b, i))
          ap = ap + b1*bp
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine multifab_saxpy_3_c

  function multifab_norm_l1_c(mf, comp, nc, mask, all) result(r)
    real(dp_t) :: r
    logical, intent(in), optional :: all
    integer, intent(in) :: comp
    integer, intent(in), optional :: nc
    type(multifab), intent(in) :: mf
    logical, pointer :: lp(:,:,:,:)
    type(lmultifab), intent(in), optional :: mask
    real(dp_t), pointer :: mp(:,:,:,:)
    integer :: i
    real(dp_t) :: r1
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0
    if ( present(mask) ) then
       !$OMP PARALLEL DO PRIVATE(i,mp) REDUCTION(+:r1)
       do i = 1, mf%nboxes
          if ( multifab_remote(mf,i) ) cycle
          if ( lall ) then
             mp => dataptr(mf, i, get_pbox(mf, i))
             lp => dataptr(mask, i, get_pbox(mask, i))
          else
             mp => dataptr(mf, i, get_ibox(mf, i))
             lp => dataptr(mask, i, get_ibox(mask, i))
          end if
          r1 = r1 + sum(abs(mp), mask = lp)
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,mp) REDUCTION(+:r1)
       do i = 1, mf%nboxes
          if ( multifab_remote(mf,i) ) cycle
          if ( lall ) then
             mp => dataptr(mf, i, get_pbox(mf, i), comp, nc)
          else
             mp => dataptr(mf, i, get_ibox(mf, i), comp, nc)
          end if
          r1 = r1 + sum(abs(mp))
       end do
       !$OMP END PARALLEL DO
    end if
    call parallel_reduce(r, r1, MPI_SUM)
  end function multifab_norm_l1_c
  function multifab_norm_l1(mf, all) result(r)
    real(dp_t) :: r
    logical, intent(in), optional :: all
    type(multifab), intent(in) :: mf
    real(dp_t), pointer :: mp(:,:,:,:)
    integer :: i
    real(dp_t) :: r1
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0
    !$OMP PARALLEL DO PRIVATE(i,mp) REDUCTION(+:r1)
    do i = 1, mf%nboxes
       if ( multifab_remote(mf,i) ) cycle
       if ( lall ) then
          mp => dataptr(mf, i, get_pbox(mf, i))
       else
          mp => dataptr(mf, i, get_ibox(mf, i))
       end if
       r1 = r1 + sum(abs(mp))
    end do
    !$OMP END PARALLEL DO
    call parallel_reduce(r, r1, MPI_SUM)
  end function multifab_norm_l1

  function multifab_norm_l2_c(mf, comp, nc, mask, all) result(r)
    real(dp_t) :: r
    integer, intent(in) :: comp
    logical, intent(in), optional :: all
    integer, intent(in), optional :: nc
    type(multifab), intent(in) :: mf
    type(lmultifab), intent(in), optional :: mask
    real(dp_t), pointer :: mp(:,:,:,:)
    logical, pointer :: lp(:,:,:,:)
    integer :: i
    real(dp_t) :: r1
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0
    if ( present(mask) ) then
       !$OMP PARALLEL DO PRIVATE(i,mp) REDUCTION(+:r1)
       do i = 1, mf%nboxes
          if ( multifab_remote(mf,i) ) cycle
          if ( lall ) then
             mp => dataptr(mf, i, get_pbox(mf, i), comp, nc)
             lp => dataptr(mask, i, get_pbox(mask, i))
          else
             mp => dataptr(mf, i, get_ibox(mf, i), comp, nc)
             lp => dataptr(mask, i, get_ibox(mask, i))
          end if
          r1 = r1 + sum(mp**2, mask=lp)
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,mp) REDUCTION(+:r1)
       do i = 1, mf%nboxes
          if ( multifab_remote(mf,i) ) cycle
          if ( lall ) then
             mp => dataptr(mf, i, get_pbox(mf, i), comp, nc)
          else
             mp => dataptr(mf, i, get_ibox(mf, i), comp, nc)
          end if
          r1 = r1 + sum(mp**2)
       end do
       !$OMP END PARALLEL DO
    end if
    call parallel_reduce(r, r1, MPI_SUM)
    r = sqrt(r)
  end function multifab_norm_l2_c
  function multifab_norm_l2(mf, mask, all) result(r)
    real(dp_t) :: r
    logical, intent(in), optional :: all
    type(multifab), intent(in) :: mf
    type(lmultifab), intent(in), optional :: mask
    real(dp_t), pointer :: mp(:,:,:,:)
    logical, pointer :: lp(:,:,:,:)
    integer :: i, nc, n
    real(dp_t) :: r1
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0
    nc = mf%nc
    if ( present(mask) ) then
       !$OMP PARALLEL DO PRIVATE(i,mp) REDUCTION(+:r1)
       do i = 1, mf%nboxes
          if ( multifab_remote(mf,i) ) cycle
          if ( lall ) then
             lp => dataptr(mask, i, get_pbox(mask, i))
          else
             lp => dataptr(mask, i, get_ibox(mask, i))
          end if
          do n = 1, nc
             if ( lall ) then
                mp => dataptr(mf, i, get_pbox(mf, i), c=n)
             else
                mp => dataptr(mf, i, get_ibox(mf, i), c=n)
             end if
             r1 = r1 + sum(mp**2, mask = lp)
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,mp) REDUCTION(+:r1)
       do i = 1, mf%nboxes
          if ( multifab_remote(mf,i) ) cycle
          if ( lall ) then
             mp => dataptr(mf, i, get_pbox(mf, i))
          else
             mp => dataptr(mf, i, get_ibox(mf, i))
          end if
          r1 = r1 + sum(mp**2)
       end do
       !$OMP END PARALLEL DO
    end if
    call parallel_reduce(r, r1, MPI_SUM)
    r = sqrt(r)
  end function multifab_norm_l2

  function multifab_norm_inf_c(mf, comp, nc, mask, all) result(r)
    real(dp_t) :: r
    logical, intent(in), optional :: all
    integer, intent(in) :: comp
    integer, intent(in), optional :: nc
    type(lmultifab), intent(in), optional :: mask
    type(multifab), intent(in) :: mf
    logical, pointer :: lp(:,:,:,:)
    real(dp_t), pointer :: mp(:,:,:,:)
    integer :: i, n
    real(dp_t) :: r1
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0
    if ( present(mask) ) then
       !$OMP PARALLEL DO PRIVATE(i,mp) REDUCTION(MAX:r1)
       do i = 1, mf%nboxes
          if ( multifab_remote(mf,i) ) cycle
          if ( lall ) then
             lp => dataptr(mask, i, get_pbox(mask, i))
          else
             lp => dataptr(mask, i, get_ibox(mask, i))
          end if
          do n = comp, comp+nc-1
             if ( lall ) then
                mp => dataptr(mf, i, get_pbox(mf, i), comp)
             else
                mp => dataptr(mf, i, get_ibox(mf, i), comp)
             end if
             r1 = max(r1, maxval(abs(mp), mask = lp))
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,mp) REDUCTION(MAX:r1)
       do i = 1, mf%nboxes
          if ( multifab_remote(mf,i) ) cycle
          if ( lall ) then
             mp => dataptr(mf, i, get_pbox(mf, i), comp, nc)
          else
             mp => dataptr(mf, i, get_ibox(mf, i), comp, nc)
          end if
          r1 = max(r1, maxval(abs(mp)))
       end do
       !$OMP END PARALLEL DO
    end if
    call parallel_reduce(r, r1, MPI_MAX)
  end function multifab_norm_inf_c
  function multifab_norm_inf(mf, mask, all) result(r)
    real(dp_t) :: r
    logical, intent(in), optional :: all
    type(lmultifab), intent(in), optional :: mask
    type(multifab), intent(in) :: mf
    real(dp_t), pointer :: mp(:,:,:,:)
    logical, pointer :: lp(:,:,:,:)
    integer :: i, n, nc
    real(dp_t) :: r1
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0
    nc = mf%nc
    if ( present(mask) ) then
       !$OMP PARALLEL DO PRIVATE(i,mp) REDUCTION(MAX:r1)
       do i = 1, mf%nboxes
          if ( multifab_remote(mf,i) ) cycle
          if ( lall ) then
             lp => dataptr(mask, i, get_pbox(mask, i))
          else
             lp => dataptr(mask, i, get_ibox(mask, i))
          end if
          do n = 1, nc
             if ( lall ) then
                mp => dataptr(mf, i, get_pbox(mf, i), n)
             else
                mp => dataptr(mf, i, get_ibox(mf, i), n)
             end if
             r1 = max(r1, maxval(abs(mp),mask=lp))
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,mp) REDUCTION(MAX:r1)
       do i = 1, mf%nboxes
          if ( multifab_remote(mf,i) ) cycle
          if ( lall ) then
             mp => dataptr(mf, i, get_pbox(mf, i))
          else
             mp => dataptr(mf, i, get_ibox(mf, i))
          end if
          r1 = max(r1, maxval(abs(mp)))
       end do
       !$OMP END PARALLEL DO
    end if
    call parallel_reduce(r, r1, MPI_MAX)
  end function multifab_norm_inf

  function imultifab_norm_inf_c(mf, comp, nc, all) result(r)
    integer :: r
    logical, intent(in), optional :: all
    integer, intent(in) :: comp
    integer, intent(in), optional :: nc
    type(imultifab), intent(in) :: mf
    integer, pointer :: mp(:,:,:,:)
    integer :: i
    integer :: r1
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0
    !$OMP PARALLEL DO PRIVATE(i,mp) REDUCTION(MAX:r1)
    do i = 1, mf%nboxes
       if ( imultifab_remote(mf,i) ) cycle
       if ( lall ) then
          mp => dataptr(mf, i, get_pbox(mf, i), comp, nc)
       else
          mp => dataptr(mf, i, get_ibox(mf, i), comp, nc)
       end if
       r1 = max(r1, maxval(abs(mp)))
    end do
    !$OMP END PARALLEL DO
    call parallel_reduce(r, r1, MPI_MAX)
  end function imultifab_norm_inf_c
  function imultifab_norm_inf(mf, all) result(r)
    integer :: r
    logical, intent(in), optional :: all
    type(imultifab), intent(in) :: mf
    integer, pointer :: mp(:,:,:,:)
    integer :: i
    integer :: r1
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0
    !$OMP PARALLEL DO PRIVATE(i,mp) REDUCTION(MAX:r1)
    do i = 1, mf%nboxes
       if ( imultifab_remote(mf,i) ) cycle
       if ( lall ) then
          mp => dataptr(mf, i, get_pbox(mf, i))
       else
          mp => dataptr(mf, i, get_ibox(mf, i))
       end if
       r1 = max(r1, maxval(abs(mp)))
    end do
    !$OMP END PARALLEL DO
    call parallel_reduce(r, r1, MPI_MAX)
  end function imultifab_norm_inf

  function imultifab_sum_c(mf, comp, nc, all) result(r)
    integer :: r
    logical, intent(in), optional :: all
    integer, intent(in) :: comp
    integer, intent(in), optional :: nc
    type(imultifab), intent(in) :: mf
    integer, pointer :: mp(:,:,:,:)
    integer :: i
    integer :: r1
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0
    !$OMP PARALLEL DO PRIVATE(i,mp) REDUCTION(+:r1)
    do i = 1, mf%nboxes
       if ( imultifab_remote(mf,i) ) cycle
       if ( lall ) then
          mp => dataptr(mf, i, get_pbox(mf, i), comp, nc)
       else
          mp => dataptr(mf, i, get_ibox(mf, i), comp, nc)
       end if
       r1 = r1 + sum(mp)
    end do
    !$OMP END PARALLEL DO
    call parallel_reduce(r, r1, MPI_SUM)
  end function imultifab_sum_c
  function imultifab_sum(mf, all) result(r)
    integer :: r
    logical, intent(in), optional :: all
    type(imultifab), intent(in) :: mf
    integer, pointer :: mp(:,:,:,:)
    integer :: i
    integer :: r1
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0
    !$OMP PARALLEL DO PRIVATE(i,mp) REDUCTION(+:r1)
    do i = 1, mf%nboxes
       if ( imultifab_remote(mf,i) ) cycle
       if ( lall ) then
          mp => dataptr(mf, i, get_pbox(mf, i))
       else
          mp => dataptr(mf, i, get_ibox(mf, i))
       end if
       r1 = r1 + sum(mp)
    end do
    !$OMP END PARALLEL DO
    call parallel_reduce(r, r1, MPI_SUM)
  end function imultifab_sum

  function lmultifab_count(mf, all) result(r)
    integer :: r
    logical, intent(in), optional :: all
    type(lmultifab), intent(in) :: mf
    logical, pointer :: mp(:,:,:,:)
    integer :: i
    integer :: r1
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0
    !$OMP PARALLEL DO PRIVATE(i,mp) REDUCTION(+:r1)
    do i = 1, mf%nboxes
       if ( lmultifab_remote(mf,i) ) cycle
       if ( lall ) then
          mp => dataptr(mf, i, get_pbox(mf, i))
       else
          mp => dataptr(mf, i, get_ibox(mf, i))
       end if
       r1 = r1 + count(mp)
    end do
    !$OMP END PARALLEL DO
    call parallel_reduce(r, r1, MPI_SUM)
  end function lmultifab_count

!   function lmultifab_any(mf, all) result(r)
!     logical :: r
!     logical, intent(in), optional :: all
!     type(lmultifab), intent(in) :: mf
!     logical, pointer :: mp(:,:,:,:)
!     integer :: i
!     logical :: r1
!     logical :: lall
!     lall = .false.; if ( present(all) ) lall = all
!     r1 = .false.
!     !$OMP PARALLEL DO PRIVATE(i,mp) REDUCTION(.and.:r1)
!     do i = 1, mf%nboxes
!        if ( lmultifab_remote(mf,i) ) cycle
!        if ( lall ) then
!           mp => dataptr(mf, i, get_pbox(mf, i))
!        else
!           mp => dataptr(mf, i, get_ibox(mf, i))
!        end if
!        r1 = r1 .and. any(mp)
!     end do
!     !$OMP END PARALLEL DO
!     call parallel_reduce(r, r1, MPI_LAND)
!   end function lmultifab_any

  subroutine multifab_div_div(a, b, all)
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    logical, intent(in), optional :: all
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,ap,bp)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i)
          bp => dataptr(b, i)
       else
          ap => dataptr(a, i, get_ibox(a, i))
          bp => dataptr(b, i, get_ibox(b, i))
       end if
       if ( any(bp == 0.0_dp_t) ) then
          call bl_error("MULTIFAB_DIV_DIV: divide by zero")
       end if
       ap = ap/bp
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_div_div
  subroutine multifab_div_div_s(a, b, all)
    type(multifab), intent(inout) :: a
    real(dp_t), intent(in)  :: b
    logical, intent(in), optional :: all
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    if ( b == 0.0_dp_t ) then
       call bl_error("MULTIFAB_DIV_DIV: divide by zero")
    end if
    !$OMP PARALLEL DO PRIVATE(i,ap)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i)
       else
          ap => dataptr(a, i, get_ibox(a, i))
       end if
       ap = ap/b
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_div_div_s
  subroutine multifab_div_div_c(a, targ, b, src, nc, all)
    integer, intent(in) :: targ, src
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,ap,bp)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i, targ, nc)
          bp => dataptr(b, i, src, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
          bp => dataptr(b, i, get_ibox(b, i), src, nc)
       end if
       if ( any(bp == 0.0_dp_t) ) then
          call bl_error("MULTIFAB_DIV_DIV: divide by zero")
       end if
       ap = ap/bp
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_div_div_c
  subroutine multifab_div_div_s_c(a, targ, b, nc, all)
    integer, intent(in) :: targ
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    type(multifab), intent(inout) :: a
    real(kind=dp_t), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    if ( b == 0.0_dp_t ) then
       call bl_error("MULTIFAB_DIV_DIV: divide by zero")
    end if
    !$OMP PARALLEL DO PRIVATE(i,ap)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i, targ, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
       end if
       ap = ap/b
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_div_div_s_c

  subroutine multifab_mult_mult(a, b, all)
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    logical, intent(in), optional :: all
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,ap,bp)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i)
          bp => dataptr(b, i)
       else
          ap => dataptr(a, i, get_ibox(a, i))
          bp => dataptr(b, i, get_ibox(b, i))
       end if
       ap = ap*bp
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_mult_mult
  subroutine multifab_mult_mult_s(a, b, all)
    type(multifab), intent(inout) :: a
    real(dp_t), intent(in)  :: b
    logical, intent(in), optional :: all
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,ap)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i)
       else
          ap => dataptr(a, i, get_ibox(a, i))
       end if
       ap = ap*b
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_mult_mult_s
  subroutine multifab_mult_mult_c(a, targ, b, src, nc, all)
    integer, intent(in) :: targ, src
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,ap,bp)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i, targ, nc)
          bp => dataptr(b, i, src, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
          bp => dataptr(b, i, get_ibox(b, i), src, nc)
       end if
       ap = ap*bp
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_mult_mult_c
  subroutine multifab_mult_mult_s_c(a, targ, b, nc, all)
    integer, intent(in) :: targ
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    type(multifab), intent(inout) :: a
    real(kind=dp_t), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,ap)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i, targ, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
       end if
       ap = ap*b
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_mult_mult_s_c

  subroutine multifab_sub_sub(a, b, all)
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    logical, intent(in), optional :: all
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,ap,bp)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i)
          bp => dataptr(b, i)
       else
          ap => dataptr(a, i, get_ibox(a, i))
          bp => dataptr(b, i, get_ibox(b, i))
       end if
       ap = ap - bp
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_sub_sub
  subroutine multifab_sub_sub_s(a, b, all)
    type(multifab), intent(inout) :: a
    real(dp_t), intent(in)  :: b
    logical, intent(in), optional :: all
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,ap)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i)
       else
          ap => dataptr(a, i, get_ibox(a, i))
       end if
       ap = ap - b
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_sub_sub_s
  subroutine multifab_sub_sub_c(a, targ, b, src, nc, all)
    integer, intent(in) :: targ, src
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,ap,bp)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i, targ, nc)
          bp => dataptr(b, i, src, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
          bp => dataptr(b, i, get_ibox(b, i), src, nc)
       end if
       ap = ap - bp
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_sub_sub_c
  subroutine multifab_sub_sub_s_c(a, targ, b, nc, all)
    integer, intent(in) :: targ
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    type(multifab), intent(inout) :: a
    real(kind=dp_t), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,ap)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i, targ, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
       end if
       ap = ap - b
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_sub_sub_s_c

  subroutine multifab_plus_plus(a, b, all)
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    logical, intent(in), optional :: all
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,ap,bp)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i)
          bp => dataptr(b, i)
       else
          ap => dataptr(a, i, get_ibox(a, i))
          bp => dataptr(b, i, get_ibox(b, i))
       end if
       ap = ap + bp
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_plus_plus
  subroutine multifab_plus_plus_s(a, b, all)
    type(multifab), intent(inout) :: a
    real(kind=dp_t), intent(in)  :: b
    logical, intent(in), optional :: all
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,ap)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i)
       else
          ap => dataptr(a, i, get_ibox(a, i))
       end if
       ap = ap + b
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_plus_plus_s
  subroutine multifab_plus_plus_c(a, targ, b, src, nc, all)
    integer, intent(in) :: targ, src
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,ap,bp)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i, targ, nc)
          bp => dataptr(b, i, src, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
          bp => dataptr(b, i, get_ibox(b, i), src, nc)
       end if
       ap = ap + bp
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_plus_plus_c
  subroutine multifab_plus_plus_s_c(a, targ, b, nc, all)
    integer, intent(in) :: targ
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    type(multifab), intent(inout) :: a
    real(kind=dp_t), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,ap)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i, targ, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
       end if
       ap = ap + b
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_plus_plus_s_c

  !! Copies valid data of MF onto FB
  !! FB: Array data
  !! LO: array of lo bound of fb
  !! CF: Start component of destination
  !! MF: Src multifab
  !! CM: start component of source
  !! NC: Number of components to copy
  subroutine multifab_fab_copy(fb, lo, cf, mf, cm, nc)
    integer, intent(in) :: cf, cm, nc, lo(:)
    real(dp_t), intent(inout) :: fb(:,:,:,:)
    type(multifab), intent(in) :: mf
    real(kind=dp_t), pointer :: mp(:,:,:,:)
    integer :: i, xo(mf%dim)
    if ( parallel_q() ) then
       call bl_error("MULTIFAB_FAB_COPY: not ready for parallel")
    end if
    if ( size(fb,dim=4) < cf + nc - 1 ) then
       call bl_error("MULTIFAB_FAB_COPY: fb extent to small")
    end if
    if ( ncomp(mf) < cm + nc - 1 ) then
       call bl_error("MULTIFAB_FAB_COPY: mf extent to small")
    end if
    do i = 1, nboxes(mf); if ( remote(mf,i) ) cycle
       mp => dataptr(mf, i, get_ibox(mf,i), cm, nc)
       xo = lwb(get_ibox(mf,i))
       select case ( mf%dim ) 
       case (1)
          call c_1d(fb(:,1,1,cf:cf+nc-1), lo, mp(:,1,1,cm:cm+nc-1), xo)
       case (2)
          call c_2d(fb(:,:,1,cf:cf+nc-1), lo, mp(:,:,1,cm:cm+nc-1), xo)
       case (3)
          call c_3d(fb(:,:,:,cf:cf+nc-1), lo, mp(:,:,:,cm:cm+nc-1), xo)
       end select
    end do
  contains
    subroutine c_1d(f, lo, x, xo)
      integer, intent(in) :: lo(:), xo(:)
      real(kind=dp_t), intent(inout)  :: f(lo(1):,:)
      real(kind=dp_t), intent(in) :: x(xo(1):,:)
      integer :: i
      do i = max(lbound(f,1),lbound(x,1)), min(ubound(f,1),ubound(x,1))
         f(i,:) = x(i,:)
      end do
    end subroutine c_1d
    subroutine c_2d(f, lo, x, xo)
      integer, intent(in) :: lo(:), xo(:)
      real(kind=dp_t), intent(inout)  :: f(lo(1):,lo(2):,:)
      real(kind=dp_t), intent(in) :: x(xo(1):,xo(2):,:)
      integer :: i, j
      do j = max(lbound(f,2),lbound(x,2)), min(ubound(f,2),ubound(x,2))
         do i = max(lbound(f,1),lbound(x,1)), min(ubound(f,1),ubound(x,1))
            f(i,j,:) = x(i,j,:)
         end do
      end do
    end subroutine c_2d
    subroutine c_3d(f, lo, x, xo)
      integer, intent(in) :: lo(:), xo(:)
      real(kind=dp_t), intent(inout) :: f(lo(1):,lo(2):,lo(3):,:)
      real(kind=dp_t), intent(in) :: x(xo(1):,xo(2):,xo(3):,:)
      integer :: i, j, k
      do k = max(lbound(f,3),lbound(x,3)), min(ubound(f,3),ubound(x,3))
         do j = max(lbound(f,2),lbound(x,2)), min(ubound(f,2),ubound(x,2))
            do i = max(lbound(f,1),lbound(x,1)), min(ubound(f,1),ubound(x,1))
               f(i,j,k,:) = x(i,j,k,:)
            end do
         end do
      end do
    end subroutine c_3d
  end subroutine multifab_fab_copy
  
end module multifab_module
