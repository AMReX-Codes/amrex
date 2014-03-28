module multifab_module

  use bl_error_module
  use layout_module
  use fab_module
  use bl_mem_stat_module
  use bl_prof_module

  implicit none

  type multifab
     !private
     integer      :: dim   = 0
     integer      :: nc    = 1
     integer      :: ng    = 0
     type(layout) :: la
     logical,   pointer :: nodal(:) => Null()
     type(fab), pointer :: fbs(:)   => Null()
  end type multifab

  type zmultifab
     !private
     integer      :: dim   = 0
     integer      :: nc    = 1
     integer      :: ng    = 0
     type(layout) :: la
     logical,    pointer :: nodal(:) => Null()
     type(zfab), pointer :: fbs(:)   => Null()
  end type zmultifab

  type imultifab
     !private
     integer      :: dim   = 0
     integer      :: nc    = 1
     integer      :: ng    = 0
     type(layout) :: la
     logical,    pointer :: nodal(:) => Null()
     type(ifab), pointer :: fbs(:)   => Null()
  end type imultifab

  type lmultifab
     !private
     integer      :: dim   = 0
     integer      :: nc    = 1
     integer      :: ng    = 0
     type(layout) :: la
     logical,    pointer :: nodal(:) => Null()
     type(lfab), pointer :: fbs(:)   => Null()
  end type lmultifab

  type mf_fb_data
     integer :: tag = 100
     logical :: sent = .false.
     logical :: rcvd = .false.
     integer, pointer :: send_request(:) => Null()
     integer, pointer :: recv_request(:) => Null()
     real(dp_t), pointer :: send_buffer(:) => Null()
     real(dp_t), pointer :: recv_buffer(:) => Null()
  end type mf_fb_data

  interface cell_centered_q
     module procedure multifab_cell_centered_q
     module procedure imultifab_cell_centered_q
     module procedure lmultifab_cell_centered_q
     module procedure zmultifab_cell_centered_q
  end interface

  interface nodal_flags
     module procedure multifab_nodal_flags
     module procedure imultifab_nodal_flags
     module procedure lmultifab_nodal_flags
     module procedure zmultifab_nodal_flags
  end interface

  interface nodal_q
     module procedure multifab_nodal_q
     module procedure imultifab_nodal_q
     module procedure lmultifab_nodal_q
     module procedure zmultifab_nodal_q
  end interface

  interface built_q
     module procedure multifab_built_q
     module procedure imultifab_built_q
     module procedure lmultifab_built_q
     module procedure zmultifab_built_q
  end interface

  interface build
     module procedure multifab_build
     module procedure multifab_build_copy

     module procedure imultifab_build
     module procedure imultifab_build_copy

     module procedure lmultifab_build
     module procedure lmultifab_build_copy

     module procedure zmultifab_build
     module procedure zmultifab_build_copy
  end interface

  interface local_index
     module procedure  multifab_local_index
     module procedure imultifab_local_index
     module procedure lmultifab_local_index
     module procedure zmultifab_local_index
  end interface

  interface global_index
     module procedure  multifab_global_index
     module procedure imultifab_global_index
     module procedure lmultifab_global_index
     module procedure zmultifab_global_index
  end interface

  interface destroy
     module procedure multifab_destroy
     module procedure imultifab_destroy
     module procedure lmultifab_destroy
     module procedure zmultifab_destroy
  end interface

  interface copy
     module procedure multifab_copy_c
     module procedure multifab_copy

     module procedure imultifab_copy_c
     module procedure imultifab_copy

     module procedure lmultifab_copy_c
     module procedure lmultifab_copy

     module procedure zmultifab_copy_c
     module procedure zmultifab_copy
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

     module procedure zmultifab_dataptr
     module procedure zmultifab_dataptr_bx
     module procedure zmultifab_dataptr_bx_c
     module procedure zmultifab_dataptr_c

  end interface

  interface print
     module procedure multifab_print
     module procedure imultifab_print
     module procedure lmultifab_print
     module procedure zmultifab_print
  end interface

  interface contains_nan
     module procedure multifab_contains_nan_allc
     module procedure multifab_contains_nan_c
     module procedure multifab_contains_nan_bx_c
  end interface contains_nan

  interface contains_inf
     module procedure multifab_contains_inf_allc
     module procedure multifab_contains_inf_c
     module procedure multifab_contains_inf_bx_c
  end interface contains_inf

  interface setval
     module procedure multifab_setval
     module procedure multifab_setval_bx
     module procedure multifab_setval_c
     module procedure multifab_setval_bx_c
     module procedure multifab_setval_ba

     module procedure imultifab_setval
     module procedure imultifab_setval_bx
     module procedure imultifab_setval_c
     module procedure imultifab_setval_bx_c
     module procedure imultifab_setval_ba

     module procedure lmultifab_setval
     module procedure lmultifab_setval_bx
     module procedure lmultifab_setval_c
     module procedure lmultifab_setval_bx_c
     module procedure lmultifab_setval_ba

     module procedure zmultifab_setval
     module procedure zmultifab_setval_bx
     module procedure zmultifab_setval_c
     module procedure zmultifab_setval_bx_c
     module procedure zmultifab_setval_ba
  end interface

  interface get_boxarray
     module procedure multifab_get_boxarray
     module procedure imultifab_get_boxarray
     module procedure lmultifab_get_boxarray
     module procedure zmultifab_get_boxarray
  end interface

  interface get_pbox
     module procedure multifab_get_pbox
     module procedure imultifab_get_pbox
     module procedure lmultifab_get_pbox
     module procedure zmultifab_get_pbox
  end interface

  interface get_ibox
     module procedure multifab_get_ibox
     module procedure imultifab_get_ibox
     module procedure lmultifab_get_ibox
     module procedure zmultifab_get_ibox
  end interface
  
  interface get_box
     module procedure multifab_get_box
     module procedure imultifab_get_box
     module procedure lmultifab_get_box
     module procedure zmultifab_get_box
  end interface

  interface get_dim
     module procedure multifab_get_dim
     module procedure imultifab_get_dim
     module procedure lmultifab_get_dim
     module procedure zmultifab_get_dim
  end interface

  interface nfabs
     module procedure multifab_nfabs
     module procedure imultifab_nfabs
     module procedure lmultifab_nfabs
     module procedure zmultifab_nfabs
  end interface

  interface nghost
     module procedure multifab_nghost
     module procedure imultifab_nghost
     module procedure lmultifab_nghost
     module procedure zmultifab_nghost
  end interface

  interface volume
     module procedure multifab_volume
     module procedure imultifab_volume
     module procedure lmultifab_volume
     module procedure zmultifab_volume
  end interface

  interface get_layout
     module procedure multifab_get_layout
     module procedure imultifab_get_layout
     module procedure lmultifab_get_layout
     module procedure zmultifab_get_layout
  end interface

  interface rescale
     module procedure multifab_rescale_c
     module procedure multifab_rescale
  end interface

  interface saxpy
     module procedure multifab_saxpy_3_c
     module procedure multifab_saxpy_3_cc
     module procedure multifab_saxpy_3
     module procedure multifab_saxpy_4
     module procedure multifab_saxpy_5
  end interface

  interface fill_boundary
     module procedure multifab_fill_boundary
     module procedure multifab_fill_boundary_c

     module procedure imultifab_fill_boundary
     module procedure imultifab_fill_boundary_c

     module procedure lmultifab_fill_boundary
     module procedure lmultifab_fill_boundary_c

     module procedure zmultifab_fill_boundary
     module procedure zmultifab_fill_boundary_c

  end interface

  interface sum_boundary
     module procedure multifab_sum_boundary
     module procedure multifab_sum_boundary_c

     module procedure lmultifab_sum_boundary
     module procedure lmultifab_sum_boundary_c
  end interface

  interface internal_sync
     module procedure multifab_internal_sync_c
     module procedure multifab_internal_sync

     module procedure lmultifab_internal_sync_c
     module procedure lmultifab_internal_sync
  end interface

  interface ncomp
     module procedure multifab_ncomp
     module procedure imultifab_ncomp
     module procedure lmultifab_ncomp
     module procedure zmultifab_ncomp
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

  interface div
     module procedure multifab_div_s_c
  end interface

  interface div_div
     module procedure multifab_div_div
     module procedure multifab_div_div_s
     module procedure multifab_div_div_c
     module procedure multifab_div_div_s_c
  end interface

  interface mult
     module procedure multifab_mult_s_c
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

  interface min_val
     module procedure multifab_min_c
     module procedure multifab_min
  end interface

  interface max_val
     module procedure multifab_max_c
     module procedure multifab_max
  end interface

  interface count
     module procedure lmultifab_count
  end interface

  interface build_nodal_dot_mask
     module procedure mf_build_nodal_dot_mask
  end interface

  type(mem_stats), private, save ::  multifab_ms
  type(mem_stats), private, save :: imultifab_ms
  type(mem_stats), private, save :: lmultifab_ms
  type(mem_stats), private, save :: zmultifab_ms

  private :: reshape_l_4_1, reshape_l_1_4
  private :: reshape_z_4_1, reshape_z_1_4

  private :: multifab_norm_inf_doit, multifab_saxpy_3_doit
  private :: multifab_div_div_c_doit, multifab_div_div_s_doit
  private :: multifab_mult_mult_c_doit, multifab_mult_mult_s_doit
  private :: multifab_sub_sub_c_doit, multifab_sub_sub_s_doit
  private :: multifab_plus_plus_c_doit, multifab_plus_plus_s_doit, multifab_norm_l2_doit

  public  :: cpy_d, cpy_i, cpy_l, cpy_z
  public  :: reshape_d_4_1, reshape_d_1_4, reshape_i_4_1, reshape_i_1_4

  private :: sum_d, logical_or, mf_internal_sync_fancy, lmf_internal_sync_fancy

  private :: mf_fb_fancy_double, mf_fb_fancy_integer, mf_fb_fancy_logical, mf_fb_fancy_z
  private :: mf_copy_fancy_double, mf_copy_fancy_integer, mf_copy_fancy_logical, mf_copy_fancy_z
  private :: mf_fb_fancy_double_nowait

contains

  pure function multifab_nodal_flags(mf) result(r)
    type(multifab), intent(in) :: mf
    logical :: r(mf%dim)
    r = mf%nodal(1:mf%dim)
  end function multifab_nodal_flags
  pure function imultifab_nodal_flags(mf) result(r)
    type(imultifab), intent(in) :: mf
    logical :: r(mf%dim)
    r = mf%nodal(1:mf%dim)
  end function imultifab_nodal_flags
  pure function zmultifab_nodal_flags(mf) result(r)
    type(zmultifab), intent(in) :: mf
    logical :: r(mf%dim)
    r = mf%nodal(1:mf%dim)
  end function zmultifab_nodal_flags
  pure function lmultifab_nodal_flags(mf) result(r)
    type(lmultifab), intent(in) :: mf
    logical :: r(mf%dim)
    r = mf%nodal(1:mf%dim)
  end function lmultifab_nodal_flags

  subroutine multifab_set_alltoallv(val)
    logical val
  end subroutine multifab_set_alltoallv

  function multifab_local_index(mf,i) result(r)
    integer,         intent(in) :: i
    type(multifab),  intent(in) :: mf
    integer                     :: r
    r = layout_local_index(mf%la,i)
  end function multifab_local_index
  function imultifab_local_index(mf,i) result(r)
    integer,         intent(in) :: i
    type(imultifab), intent(in) :: mf
    integer                     :: r
    r = layout_local_index(mf%la,i)
  end function imultifab_local_index
  function lmultifab_local_index(mf,i) result(r)
    integer,         intent(in) :: i
    type(lmultifab), intent(in) :: mf
    integer                     :: r
    r = layout_local_index(mf%la,i)
  end function lmultifab_local_index
  function zmultifab_local_index(mf,i) result(r)
    integer,         intent(in) :: i
    type(zmultifab), intent(in) :: mf
    integer                     :: r
    r = layout_local_index(mf%la,i)
  end function zmultifab_local_index

  pure function multifab_global_index(mf,i) result(r)
    integer,         intent(in) :: i
    type(multifab),  intent(in) :: mf
    integer                     :: r
    r = layout_global_index(mf%la,i)
  end function multifab_global_index
  pure function imultifab_global_index(mf,i) result(r)
    integer,         intent(in) :: i
    type(imultifab), intent(in) :: mf
    integer                     :: r
    r = layout_global_index(mf%la,i)
  end function imultifab_global_index
  pure function lmultifab_global_index(mf,i) result(r)
    integer,         intent(in) :: i
    type(lmultifab), intent(in) :: mf
    integer                     :: r
    r = layout_global_index(mf%la,i)
  end function lmultifab_global_index
  pure function zmultifab_global_index(mf,i) result(r)
    integer,         intent(in) :: i
    type(zmultifab), intent(in) :: mf
    integer                     :: r
    r = layout_global_index(mf%la,i)
  end function zmultifab_global_index

  pure function multifab_ncomp(mf) result(r)
    integer :: r
    type(multifab), intent(in) :: mf
    r = mf%nc
  end function multifab_ncomp
  pure function imultifab_ncomp(mf) result(r)
    integer :: r
    type(imultifab), intent(in) :: mf
    r = mf%nc
  end function imultifab_ncomp
  pure function lmultifab_ncomp(mf) result(r)
    integer :: r
    type(lmultifab), intent(in) :: mf
    r = mf%nc
  end function lmultifab_ncomp
  pure function zmultifab_ncomp(mf) result(r)
    integer :: r
    type(zmultifab), intent(in) :: mf
    r = mf%nc
  end function zmultifab_ncomp

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
  subroutine zmultifab_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    zmultifab_ms = ms
  end subroutine zmultifab_set_mem_stats

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
  function zmultifab_mem_stats() result(r)
    type(mem_stats) :: r
    r = zmultifab_ms
  end function zmultifab_mem_stats

  pure function multifab_cell_centered_q(mf) result(r)
    logical :: r
    type(multifab), intent(in) :: mf
    r = .not. any(mf%nodal)
  end function multifab_cell_centered_q
  pure function imultifab_cell_centered_q(mf) result(r)
    logical :: r
    type(imultifab), intent(in) :: mf
    r = .not. any(mf%nodal)
  end function imultifab_cell_centered_q
  pure function lmultifab_cell_centered_q(mf) result(r)
    logical :: r
    type(lmultifab), intent(in) :: mf
    r = .not. any(mf%nodal)
  end function lmultifab_cell_centered_q
  pure function zmultifab_cell_centered_q(mf) result(r)
    logical :: r
    type(zmultifab), intent(in) :: mf
    r = .not. any(mf%nodal)
  end function zmultifab_cell_centered_q
  
  pure function multifab_nodal_q(mf) result(r)
    logical :: r
    type(multifab), intent(in) :: mf
    r = all(mf%nodal)
  end function multifab_nodal_q
  pure function imultifab_nodal_q(mf) result(r)
    logical :: r
    type(imultifab), intent(in) :: mf
    r = all(mf%nodal)
  end function imultifab_nodal_q
  pure function lmultifab_nodal_q(mf) result(r)
    logical :: r
    type(lmultifab), intent(in) :: mf
    r = all(mf%nodal)
  end function lmultifab_nodal_q
  pure function zmultifab_nodal_q(mf) result(r)
    logical :: r
    type(zmultifab), intent(in) :: mf
    r = all(mf%nodal)
  end function zmultifab_nodal_q
  
  pure function multifab_built_q(mf) result(r)
    logical :: r
    type(multifab), intent(in) :: mf
    r = mf%dim /= 0
  end function multifab_built_q
  pure function imultifab_built_q(mf) result(r)
    logical :: r
    type(imultifab), intent(in) :: mf
    r = mf%dim /= 0
  end function imultifab_built_q
  pure function lmultifab_built_q(mf) result(r)
    logical :: r
    type(lmultifab), intent(in) :: mf
    r = mf%dim /= 0
  end function lmultifab_built_q
  pure function zmultifab_built_q(mf) result(r)
    logical :: r
    type(zmultifab), intent(in) :: mf
    r = mf%dim /= 0
  end function zmultifab_built_q

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
  function zmultifab_get_layout(mf) result(r)
    type(layout) :: r
    type(zmultifab), intent(in) :: mf
    r = mf%la
  end function zmultifab_get_layout

  subroutine multifab_build_edge(mf, la, nc, ng, dir, stencil)
    type(multifab), intent(  out) :: mf
    type(layout)  , intent(in   ) :: la
    integer       , intent(in   ) :: nc, ng, dir
    logical       , intent(in), optional :: stencil
    logical :: nodal(MAX_SPACEDIM)
    nodal      = .false.
    nodal(dir) = .true.
    call multifab_build(mf, la, nc, ng, nodal, stencil)
  end subroutine multifab_build_edge

  subroutine multifab_build_nodal(mf, la, nc, ng, stencil)
    type(multifab), intent(  out) :: mf
    type(layout)  , intent(in   ) :: la
    integer       , intent(in   ) :: nc, ng
    logical       , intent(in), optional :: stencil
    logical :: nodal(MAX_SPACEDIM)
    nodal = .true.
    call multifab_build(mf, la, nc, ng, nodal, stencil)
  end subroutine multifab_build_nodal

  subroutine multifab_build(mf, la, nc, ng, nodal, stencil)
    type(multifab), intent(out)   :: mf
    type(layout),   intent(in )   :: la
    integer, intent(in), optional :: nc, ng
    logical, intent(in), optional :: nodal(:), stencil
    integer :: i, lnc, lng
    if ( built_q(mf) ) call bl_error("MULTIFAB_BUILD: already built")
    lng = 0; if ( present(ng) ) lng = ng
    lnc = 1; if ( present(nc) ) lnc = nc
    mf%dim = get_dim(la)
    mf%la  = la
    mf%nc  = lnc
    mf%ng  = lng
    allocate(mf%nodal(mf%dim))
    mf%nodal = .False.; if ( present(nodal) ) mf%nodal = nodal(1:mf%dim)
    allocate(mf%fbs(nlocal(mf%la)))

    do i = 1, nlocal(mf%la)
      call build( &
           mf%fbs(i), get_box(mf%la, global_index(mf%la,i)), &
           mf%nc, mf%ng, mf%nodal,  &
           alloc = .true., stencil = stencil)
    end do
    call mem_stats_alloc(multifab_ms, volume(mf, all = .TRUE.))
  end subroutine multifab_build

  subroutine imultifab_build(mf, la, nc, ng, nodal)
    type(imultifab), intent(out)  :: mf
    type(layout),    intent(in )  :: la
    integer, intent(in), optional :: nc, ng
    logical, intent(in), optional :: nodal(:)
    integer :: i, lnc, lng
    if ( built_q(mf) ) call bl_error("IMULTIFAB_BUILD: already built")
    lng = 0; if ( present(ng) ) lng = ng
    lnc = 1; if ( present(nc) ) lnc = nc
    mf%dim = get_dim(la)
    mf%la  = la
    mf%nc  = lnc
    mf%ng  = lng
    allocate(mf%nodal(mf%dim))
    mf%nodal = .False.; if ( present(nodal) ) mf%nodal = nodal(1:mf%dim)
    allocate(mf%fbs(nlocal(mf%la)))
    do i = 1, nlocal(mf%la)
       call build(mf%fbs(i), get_box(mf%la, global_index(mf%la,i)), &
            mf%nc, mf%ng, mf%nodal, alloc = .true.)
    end do
    call mem_stats_alloc(imultifab_ms, volume(mf, all = .TRUE.))
  end subroutine imultifab_build

  subroutine lmultifab_build(mf, la, nc, ng, nodal)
    type(lmultifab), intent(out) :: mf
    type(layout),    intent(in ) :: la
    integer, intent(in),optional :: nc, ng
    logical, intent(in),optional :: nodal(:)
    integer :: i, lnc, lng
    if ( built_q(mf) ) call bl_error("LMULTIFAB_BUILD: already built")
    lng = 0; if ( present(ng) ) lng = ng
    lnc = 1; if ( present(nc) ) lnc = nc
    mf%dim = get_dim(la)
    mf%la  = la
    mf%nc  = lnc
    mf%ng  = lng
    allocate(mf%nodal(mf%dim))
    mf%nodal = .False.; if ( present(nodal) ) mf%nodal = nodal(1:mf%dim)
    allocate(mf%fbs(nlocal(mf%la)))
    do i = 1, nlocal(mf%la)
       call build(mf%fbs(i), get_box(mf%la, global_index(mf%la,i)), &
            mf%nc, mf%ng, mf%nodal, alloc = .true.)
    end do
    call mem_stats_alloc(lmultifab_ms, volume(mf, all = .TRUE.))
  end subroutine lmultifab_build

  subroutine lmultifab_build_edge(mf, la, nc, ng, dir)
    type(lmultifab), intent(  out) :: mf
    type(layout)   , intent(in   ) :: la
    integer        , intent(in   ) :: nc, ng, dir

    logical :: nodal(MAX_SPACEDIM)

    nodal      = .false.
    nodal(dir) = .true.

    call lmultifab_build(mf, la, nc, ng, nodal)

  end subroutine lmultifab_build_edge

  subroutine zmultifab_build(mf, la, nc, ng, nodal)
    type(zmultifab), intent(out) :: mf
    type(layout),    intent(in ) :: la
    integer, intent(in),optional :: nc, ng
    logical, intent(in),optional :: nodal(:)
    integer :: i, lnc, lng
    if ( built_q(mf) ) call bl_error("ZMULTIFAB_BUILD: already built")
    lng = 0; if ( present(ng) ) lng = ng
    lnc = 1; if ( present(nc) ) lnc = nc
    mf%dim = get_dim(la)
    mf%la  = la
    mf%nc  = lnc
    mf%ng  = lng
    allocate(mf%nodal(mf%dim))
    mf%nodal = .False.; if ( present(nodal) ) mf%nodal = nodal(1:mf%dim)
    allocate(mf%fbs(nlocal(mf%la)))
    do i = 1, nlocal(mf%la)
       call build(mf%fbs(i), get_box(mf%la, global_index(mf%la,i)), &
            mf%nc, mf%ng, mf%nodal, alloc = .true.)
    end do
    call mem_stats_alloc(zmultifab_ms, volume(mf, all = .TRUE.))
  end subroutine zmultifab_build

  subroutine multifab_build_copy(m1, m2)
    type(multifab), intent(inout) :: m1
    type(multifab), intent(in) :: m2
    real(dp_t), pointer :: m1p(:,:,:,:)
    real(dp_t), pointer :: m2p(:,:,:,:)
    integer :: i
    if ( built_q(m1) ) call bl_error("MULTIFAB_BUILD_COPY: already built")
    if ( built_q(m1) ) call destroy(m1)
    m1%dim = m2%dim
    m1%la  = m2%la
    m1%nc  = m2%nc
    m1%ng  = m2%ng
    allocate(m1%nodal(m1%dim))
    m1%nodal = m2%nodal
    allocate(m1%fbs(nlocal(m1%la)))
    do i = 1, nlocal(m1%la)
       call build(m1%fbs(i), get_box(m2%fbs(i)), m1%nc, m1%ng, m1%nodal)
       m1p => dataptr(m1,i)
       m2p => dataptr(m2,i)
       call cpy_d(m1p,m2p)
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
    if ( built_q(m1) ) call destroy(m1)
    m1%dim = m2%dim
    m1%la  = m2%la
    m1%nc  = m2%nc
    m1%ng  = m2%ng
    allocate(m1%nodal(m1%dim))
    m1%nodal = m2%nodal
    allocate(m1%fbs(nlocal(m1%la)))
    do i = 1, nlocal(m1%la)
       call build(m1%fbs(i), get_box(m2%fbs(i)), m1%nc, m1%ng, m1%nodal)
       m1p => dataptr(m1,i)
       m2p => dataptr(m2,i)
       call cpy_i(m1p,m2p)
    end do
    call mem_stats_alloc(imultifab_ms, volume(m1, all = .TRUE.))
  end subroutine imultifab_build_copy
  subroutine lmultifab_build_copy(m1, m2)
    type(lmultifab), intent(inout) :: m1
    type(lmultifab), intent(in) :: m2
    logical, pointer :: m1p(:,:,:,:)
    logical, pointer :: m2p(:,:,:,:)
    integer :: i
    if ( built_q(m1) ) call bl_error("LMULTIFAB_BUILD_COPY: already built")
    if ( built_q(m1) ) call destroy(m1)
    m1%dim = m2%dim
    m1%la  = m2%la
    m1%nc  = m2%nc
    m1%ng  = m2%ng
    allocate(m1%nodal(m1%dim))
    m1%nodal = m2%nodal
    allocate(m1%fbs(nlocal(m1%la)))
    do i = 1, nlocal(m1%la)
       call build(m1%fbs(i), get_box(m2%fbs(i)), m1%nc, m1%ng, m1%nodal)
       m1p => dataptr(m1,i)
       m2p => dataptr(m2,i)
       call cpy_l(m1p,m2p)
    end do
    call mem_stats_alloc(lmultifab_ms, volume(m1, all = .TRUE.))
  end subroutine lmultifab_build_copy
  subroutine zmultifab_build_copy(m1, m2)
    type(zmultifab), intent(inout) :: m1
    type(zmultifab), intent(in) :: m2
    complex(dp_t), pointer :: m1p(:,:,:,:)
    complex(dp_t), pointer :: m2p(:,:,:,:)
    integer :: i
    if ( built_q(m1) ) call bl_error("ZMULTIFAB_BUILD_COPY: already built")
    if ( built_q(m1) ) call destroy(m1)
    m1%dim = m2%dim
    m1%la  = m2%la
    m1%nc  = m2%nc
    m1%ng  = m2%ng
    allocate(m1%nodal(m1%dim))
    m1%nodal = m2%nodal
    allocate(m1%fbs(nlocal(m1%la)))
    do i = 1, nlocal(m1%la)
       call build(m1%fbs(i), get_box(m2%fbs(i)), m1%nc, m1%ng, m1%nodal)
       m1p => dataptr(m1,i)
       m2p => dataptr(m2,i)
       call cpy_z(m1p,m2p)
    end do
    call mem_stats_alloc(zmultifab_ms, volume(m1, all = .TRUE.))
  end subroutine zmultifab_build_copy

  subroutine multifab_destroy(mf)
    type(multifab), intent(inout) :: mf
    integer :: i
    call mem_stats_dealloc(multifab_ms, volume(mf, all = .TRUE.))
    do i = 1, nlocal(mf%la)
       call destroy(mf%fbs(i))
    end do
    deallocate(mf%fbs)
    deallocate(mf%nodal)
    mf%dim = 0
    mf%nc  = 0
    mf%ng  = 0
  end subroutine multifab_destroy
  subroutine imultifab_destroy(mf)
    type(imultifab), intent(inout) :: mf
    integer :: i
    call mem_stats_dealloc(imultifab_ms, volume(mf, all = .TRUE.))
    do i = 1, nlocal(mf%la)
       call destroy(mf%fbs(i))
    end do
    deallocate(mf%fbs)
    deallocate(mf%nodal)
    mf%dim = 0
    mf%nc  = 0
    mf%ng  = 0
  end subroutine imultifab_destroy
  subroutine lmultifab_destroy(mf)
    type(lmultifab), intent(inout) :: mf
    integer :: i
    call mem_stats_dealloc(lmultifab_ms, volume(mf, all = .TRUE.))
    do i = 1, nlocal(mf%la)
       call destroy(mf%fbs(i))
    end do
    deallocate(mf%fbs)
    deallocate(mf%nodal)
    mf%dim = 0
    mf%nc  = 0
    mf%ng  = 0
  end subroutine lmultifab_destroy
  subroutine zmultifab_destroy(mf)
    type(zmultifab), intent(inout) :: mf
    integer :: i
    call mem_stats_dealloc(zmultifab_ms, volume(mf, all = .TRUE.))
    do i = 1, nlocal(mf%la)
       call destroy(mf%fbs(i))
    end do
    deallocate(mf%fbs)
    deallocate(mf%nodal)
    mf%dim = 0
    mf%nc  = 0
    mf%ng  = 0
  end subroutine zmultifab_destroy

  function multifab_volume(mf, all) result(r)
    integer(kind=ll_t) :: r
    type(multifab), intent(in) :: mf
    logical, optional :: all
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    if ( lall ) then
       r = 0_ll_t
       do i = 1, nboxes(mf%la)
          r = r + volume(grow(box_nodalize(get_box(mf%la,i),mf%nodal),mf%ng))
       end do
    else
       r = volume(get_boxarray(mf))
    end if
    r = r * mf%nc
  end function multifab_volume
  function imultifab_volume(mf, all) result(r)
    integer(kind=ll_t) :: r
    type(imultifab), intent(in) :: mf
    logical, optional :: all
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    if ( lall ) then
       r = 0_ll_t
       do i = 1, nboxes(mf%la)
          r = r + volume(grow(box_nodalize(get_box(mf%la,i),mf%nodal),mf%ng))
       end do
    else
       r = volume(get_boxarray(mf))
    end if
    r = r * mf%nc
  end function imultifab_volume
  function lmultifab_volume(mf, all) result(r)
    integer(kind=ll_t) :: r
    type(lmultifab), intent(in) :: mf
    logical, optional :: all
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    if ( lall ) then
       r = 0_ll_t
       do i = 1, nboxes(mf%la)
          r = r + volume(grow(box_nodalize(get_box(mf%la,i),mf%nodal),mf%ng))
       end do
    else
       r = volume(get_boxarray(mf))
    end if
    r = r * mf%nc
  end function lmultifab_volume
  function zmultifab_volume(mf, all) result(r)
    integer(kind=ll_t) :: r
    type(zmultifab), intent(in) :: mf
    logical, optional :: all
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    if ( lall ) then
       r = 0_ll_t
       do i = 1, nboxes(mf%la)
          r = r + volume(grow(box_nodalize(get_box(mf%la,i),mf%nodal),mf%ng))
       end do
    else
       r = volume(get_boxarray(mf))
    end if
    r = r * mf%nc
  end function zmultifab_volume

  pure function multifab_get_dim(mf) result(r)
    type(multifab), intent(in) :: mf
    integer :: r
    r = mf%dim
  end function multifab_get_dim
  pure function imultifab_get_dim(mf) result(r)
    type(imultifab), intent(in) :: mf
    integer :: r
    r = mf%dim
  end function imultifab_get_dim
  pure function lmultifab_get_dim(mf) result(r)
    type(lmultifab), intent(in) :: mf
    integer :: r
    r = mf%dim
  end function lmultifab_get_dim
  pure function zmultifab_get_dim(mf) result(r)
    type(zmultifab), intent(in) :: mf
    integer :: r
    r = mf%dim
  end function zmultifab_get_dim

  pure function multifab_nfabs(mf) result(r)
    type(multifab), intent(in) :: mf
    integer :: r
    r = nlocal(mf%la)
  end function multifab_nfabs
  pure function imultifab_nfabs(mf) result(r)
    type(imultifab), intent(in) :: mf
    integer :: r
    r = nlocal(mf%la)
  end function imultifab_nfabs
  pure function lmultifab_nfabs(mf) result(r)
    type(lmultifab), intent(in) :: mf
    integer :: r
    r = nlocal(mf%la)
  end function lmultifab_nfabs
  pure function zmultifab_nfabs(mf) result(r)
    type(zmultifab), intent(in) :: mf
    integer :: r
    r = nlocal(mf%la)
  end function zmultifab_nfabs

  pure function multifab_nghost(mf) result(r)
    type(multifab), intent(in) :: mf
    integer :: r
    r = mf%ng
  end function multifab_nghost
  pure function imultifab_nghost(mf) result(r)
    type(imultifab), intent(in) :: mf
    integer :: r
    r = mf%ng
  end function imultifab_nghost
  pure function lmultifab_nghost(mf) result(r)
    type(lmultifab), intent(in) :: mf
    integer :: r
    r = mf%ng
  end function lmultifab_nghost
  pure function zmultifab_nghost(mf) result(r)
    type(zmultifab), intent(in) :: mf
    integer :: r
    r = mf%ng
  end function zmultifab_nghost

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
  function zmultifab_get_boxarray(mf) result(r)
    type(boxarray) :: r
    type(zmultifab), intent(in) :: mf
    r = get_boxarray(mf%la)
  end function zmultifab_get_boxarray

  pure function multifab_get_box(mf, i) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = get_box(mf%fbs(i))
  end function multifab_get_box
  pure function imultifab_get_box(mf, i) result(r)
    type(imultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = get_box(mf%fbs(i))
  end function imultifab_get_box
  pure function lmultifab_get_box(mf, i) result(r)
    type(lmultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = get_box(mf%fbs(i))
  end function lmultifab_get_box
  pure function zmultifab_get_box(mf, i) result(r)
    type(zmultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = get_box(mf%fbs(i))
  end function zmultifab_get_box

  pure function multifab_get_ibox(mf, i) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = get_ibox(mf%fbs(i))
  end function multifab_get_ibox
  pure function imultifab_get_ibox(mf, i) result(r)
    type(imultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = get_ibox(mf%fbs(i))
  end function imultifab_get_ibox
  pure function lmultifab_get_ibox(mf, i) result(r)
    type(lmultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = get_ibox(mf%fbs(i))
  end function lmultifab_get_ibox
  pure function zmultifab_get_ibox(mf, i) result(r)
    type(zmultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = get_ibox(mf%fbs(i))
  end function zmultifab_get_ibox

  pure function multifab_get_pbox(mf, i) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = get_pbox(mf%fbs(i))
  end function multifab_get_pbox
  pure function imultifab_get_pbox(mf, i) result(r)
    type(imultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = get_pbox(mf%fbs(i))
  end function imultifab_get_pbox
  pure function lmultifab_get_pbox(mf, i) result(r)
    type(lmultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = get_pbox(mf%fbs(i))
  end function lmultifab_get_pbox
  pure function zmultifab_get_pbox(mf, i) result(r)
    type(zmultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = get_pbox(mf%fbs(i))
  end function zmultifab_get_pbox

  function multifab_dataptr(mf, i) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i
    real(dp_t), pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i))
  end function multifab_dataptr
  function imultifab_dataptr(mf, i) result(r)
    type(imultifab), intent(in) :: mf
    integer, intent(in) :: i
    integer, pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i))
  end function imultifab_dataptr
  function lmultifab_dataptr(mf, i) result(r)
    type(lmultifab), intent(in) :: mf
    integer, intent(in) :: i
    logical, pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i))
  end function lmultifab_dataptr
  function zmultifab_dataptr(mf, i) result(r)
    type(zmultifab), intent(in) :: mf
    integer, intent(in) :: i
    complex(dp_t), pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i))
  end function zmultifab_dataptr

  function multifab_dataptr_c(mf, i, c, nc) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i, c
    integer, intent(in), optional :: nc
    real(dp_t), pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i), c, nc)
  end function multifab_dataptr_c
  function imultifab_dataptr_c(mf, i, c, nc) result(r)
    type(imultifab), intent(in) :: mf
    integer, intent(in) :: i, c
    integer, intent(in), optional :: nc
    integer, pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i), c, nc)
  end function imultifab_dataptr_c
  function lmultifab_dataptr_c(mf, i, c, nc) result(r)
    type(lmultifab), intent(in) :: mf
    integer, intent(in) :: i, c
    integer, intent(in), optional :: nc
    logical, pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i), c, nc)
  end function lmultifab_dataptr_c
  function zmultifab_dataptr_c(mf, i, c, nc) result(r)
    type(zmultifab), intent(in) :: mf
    integer, intent(in) :: i, c
    integer, intent(in), optional :: nc
    complex(dp_t), pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i), c, nc)
  end function zmultifab_dataptr_c

  function multifab_dataptr_bx_c(mf, i, bx, c, nc) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i, c
    integer, intent(in), optional :: nc
    type(box), intent(in) :: bx
    real(dp_t), pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i), bx, c, nc)
  end function multifab_dataptr_bx_c
  function imultifab_dataptr_bx_c(mf, i, bx, c, nc) result(r)
    type(imultifab), intent(in) :: mf
    integer, intent(in) :: i, c
    integer, intent(in), optional :: nc
    type(box), intent(in) :: bx
    integer, pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i), bx, c, nc)
  end function imultifab_dataptr_bx_c
  function lmultifab_dataptr_bx_c(mf, i, bx, c, nc) result(r)
    type(lmultifab), intent(in) :: mf
    integer, intent(in) :: i, c
    integer, intent(in), optional :: nc
    type(box), intent(in) :: bx
    logical, pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i), bx, c, nc)
  end function lmultifab_dataptr_bx_c
  function zmultifab_dataptr_bx_c(mf, i, bx, c, nc) result(r)
    type(zmultifab), intent(in) :: mf
    integer, intent(in) :: i, c
    integer, intent(in), optional :: nc
    type(box), intent(in) :: bx
    complex(dp_t), pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i), bx, c, nc)
  end function zmultifab_dataptr_bx_c

  function multifab_dataptr_bx(mf, i, bx) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box), intent(in) :: bx
    real(dp_t), pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i), bx)
  end function multifab_dataptr_bx
  function imultifab_dataptr_bx(mf, i, bx) result(r)
    type(imultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box), intent(in) :: bx
    integer, pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i), bx)
  end function imultifab_dataptr_bx
  function lmultifab_dataptr_bx(mf, i, bx) result(r)
    type(lmultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box), intent(in) :: bx
    logical, pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i), bx)
  end function lmultifab_dataptr_bx
  function zmultifab_dataptr_bx(mf, i, bx) result(r)
    type(zmultifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box), intent(in) :: bx
    complex(dp_t), pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i), bx)
  end function zmultifab_dataptr_bx

  function multifab_contains_nan_c(mf,c,nc) result(r)
    logical                    :: r
    type(multifab), intent(in) :: mf
    integer,        intent(in) :: c, nc
    integer i
    r = .false.
    do i = 1, nfabs(mf)
       r = r .or. contains_nan_c(mf%fbs(i),c,nc)
    enddo
  end function multifab_contains_nan_c

  function multifab_contains_nan_allc(mf) result(r)
    logical                    :: r
    type(multifab), intent(in) :: mf
    r = multifab_contains_nan_c(mf,1,ncomp(mf))
  end function multifab_contains_nan_allc

  function multifab_contains_nan_bx_c(mf,bx,c,nc) result(r)
    logical                    :: r
    type(multifab), intent(in) :: mf
    type(box),      intent(in) :: bx
    integer,        intent(in) :: c, nc
    integer i
    r = .false.
    do i = 1, nfabs(mf)
       r = r .or. contains_nan_bx_c(mf%fbs(i),bx,c,nc)
    enddo
  end function multifab_contains_nan_bx_c

  function multifab_contains_inf_c(mf,c,nc) result(r)
    logical                    :: r
    type(multifab), intent(in) :: mf
    integer,        intent(in) :: c, nc
    integer i
    r = .false.
    do i = 1, nfabs(mf)
       r = r .or. contains_inf_c(mf%fbs(i),c,nc)
    enddo
  end function multifab_contains_inf_c

  function multifab_contains_inf_allc(mf) result(r)
    logical                    :: r
    type(multifab), intent(in) :: mf
    r = multifab_contains_inf_c(mf,1,ncomp(mf))
  end function multifab_contains_inf_allc

  function multifab_contains_inf_bx_c(mf,bx,c,nc) result(r)
    logical                    :: r
    type(multifab), intent(in) :: mf
    type(box),      intent(in) :: bx
    integer,        intent(in) :: c, nc
    integer i
    r = .false.
    do i = 1, nfabs(mf)
       r = r .or. contains_inf_bx_c(mf%fbs(i),bx,c,nc)
    enddo
  end function multifab_contains_inf_bx_c

  subroutine multifab_setval(mf, val, all, allow_empty)
    type(multifab), intent(inout) :: mf
    real(dp_t), intent(in) :: val
    logical, intent(in), optional :: all, allow_empty
    call multifab_setval_c(mf, val, 1, mf%nc, all, allow_empty)
  end subroutine multifab_setval
  subroutine imultifab_setval(mf, val, all)
    type(imultifab), intent(inout) :: mf
    integer, intent(in) :: val
    logical, intent(in), optional :: all
    call imultifab_setval_c(mf, val, 1, mf%nc, all)
  end subroutine imultifab_setval
  subroutine lmultifab_setval(mf, val, all)
    type(lmultifab), intent(inout) :: mf
    logical, intent(in) :: val
    logical, intent(in), optional :: all
    call lmultifab_setval_c(mf, val, 1, mf%nc, all)
  end subroutine lmultifab_setval
  subroutine zmultifab_setval(mf, val, all)
    type(zmultifab), intent(inout) :: mf
    complex(dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    call zmultifab_setval_c(mf, val, 1, mf%nc, all)
  end subroutine zmultifab_setval

  subroutine multifab_setval_bx(mf, val, bx, all)
    type(multifab), intent(inout) :: mf
    type(box), intent(in) :: bx
    real(dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    call multifab_setval_bx_c(mf, val, bx, 1, mf%nc, all)
  end subroutine multifab_setval_bx
  subroutine imultifab_setval_bx(mf, val, bx, all)
    type(imultifab), intent(inout) :: mf
    type(box), intent(in) :: bx
    integer, intent(in) :: val
    logical, intent(in), optional :: all
    call imultifab_setval_bx_c(mf, val, bx, 1, mf%nc, all)
  end subroutine imultifab_setval_bx
  subroutine lmultifab_setval_bx(mf, val, bx, all)
    type(lmultifab), intent(inout) :: mf
    type(box), intent(in) :: bx
    logical, intent(in) :: val
    logical, intent(in), optional :: all
    call lmultifab_setval_bx_c(mf, val, bx, 1, mf%nc, all)
  end subroutine lmultifab_setval_bx
  subroutine zmultifab_setval_bx(mf, val, bx, all)
    type(zmultifab), intent(inout) :: mf
    type(box), intent(in) :: bx
    complex(dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    call zmultifab_setval_bx_c(mf, val, bx, 1, mf%nc, all)
  end subroutine zmultifab_setval_bx

  subroutine multifab_setval_ba(mf, val, ba, all)
    type(multifab), intent(inout) :: mf
    type(boxarray), intent(in) :: ba
    real(dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i, n
    type(box) :: bx, bx1
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,n,bx,bx1)
    do n = 1, nboxes(ba)
       bx = get_box(ba,n)
       do i = 1, nlocal(mf%la)
          if ( lall ) then
             bx1 = intersection(bx, get_pbox(mf, i))
          else
             bx1 = intersection(bx, get_ibox(mf, i))
          end if
          if ( .not. empty(bx1) ) call setval(mf%fbs(i), val, bx1)
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_setval_ba
  subroutine lmultifab_setval_ba(mf, val, ba, all)
    type(lmultifab), intent(inout) :: mf
    type(boxarray), intent(in) :: ba
    logical, intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i, n
    type(box) :: bx, bx1
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,n,bx,bx1)
    do n = 1, nboxes(ba)
       bx = get_box(ba,n)
       do i = 1, nlocal(mf%la)
          if ( lall ) then
             bx1 = intersection(bx, get_pbox(mf, i))
          else
             bx1 = intersection(bx, get_ibox(mf, i))
          end if
          if ( .not. empty(bx1) ) call setval(mf%fbs(i), val, bx1)
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine lmultifab_setval_ba
  subroutine imultifab_setval_ba(mf, val, ba, all)
    type(imultifab), intent(inout) :: mf
    type(boxarray), intent(in) :: ba
    integer, intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i, n
    type(box) :: bx, bx1
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,n,bx,bx1)
    do n = 1, nboxes(ba)
       bx = get_box(ba,n)
       do i = 1, nlocal(mf%la)
          if ( lall ) then
             bx1 = intersection(bx, get_pbox(mf, i))
          else
             bx1 = intersection(bx, get_ibox(mf, i))
          end if
          if ( .not. empty(bx1) ) call setval(mf%fbs(i), val, bx1)
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine imultifab_setval_ba
  subroutine zmultifab_setval_ba(mf, val, ba, all)
    type(zmultifab), intent(inout) :: mf
    type(boxarray), intent(in) :: ba
    complex(dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i, n
    type(box) :: bx, bx1
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    do n = 1, nboxes(ba)
       bx = get_box(ba,n)
       do i = 1, nlocal(mf%la)
          if ( lall ) then
             bx1 = intersection(bx, get_pbox(mf, i))
          else
             bx1 = intersection(bx, get_ibox(mf, i))
          end if
          if ( .not. empty(bx1) ) call setval(mf%fbs(i), val, bx1)
       end do
    end do
  end subroutine zmultifab_setval_ba

  subroutine multifab_setval_bx_c(mf, val, bx, c, nc, all)
    type(multifab), intent(inout) :: mf
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    real(dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    type(box) :: bx1
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,bx1)
    do i = 1, nlocal(mf%la)
       if ( lall ) then
          bx1 = intersection(bx, get_pbox(mf, i))
       else
          bx1 = intersection(bx, get_ibox(mf, i))
       end if
       if ( .not. empty(bx1) ) call setval(mf%fbs(i), val, bx1, c, nc)
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
    !$OMP PARALLEL DO PRIVATE(i,bx1)
    do i = 1, nlocal(mf%la)
       if ( lall ) then
          bx1 = intersection(bx, get_pbox(mf, i))
       else
          bx1 = intersection(bx, get_ibox(mf, i))
       end if
       if ( .not. empty(bx1) ) call setval(mf%fbs(i), val, bx1, c, nc)
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
    !$OMP PARALLEL DO PRIVATE(i,bx1)
    do i = 1, nlocal(mf%la)
       if ( lall ) then
          bx1 = intersection(bx, get_pbox(mf, i))
       else
          bx1 = intersection(bx, get_ibox(mf, i))
       end if
       if ( .not. empty(bx1) ) call setval(mf%fbs(i), val, bx1, c, nc)
    end do
    !$OMP END PARALLEL DO
  end subroutine lmultifab_setval_bx_c
  subroutine zmultifab_setval_bx_c(mf, val, bx, c, nc, all)
    type(zmultifab), intent(inout) :: mf
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    complex(dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    type(box) :: bx1
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    do i = 1, nlocal(mf%la)
       if ( lall ) then
          bx1 = intersection(bx, get_pbox(mf, i))
       else
          bx1 = intersection(bx, get_ibox(mf, i))
       end if
       if ( .not. empty(bx1) ) call setval(mf%fbs(i), val, bx1, c, nc)
    end do
  end subroutine zmultifab_setval_bx_c

  subroutine multifab_setval_c(mf, val, c, nc, all, allow_empty)
    type(multifab), intent(inout) :: mf
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    real(dp_t), intent(in) :: val
    logical, intent(in), optional :: all, allow_empty
    integer :: i
    logical lall, lallow
    type(bl_prof_timer), save :: bpt
    call build(bpt, "mf_setval_c")
    lall = .FALSE.; if ( present(all) ) lall = all
    lallow = .FALSE.; if ( present(allow_empty) ) lallow = allow_empty
    if ( lallow ) then
      !$OMP PARALLEL DO PRIVATE(i)
      do i = 1, nlocal(mf%la)
         if ( lall ) then
            call setval(mf%fbs(i), val, c, nc)
         else
            if ( empty(get_ibox(mf, i)) ) cycle
            call setval(mf%fbs(i), val, get_ibox(mf, i), c, nc)
         end if
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(i)
      do i = 1, nlocal(mf%la)
         if ( lall ) then
            call setval(mf%fbs(i), val, c, nc)
         else
            call setval(mf%fbs(i), val, get_ibox(mf, i), c, nc)
         end if
      end do
      !$OMP END PARALLEL DO
    endif
    call destroy(bpt)
  end subroutine multifab_setval_c
  subroutine imultifab_setval_c(mf, val, c, nc, all)
    type(imultifab), intent(inout) :: mf
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    integer, intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    logical lall
    type(bl_prof_timer), save :: bpt

    call build(bpt, "imf_setval_c")

    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i)
    do i = 1, nlocal(mf%la)
       if ( lall ) then
          call setval(mf%fbs(i), val, c, nc)
       else
          call setval(mf%fbs(i), val, get_ibox(mf, i), c, nc)
       end if
    end do
    !$OMP END PARALLEL DO
    call destroy(bpt)
  end subroutine imultifab_setval_c
  subroutine lmultifab_setval_c(mf, val, c, nc, all)
    type(lmultifab), intent(inout) :: mf
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    logical, intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    logical lall
    type(bl_prof_timer), save :: bpt

    call build(bpt, "lmf_setval_c")
    lall = .FALSE.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i)
    do i = 1, nlocal(mf%la)
       if ( lall ) then
          call setval(mf%fbs(i), val, c, nc)
       else
          call setval(mf%fbs(i), val, get_ibox(mf, i), c, nc)
       end if
    end do
    !$OMP END PARALLEL DO
    call destroy(bpt)
  end subroutine lmultifab_setval_c
  subroutine zmultifab_setval_c(mf, val, c, nc, all)
    type(zmultifab), intent(inout) :: mf
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    complex(dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    logical lall
    lall = .FALSE.; if ( present(all) ) lall = all
    do i = 1, nlocal(mf%la)
       if ( lall ) then
          call setval(mf%fbs(i), val, c, nc)
       else
          call setval(mf%fbs(i), val, get_ibox(mf, i), c, nc)
       end if
    end do
  end subroutine zmultifab_setval_c

  subroutine logical_or(out, in)
     logical, intent(inout) :: out(:,:,:,:)
     logical, intent(in   ) ::  in(:,:,:,:)
     integer                :: i, j, k, n, nx, ny, nz, nc
     !
     ! out = out + in 
     !
     nx = size(out,1)
     ny = size(out,2)
     nz = size(out,3)
     nc = size(out,4)
     do n = 1, nc
        do k = 1, nz
           do j = 1, ny
              do i = 1, nx
                 out(i,j,k,n) = out(i,j,k,n) .or. in(i,j,k,n)
              end do
           end do
        end do
     end do
   end subroutine logical_or

  subroutine sum_d(out, in)
     use bl_types
     real(dp_t), intent(inout) :: out(:,:,:,:)
     real(dp_t), intent(in   ) ::  in(:,:,:,:)
     integer                   :: i, j, k, n, nx, ny, nz, nc
     !
     ! out = out + in 
     !
     nx = size(out,1)
     ny = size(out,2)
     nz = size(out,3)
     nc = size(out,4)
     do n = 1, nc
        do k = 1, nz
           do j = 1, ny
              do i = 1, nx
                 out(i,j,k,n) = out(i,j,k,n) + in(i,j,k,n)
              end do
           end do
        end do
     end do
   end subroutine sum_d

  subroutine cpy_d(out, in, filter)
     real(dp_t), intent(inout) :: out(:,:,:,:)
     real(dp_t), intent(in   ) ::  in(:,:,:,:)
     integer                   :: i, j, k, n, nx, ny, nz, nc
    interface
       subroutine filter(out, in)
         use bl_types
         real(dp_t), intent(inout) :: out(:,:,:,:)
         real(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    if ( present(filter) ) then
       call filter(out,in)
    else
       nx = size(out,1)
       ny = size(out,2)
       nz = size(out,3)
       nc = size(out,4)
       do n = 1, nc
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx
                   out(i,j,k,n) = in(i,j,k,n)
                end do
             end do
          end do
       end do
    end if
  end subroutine cpy_d

  subroutine cpy_i(out, in, filter)
     integer, intent(inout) :: out(:,:,:,:)
     integer, intent(in   ) ::  in(:,:,:,:)
     integer                :: i, j, k, n, nx, ny, nz, nc
    interface
       subroutine filter(out, in)
         integer, intent(inout) :: out(:,:,:,:)
         integer, intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    if ( present(filter) ) then
       call filter(out,in)
    else
       nx = size(out,1)
       ny = size(out,2)
       nz = size(out,3)
       nc = size(out,4)
       do n = 1, nc
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx
                   out(i,j,k,n) = in(i,j,k,n)
                end do
             end do
          end do
       end do
    end if
  end subroutine cpy_i

  subroutine cpy_l(out, in, filter)
     logical, intent(inout) :: out(:,:,:,:)
     logical, intent(in   ) ::  in(:,:,:,:)
     integer                :: i, j, k, n, nx, ny, nz, nc
    interface
       subroutine filter(out, in)
         logical, intent(inout) :: out(:,:,:,:)
         logical, intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    if ( present(filter) ) then
       call filter(out,in)
    else
       nx = size(out,1)
       ny = size(out,2)
       nz = size(out,3)
       nc = size(out,4)
       do n = 1, nc
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx
                   out(i,j,k,n) = in(i,j,k,n)
                end do
             end do
          end do
       end do
    end if
  end subroutine cpy_l

  subroutine cpy_z(out, in, filter)
     complex(dp_t), intent(inout) :: out(:,:,:,:)
     complex(dp_t), intent(in   ) ::  in(:,:,:,:)
     integer                   :: i, j, k, n, nx, ny, nz, nc
    interface
       subroutine filter(out, in)
         use bl_types
         complex(dp_t), intent(inout) :: out(:,:,:,:)
         complex(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    if ( present(filter) ) then
       call filter(out,in)
    else
       nx = size(out,1)
       ny = size(out,2)
       nz = size(out,3)
       nc = size(out,4)
       do n = 1, nc
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx
                   out(i,j,k,n) = in(i,j,k,n)
                end do
             end do
          end do
       end do
    end if
  end subroutine cpy_z

  subroutine reshape_d_4_1(dst,ic,src)
    real(dp_t),intent(in)    :: src(:,:,:,:)
    real(dp_t),intent(inout) :: dst(:)
    integer,intent(in)       :: ic
    integer                  :: i,j,k,n,c,nx,ny,nz,nc
    nx = size(src,1)
    ny = size(src,2)
    nz = size(src,3)
    nc = size(src,4)
    c  = ic
    do n = 1, nc
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                dst(c) = src(i,j,k,n)
                c = c + 1
             end do
          end do
       end do
    end do
  end subroutine reshape_d_4_1

  subroutine reshape_d_1_4(dst,src,ic,sh,filter)
    real(dp_t),intent(in)    :: src(:)
    real(dp_t),intent(inout) :: dst(:,:,:,:)
    integer,intent(in)       :: ic
    integer,intent(in)       :: sh(:)
    integer                  :: i,j,k,n,c,nx,ny,nz,nc
    real(dp_t), allocatable  :: ptmp(:,:,:,:)
    interface
       subroutine filter(out, in)
         use bl_types
         real(dp_t), intent(inout) :: out(:,:,:,:)
         real(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    if ( size(sh) /= 4 ) call bl_error("reshape_d_1_4: how did this happen?")
    c  = ic
    nx = size(dst,1)
    ny = size(dst,2)
    nz = size(dst,3)
    nc = size(dst,4)
    if ( present(filter) ) then
       allocate(ptmp(sh(1),sh(2),sh(3),sh(4)))
       do n = 1, nc
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx
                   ptmp(i,j,k,n) = src(c)
                   c = c + 1
                end do
             end do
          end do
       end do
       call filter(dst, ptmp)
       deallocate(ptmp)
    else
       do n = 1, nc
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx
                   dst(i,j,k,n) = src(c)
                   c = c + 1
                end do
             end do
          end do
       end do
    end if
  end subroutine reshape_d_1_4

  subroutine reshape_i_4_1(dst,ic,src)
    integer,intent(in)    :: src(:,:,:,:)
    integer,intent(inout) :: dst(:)
    integer,intent(in)    :: ic
    integer               :: i,j,k,n,c,nx,ny,nz,nc
    nx = size(src,1)
    ny = size(src,2)
    nz = size(src,3)
    nc = size(src,4)
    c  = ic
    do n = 1, nc
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                dst(c) = src(i,j,k,n)
                c = c + 1
             end do
          end do
       end do
    end do
  end subroutine reshape_i_4_1

  subroutine reshape_i_1_4(dst,src,ic,sh,filter)
    integer,intent(in)    :: src(:)
    integer,intent(inout) :: dst(:,:,:,:)
    integer,intent(in)    :: ic
    integer,intent(in)    :: sh(:)
    integer               :: i,j,k,n,c,nx,ny,nz,nc
    integer, allocatable  :: ptmp(:,:,:,:)
    interface
       subroutine filter(out, in)
         integer, intent(inout) :: out(:,:,:,:)
         integer, intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    if ( size(sh) /= 4 ) call bl_error("reshape_i_1_4: how did this happen?")
    nx = size(dst,1)
    ny = size(dst,2)
    nz = size(dst,3)
    nc = size(dst,4)
    c  = ic
    if ( present(filter) ) then
       allocate(ptmp(sh(1),sh(2),sh(3),sh(4)))
       do n = 1, nc
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx
                   ptmp(i,j,k,n) = src(c)
                   c = c + 1
                end do
             end do
          end do
       end do
       call filter(dst, ptmp)
       deallocate(ptmp)
    else
       do n = 1, nc
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx
                   dst(i,j,k,n) = src(c)
                   c = c + 1
                end do
             end do
          end do
       end do
    end if
  end subroutine reshape_i_1_4

  subroutine reshape_l_4_1(dst,ic,src)
    logical,intent(in)    :: src(:,:,:,:)
    logical,intent(inout) :: dst(:)
    integer,intent(in)    :: ic
    integer               :: i,j,k,n,c,nx,ny,nz,nc
    nx = size(src,1)
    ny = size(src,2)
    nz = size(src,3)
    nc = size(src,4)
    c  = ic
    do n = 1, nc
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                dst(c) = src(i,j,k,n)
                c = c + 1
             end do
          end do
       end do
    end do
  end subroutine reshape_l_4_1

  subroutine reshape_l_1_4(dst,src,ic,sh,filter)
    logical,intent(in)    :: src(:)
    logical,intent(inout) :: dst(:,:,:,:)
    integer,intent(in)    :: ic
    integer,intent(in)    :: sh(:)
    integer               :: i,j,k,n,c,nx,ny,nz,nc
    logical, allocatable  :: ptmp(:,:,:,:)
    interface
       subroutine filter(out, in)
         logical, intent(inout) :: out(:,:,:,:)
         logical, intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    if ( size(sh) /= 4 ) call bl_error("reshape_l_1_4: how did this happen?")
    nx = size(dst,1)
    ny = size(dst,2)
    nz = size(dst,3)
    nc = size(dst,4)
    c  = ic
    if ( present(filter) ) then
       allocate(ptmp(sh(1),sh(2),sh(3),sh(4)))
       do n = 1, nc
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx
                   ptmp(i,j,k,n) = src(c)
                   c = c + 1
                end do
             end do
          end do
       end do
       call filter(dst, ptmp)
       deallocate(ptmp)
    else
       do n = 1, nc
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx
                   dst(i,j,k,n) = src(c)
                   c = c + 1
                end do
             end do
          end do
       end do
    end if
  end subroutine reshape_l_1_4

  subroutine reshape_z_4_1(dst,ic,src)
    complex(dp_t),intent(in)    :: src(:,:,:,:)
    complex(dp_t),intent(inout) :: dst(:)
    integer,intent(in)          :: ic
    integer                     :: i,j,k,n,c,nx,ny,nz,nc
    nx = size(src,1)
    ny = size(src,2)
    nz = size(src,3)
    nc = size(src,4)
    c  = ic
    do n = 1, nc
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                dst(c) = src(i,j,k,n);
                c = c + 1
             end do
          end do
       end do
    end do
  end subroutine reshape_z_4_1

  subroutine reshape_z_1_4(dst,src,ic,sh,filter)
    complex(dp_t),intent(in)    :: src(:)
    complex(dp_t),intent(inout) :: dst(:,:,:,:)
    integer,intent(in)          :: ic
    integer,intent(in)          :: sh(:)
    integer                     :: i,j,k,n,c,nx,ny,nz,nc
    complex(dp_t), allocatable  :: ptmp(:,:,:,:)
    interface
       subroutine filter(out, in)
         use bl_types
         complex(dp_t), intent(inout) :: out(:,:,:,:)
         complex(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    if ( size(sh) /= 4 ) call bl_error("reshape_z_1_4: how did this happen?")
    nx = size(dst,1)
    ny = size(dst,2)
    nz = size(dst,3)
    nc = size(dst,4)
    c  = ic
    if ( present(filter) ) then
       allocate(ptmp(sh(1),sh(2),sh(3),sh(4)))
       do n = 1, nc
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx
                   ptmp(i,j,k,n) = src(c)
                   c = c + 1
                end do
             end do
          end do
       end do
       call filter(dst, ptmp)
       deallocate(ptmp)
    else
       do n = 1, nc
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx
                   dst(i,j,k,n) = src(c)
                   c = c + 1
                end do
             end do
          end do
       end do
    end if
  end subroutine reshape_z_1_4

  subroutine mf_fb_fancy_double(mf, c, nc, ng, lcross, idim)
    type(multifab), intent(inout) :: mf
    integer,        intent(in)    :: c, nc, ng
    logical,        intent(in)    :: lcross
    integer, intent(in), optional :: idim

    real(dp_t), pointer     :: p(:,:,:,:), p1(:,:,:,:), p2(:,:,:,:)
    integer,    allocatable :: rst(:)
    integer,    parameter   :: tag = 1102
    integer                 :: i, ii, jj, np, sh(MAX_SPACEDIM+1)
    type(boxassoc)          :: bxasc
    real(dp_t), allocatable :: g_snd_d(:), g_rcv_d(:)
    logical                 :: cc

    bxasc = layout_boxassoc(mf%la, ng, mf%nodal, lcross, idim)

    cc = multifab_cell_centered_q(mf)
    !$OMP PARALLEL DO PRIVATE(i,ii,jj,p1,p2) if (cc)
    do i = 1, bxasc%l_con%ncpy
       ii  =  local_index(mf,bxasc%l_con%cpy(i)%nd)
       jj  =  local_index(mf,bxasc%l_con%cpy(i)%ns)
       p1  => dataptr(mf%fbs(ii), bxasc%l_con%cpy(i)%dbx, c, nc)
       p2  => dataptr(mf%fbs(jj), bxasc%l_con%cpy(i)%sbx, c, nc)
       call cpy_d(p1,p2)
    end do
    !$OMP END PARALLEL DO

    np = parallel_nprocs()

    if (np == 1) return

    allocate(g_snd_d(nc*bxasc%r_con%svol))
    allocate(g_rcv_d(nc*bxasc%r_con%rvol))

    !$OMP PARALLEL DO PRIVATE(i,p)
    do i = 1, bxasc%r_con%nsnd
       p => dataptr(mf, local_index(mf,bxasc%r_con%snd(i)%ns), bxasc%r_con%snd(i)%sbx, c, nc)
       call reshape_d_4_1(g_snd_d, 1 + nc*bxasc%r_con%snd(i)%pv, p)
    end do
    !$OMP END PARALLEL DO

    allocate(rst(bxasc%r_con%nrp))
    do i = 1, bxasc%r_con%nrp
       rst(i) = parallel_irecv_dv(g_rcv_d(1+nc*bxasc%r_con%rtr(i)%pv:), &
            nc*bxasc%r_con%rtr(i)%sz, bxasc%r_con%rtr(i)%pr, tag)
    end do
    do i = 1, bxasc%r_con%nsp
       call parallel_send_dv(g_snd_d(1+nc*bxasc%r_con%str(i)%pv), &
            nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
    end do
    call parallel_wait(rst)

    !$OMP PARALLEL DO PRIVATE(i,sh,p) if (cc)
    do i = 1, bxasc%r_con%nrcv
       sh = bxasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mf, local_index(mf,bxasc%r_con%rcv(i)%nd), bxasc%r_con%rcv(i)%dbx, c, nc)
       call reshape_d_1_4(p, g_rcv_d, 1 + nc*bxasc%r_con%rcv(i)%pv, sh)
    end do
    !$OMP END PARALLEL DO

  end subroutine mf_fb_fancy_double

  subroutine mf_fb_fancy_double_nowait(mf, fb_data, c, nc, ng, lcross, idim)
    type(multifab), intent(inout) :: mf
    type(mf_fb_data), intent(inout) :: fb_data
    integer,        intent(in)    :: c, nc, ng
    logical,        intent(in)    :: lcross
    integer, intent(in), optional :: idim

    real(dp_t), pointer     :: p(:,:,:,:), p1(:,:,:,:), p2(:,:,:,:)
    integer                 :: i, ii, jj, np, istart, iend, nsize
    type(boxassoc)          :: bxasc

    ! make sure fb_data is clean
    fb_data%sent = .false.
    fb_data%rcvd = .false.
    ! This shouln't happen, unless you forget to call multifab_fill_boundary_finish.
    if (associated(fb_data%send_request)) deallocate(fb_data%send_request)
    if (associated(fb_data%recv_request)) deallocate(fb_data%recv_request)
    if (associated(fb_data%send_buffer )) deallocate(fb_data%send_buffer)
    if (associated(fb_data%recv_buffer )) deallocate(fb_data%recv_buffer)

    bxasc = layout_boxassoc(mf%la, ng, mf%nodal, lcross, idim)

    !$OMP PARALLEL DO PRIVATE(i,ii,jj,p1,p2)
    do i = 1, bxasc%l_con%ncpy
       ii  =  local_index(mf,bxasc%l_con%cpy(i)%nd)
       jj  =  local_index(mf,bxasc%l_con%cpy(i)%ns)
       p1  => dataptr(mf%fbs(ii), bxasc%l_con%cpy(i)%dbx, c, nc)
       p2  => dataptr(mf%fbs(jj), bxasc%l_con%cpy(i)%sbx, c, nc)
       call cpy_d(p1,p2)
    end do
    !$OMP END PARALLEL DO

    np = parallel_nprocs()

    if (np == 1) then
       fb_data%sent = .true. 
       fb_data%rcvd = .true.
       return
    end if

    allocate(fb_data%send_buffer(nc*bxasc%r_con%svol))
    allocate(fb_data%recv_buffer(nc*bxasc%r_con%rvol))

    !$OMP PARALLEL DO PRIVATE(i,p)
    do i = 1, bxasc%r_con%nsnd
       p => dataptr(mf, local_index(mf,bxasc%r_con%snd(i)%ns), bxasc%r_con%snd(i)%sbx, c, nc)
       call reshape_d_4_1(fb_data%send_buffer, 1 + nc*bxasc%r_con%snd(i)%pv, p)
    end do
    !$OMP END PARALLEL DO

    allocate(fb_data%send_request(bxasc%r_con%nsp))
    allocate(fb_data%recv_request(bxasc%r_con%nrp))

    if (bxasc%r_con%nsp .le. 0) then
       fb_data%sent = .true. ! nothing to send
       deallocate(fb_data%send_request)
       deallocate(fb_data%send_buffer)
    end if
    if (bxasc%r_con%nrp .le. 0) then
       fb_data%rcvd = .true. ! nothing to receive
       deallocate(fb_data%recv_request)
       deallocate(fb_data%recv_buffer)
    end if

    do i = 1, bxasc%r_con%nrp
       istart = nc*bxasc%r_con%rtr(i)%pv + 1
       nsize = nc*bxasc%r_con%rtr(i)%sz
       iend = istart + nsize - 1
       fb_data%recv_request(i) = parallel_irecv_dv(fb_data%recv_buffer(istart:iend), &
            nsize, bxasc%r_con%rtr(i)%pr, fb_data%tag)
    end do

    do i = 1, bxasc%r_con%nsp
       istart = nc*bxasc%r_con%str(i)%pv + 1
       nsize = nc*bxasc%r_con%str(i)%sz
       iend = istart + nsize - 1
       fb_data%send_request(i) = parallel_isend_dv(fb_data%send_buffer(istart:iend), &
            nsize, bxasc%r_con%str(i)%pr, fb_data%tag)
    end do

  end subroutine mf_fb_fancy_double_nowait

  subroutine mf_fb_fancy_double_finish(mf, fb_data, c, nc, ng, lcross, idim)
    type(multifab), intent(inout) :: mf
    type(mf_fb_data), intent(inout) :: fb_data
    integer,        intent(in)    :: c, nc, ng
    logical,        intent(in)    :: lcross
    integer, intent(in), optional :: idim

    real(dp_t), pointer :: p(:,:,:,:)
    integer :: i, sh(MAX_SPACEDIM+1)
    type(boxassoc) :: bxasc

    if (fb_data%sent .and. fb_data%rcvd) return

    if (.not. fb_data%sent) then
       call parallel_wait(fb_data%send_request)
       fb_data%sent = .true.
       deallocate(fb_data%send_request)
       deallocate(fb_data%send_buffer)
    end if

    if (.not. fb_data%rcvd) then
       bxasc = layout_boxassoc(mf%la, ng, mf%nodal, lcross, idim)

       call parallel_wait(fb_data%recv_request)

       !$omp parallel do private(i,sh,p)
       do i = 1, bxasc%r_con%nrcv
          sh = bxasc%r_con%rcv(i)%sh
          sh(4) = nc
          p => dataptr(mf, local_index(mf,bxasc%r_con%rcv(i)%nd), bxasc%r_con%rcv(i)%dbx, c, nc)
          call reshape_d_1_4(p, fb_data%recv_buffer, 1 + nc*bxasc%r_con%rcv(i)%pv, sh)
       end do
       !$omp end parallel do

       fb_data%rcvd = .true.
       deallocate(fb_data%recv_request)
       deallocate(fb_data%recv_buffer)
    end if

  end subroutine mf_fb_fancy_double_finish

  subroutine mf_fb_fancy_double_waitrecv(mf, fb_data, c, nc, ng, lcross, idim)
    type(multifab), intent(inout) :: mf
    type(mf_fb_data), intent(inout) :: fb_data
    integer,        intent(in)    :: c, nc, ng
    logical,        intent(in)    :: lcross
    integer, intent(in), optional :: idim

    real(dp_t), pointer :: p(:,:,:,:)
    integer :: i, sh(MAX_SPACEDIM+1)
    type(boxassoc) :: bxasc

    if (fb_data%rcvd) return

    if (.not. fb_data%rcvd) then
       bxasc = layout_boxassoc(mf%la, ng, mf%nodal, lcross, idim)

       call parallel_wait(fb_data%recv_request)

       !$omp parallel do private(i,sh,p)
       do i = 1, bxasc%r_con%nrcv
          sh = bxasc%r_con%rcv(i)%sh
          sh(4) = nc
          p => dataptr(mf, local_index(mf,bxasc%r_con%rcv(i)%nd), bxasc%r_con%rcv(i)%dbx, c, nc)
          call reshape_d_1_4(p, fb_data%recv_buffer, 1 + nc*bxasc%r_con%rcv(i)%pv, sh)
       end do
       !$omp end parallel do

       fb_data%rcvd = .true.
       deallocate(fb_data%recv_request)
       deallocate(fb_data%recv_buffer)
    end if

  end subroutine mf_fb_fancy_double_waitrecv

  subroutine mf_fb_fancy_double_test(mf, fb_data, c, nc, ng, lcross, idim)
    type(multifab), intent(inout) :: mf
    type(mf_fb_data), intent(inout) :: fb_data
    integer,        intent(in)    :: c, nc, ng
    logical,        intent(in)    :: lcross
    integer, intent(in), optional :: idim

    real(dp_t), pointer :: p(:,:,:,:)
    integer :: i, sh(MAX_SPACEDIM+1)
    type(boxassoc) :: bxasc

    if (fb_data%sent .and. fb_data%rcvd) return

    if (.not. fb_data%sent) then
       fb_data%sent = parallel_test(fb_data%send_request)
    end if

    if (fb_data%sent .and. associated(fb_data%send_buffer)) then
       deallocate(fb_data%send_request)
       deallocate(fb_data%send_buffer)
    end if

    if (.not. fb_data%rcvd) then
       fb_data%rcvd = parallel_test(fb_data%recv_request)
    end if

    if (fb_data%rcvd .and. associated(fb_data%recv_buffer)) then
       bxasc = layout_boxassoc(mf%la, ng, mf%nodal, lcross, idim)

       !$omp parallel do private(i,sh,p)
       do i = 1, bxasc%r_con%nrcv
          sh = bxasc%r_con%rcv(i)%sh
          sh(4) = nc
          p => dataptr(mf, local_index(mf,bxasc%r_con%rcv(i)%nd), bxasc%r_con%rcv(i)%dbx, c, nc)
          call reshape_d_1_4(p, fb_data%recv_buffer, 1 + nc*bxasc%r_con%rcv(i)%pv, sh)
       end do
       !$omp end parallel do

       deallocate(fb_data%recv_request)
       deallocate(fb_data%recv_buffer)
    end if

  end subroutine mf_fb_fancy_double_test

  subroutine mf_fb_fancy_integer(mf, c, nc, ng, lcross)
    type(imultifab), intent(inout) :: mf
    integer,         intent(in)    :: c, nc, ng
    logical,         intent(in)    :: lcross

    integer, pointer     :: p(:,:,:,:), p1(:,:,:,:), p2(:,:,:,:)
    integer, allocatable :: rst(:)
    integer, parameter   :: tag = 1102
    integer              :: i, ii, jj, np, sh(MAX_SPACEDIM+1)
    type(boxassoc)       :: bxasc
    integer, allocatable :: g_snd_i(:), g_rcv_i(:)

    bxasc = layout_boxassoc(mf%la, ng, mf%nodal, lcross)

    !$OMP PARALLEL DO PRIVATE(i,ii,jj,p1,p2)
    do i = 1, bxasc%l_con%ncpy
       ii  =  local_index(mf,bxasc%l_con%cpy(i)%nd)
       jj  =  local_index(mf,bxasc%l_con%cpy(i)%ns)
       p1  => dataptr(mf%fbs(ii), bxasc%l_con%cpy(i)%dbx, c, nc)
       p2  => dataptr(mf%fbs(jj), bxasc%l_con%cpy(i)%sbx, c, nc)
       call cpy_i(p1,p2)
    end do
    !$OMP END PARALLEL DO

    np = parallel_nprocs()

    if (np == 1) return

    allocate(g_snd_i(nc*bxasc%r_con%svol))
    allocate(g_rcv_i(nc*bxasc%r_con%rvol))

    !$OMP PARALLEL DO PRIVATE(i,p)
    do i = 1, bxasc%r_con%nsnd
       p => dataptr(mf, local_index(mf,bxasc%r_con%snd(i)%ns), bxasc%r_con%snd(i)%sbx, c, nc)
       call reshape_i_4_1(g_snd_i, 1 + nc*bxasc%r_con%snd(i)%pv, p)
    end do
    !$OMP END PARALLEL DO

    allocate(rst(bxasc%r_con%nrp))
    do i = 1, bxasc%r_con%nrp
       rst(i) = parallel_irecv_iv(g_rcv_i(1+nc*bxasc%r_con%rtr(i)%pv:), &
            nc*bxasc%r_con%rtr(i)%sz, bxasc%r_con%rtr(i)%pr, tag)
    end do
    do i = 1, bxasc%r_con%nsp
       call parallel_send_iv(g_snd_i(1+nc*bxasc%r_con%str(i)%pv), &
            nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
    end do
    call parallel_wait(rst)

    !$OMP PARALLEL DO PRIVATE(i,sh,p)
    do i = 1, bxasc%r_con%nrcv
       sh = bxasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mf, local_index(mf,bxasc%r_con%rcv(i)%nd), bxasc%r_con%rcv(i)%dbx, c, nc)
       call reshape_i_1_4(p, g_rcv_i, 1 + nc*bxasc%r_con%rcv(i)%pv, sh)
    end do
    !$OMP END PARALLEL DO

  end subroutine mf_fb_fancy_integer

  subroutine mf_fb_fancy_logical(mf, c, nc, ng, lcross)
    type(lmultifab), intent(inout) :: mf
    integer,         intent(in)    :: c, nc, ng
    logical,         intent(in)    :: lcross

    logical, pointer     :: p(:,:,:,:), p1(:,:,:,:), p2(:,:,:,:)
    integer, allocatable :: rst(:)
    integer, parameter   :: tag = 1102
    integer              :: i, ii, jj, np, sh(MAX_SPACEDIM+1)
    type(boxassoc)       :: bxasc
    logical, allocatable :: g_snd_l(:), g_rcv_l(:)

    bxasc = layout_boxassoc(mf%la, ng, mf%nodal, lcross)

    !$OMP PARALLEL DO PRIVATE(i,ii,jj,p1,p2)
    do i = 1, bxasc%l_con%ncpy
       ii  =  local_index(mf,bxasc%l_con%cpy(i)%nd)
       jj  =  local_index(mf,bxasc%l_con%cpy(i)%ns)
       p1  => dataptr(mf%fbs(ii), bxasc%l_con%cpy(i)%dbx, c, nc)
       p2  => dataptr(mf%fbs(jj), bxasc%l_con%cpy(i)%sbx, c, nc)
       call cpy_l(p1,p2)
    end do
    !$OMP END PARALLEL DO

    np = parallel_nprocs()

    if (np == 1) return

    allocate(g_snd_l(nc*bxasc%r_con%svol))
    allocate(g_rcv_l(nc*bxasc%r_con%rvol))

    !$OMP PARALLEL DO PRIVATE(i,p)
    do i = 1, bxasc%r_con%nsnd
       p => dataptr(mf, local_index(mf,bxasc%r_con%snd(i)%ns), bxasc%r_con%snd(i)%sbx, c, nc)
       call reshape_l_4_1(g_snd_l, 1 + nc*bxasc%r_con%snd(i)%pv, p)
    end do
    !$OMP END PARALLEL DO

    allocate(rst(bxasc%r_con%nrp))
    do i = 1, bxasc%r_con%nrp
       rst(i) = parallel_irecv_lv(g_rcv_l(1+nc*bxasc%r_con%rtr(i)%pv:), &
            nc*bxasc%r_con%rtr(i)%sz, bxasc%r_con%rtr(i)%pr, tag)
    end do
    do i = 1, bxasc%r_con%nsp
       call parallel_send_lv(g_snd_l(1+nc*bxasc%r_con%str(i)%pv), &
            nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
    end do
    call parallel_wait(rst)

    !$OMP PARALLEL DO PRIVATE(i,sh,p)
    do i = 1, bxasc%r_con%nrcv
       sh = bxasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mf, local_index(mf,bxasc%r_con%rcv(i)%nd), bxasc%r_con%rcv(i)%dbx, c, nc)
       call reshape_l_1_4(p, g_rcv_l, 1 + nc*bxasc%r_con%rcv(i)%pv, sh)
    end do
    !$OMP END PARALLEL DO

  end subroutine mf_fb_fancy_logical

  subroutine mf_fb_fancy_z(mf, c, nc, ng, lcross)
    type(zmultifab), intent(inout) :: mf
    integer,         intent(in)    :: c, nc, ng
    logical,         intent(in)    :: lcross

    complex(dp_t), pointer     :: p(:,:,:,:), p1(:,:,:,:), p2(:,:,:,:)
    integer, allocatable       :: rst(:)
    integer, parameter         :: tag = 1102
    integer                    :: i, ii, jj, sh(MAX_SPACEDIM+1)
    type(boxassoc)             :: bxasc
    complex(dp_t), allocatable :: g_snd_z(:), g_rcv_z(:)

    bxasc = layout_boxassoc(mf%la, ng, mf%nodal, lcross)

    do i = 1, bxasc%l_con%ncpy
       ii  =  local_index(mf,bxasc%l_con%cpy(i)%nd)
       jj  =  local_index(mf,bxasc%l_con%cpy(i)%ns)
       p1  => dataptr(mf%fbs(ii), bxasc%l_con%cpy(i)%dbx, c, nc)
       p2  => dataptr(mf%fbs(jj), bxasc%l_con%cpy(i)%sbx, c, nc)
       call cpy_z(p1,p2)
    end do

    allocate(g_snd_z(nc*bxasc%r_con%svol))
    allocate(g_rcv_z(nc*bxasc%r_con%rvol))

    do i = 1, bxasc%r_con%nsnd
       p => dataptr(mf, local_index(mf,bxasc%r_con%snd(i)%ns), bxasc%r_con%snd(i)%sbx, c, nc)
       call reshape_z_4_1(g_snd_z, 1 + nc*bxasc%r_con%snd(i)%pv, p)
    end do

    allocate(rst(bxasc%r_con%nrp))
    do i = 1, bxasc%r_con%nrp
       rst(i) = parallel_irecv_zv(g_rcv_z(1+nc*bxasc%r_con%rtr(i)%pv:), &
            nc*bxasc%r_con%rtr(i)%sz, bxasc%r_con%rtr(i)%pr, tag)
    end do
    do i = 1, bxasc%r_con%nsp
       call parallel_send_zv(g_snd_z(1+nc*bxasc%r_con%str(i)%pv), &
            nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
    end do
    call parallel_wait(rst)

    do i = 1, bxasc%r_con%nrcv
       sh = bxasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mf, local_index(mf,bxasc%r_con%rcv(i)%nd), bxasc%r_con%rcv(i)%dbx, c, nc)
       call reshape_z_1_4(p, g_rcv_z, 1 + nc*bxasc%r_con%rcv(i)%pv, sh)
    end do

  end subroutine mf_fb_fancy_z

  subroutine multifab_fill_boundary_c(mf, c, nc, ng, cross, idim)
    type(multifab), intent(inout) :: mf
    integer, intent(in)           :: c, nc
    integer, intent(in), optional :: ng, idim
    logical, intent(in), optional :: cross
    integer :: lng
    logical :: lcross
    type(bl_prof_timer), save :: bpt
    lcross  = .false.; if ( present(cross)  ) lcross  = cross
    lng     = mf%ng;   if ( present(ng)     ) lng     = ng
    if ( lng > mf%ng      ) call bl_error("MULTIFAB_FILL_BOUNDARY_C: ng too large", lng)
    if ( mf%nc < (c+nc-1) ) call bl_error('MULTIFAB_FILL_BOUNDARY_C: nc too large', nc)
    if ( present(idim) ) then
       if (idim > 0) lcross = .true. 
    end if
   
    ! If the boxarray is contained in the domain, then this made sense because nothing will
    !  be done if ng = 0.  However, sometimes fillpatch calls this with a boxarray that is 
    !  not contained in the domain, and we need to use fill_boundary to fill regions of the 
    !  boxarray that are "valid" (i.e. not ghost cells) but that are outside the domain.
    ! if ( lng < 1          ) return

    call build(bpt, "mf_fill_boundary_c")
    call mf_fb_fancy_double(mf, c, nc, lng, lcross, idim)
    call destroy(bpt)
  end subroutine multifab_fill_boundary_c

  subroutine multifab_fill_boundary(mf, ng, cross, idim)
    type(multifab), intent(inout) :: mf
    integer, intent(in), optional :: ng, idim
    logical, intent(in), optional :: cross
    call multifab_fill_boundary_c(mf, 1, mf%nc, ng, cross, idim)
  end subroutine multifab_fill_boundary

  subroutine multifab_fill_boundary_nowait_c(mf, fb_data, c, nc, ng, cross, idim)
    type(multifab), intent(inout) :: mf
    type(mf_fb_data), intent(inout) :: fb_data
    integer, intent(in)           :: c, nc
    integer, intent(in), optional :: ng, idim
    logical, intent(in), optional :: cross
    integer :: lng
    logical :: lcross
    type(bl_prof_timer), save :: bpt
    lcross  = .false.; if ( present(cross)  ) lcross  = cross
    lng     = mf%ng;   if ( present(ng)     ) lng     = ng
    if ( lng > mf%ng      ) call bl_error("MULTIFAB_FILL_BOUNDARY_NOWAIT_C: ng too large", lng)
    if ( mf%nc < (c+nc-1) ) call bl_error('MULTIFAB_FILL_BOUNDARY_NOWAIT_C: nc too large', nc)
    if ( present(idim) ) then
       if (idim > 0) lcross = .true. 
    end if
   
    ! If the boxarray is contained in the domain, then this made sense because nothing will
    !  be done if ng = 0.  However, sometimes fillpatch calls this with a boxarray that is 
    !  not contained in the domain, and we need to use fill_boundary to fill regions of the 
    !  boxarray that are "valid" (i.e. not ghost cells) but that are outside the domain.
    ! if ( lng < 1          ) return

    call build(bpt, "mf_fill_boundary_nowait_c")
    call mf_fb_fancy_double_nowait(mf, fb_data, c, nc, lng, lcross, idim)
    call destroy(bpt)
  end subroutine multifab_fill_boundary_nowait_c

  subroutine multifab_fill_boundary_nowait(mf, fb_data, ng, cross, idim)
    type(multifab), intent(inout) :: mf
    type(mf_fb_data), intent(inout) :: fb_data
    integer, intent(in), optional :: ng, idim
    logical, intent(in), optional :: cross
    call multifab_fill_boundary_nowait_c(mf, fb_data, 1, mf%nc, ng, cross, idim)
  end subroutine multifab_fill_boundary_nowait

  subroutine multifab_fill_boundary_finish_c(mf, fb_data, c, nc, ng, cross, idim)
    type(multifab), intent(inout) :: mf
    type(mf_fb_data), intent(inout) :: fb_data
    integer, intent(in)           :: c, nc
    integer, intent(in), optional :: ng, idim
    logical, intent(in), optional :: cross
    integer :: lng
    logical :: lcross
    type(bl_prof_timer), save :: bpt
    lcross  = .false.; if ( present(cross)  ) lcross  = cross
    lng     = mf%ng;   if ( present(ng)     ) lng     = ng
    if ( lng > mf%ng      ) call bl_error("MULTIFAB_FILL_BOUNDARY_FINISH_C: ng too large", lng)
    if ( mf%nc < (c+nc-1) ) call bl_error('MULTIFAB_FILL_BOUNDARY_FINISH_C: nc too large', nc)
    if ( present(idim) ) then
       if (idim > 0) lcross = .true. 
    end if
    call build(bpt, "mf_fill_boundary_finish_c")
    call mf_fb_fancy_double_finish(mf, fb_data, c, nc, lng, lcross, idim)
    call destroy(bpt)
  end subroutine multifab_fill_boundary_finish_c

  subroutine multifab_fill_boundary_waitrecv_c(mf, fb_data, c, nc, ng, cross, idim)
    type(multifab), intent(inout) :: mf
    type(mf_fb_data), intent(inout) :: fb_data
    integer, intent(in)           :: c, nc
    integer, intent(in), optional :: ng, idim
    logical, intent(in), optional :: cross
    integer :: lng
    logical :: lcross
    type(bl_prof_timer), save :: bpt
    lcross  = .false.; if ( present(cross)  ) lcross  = cross
    lng     = mf%ng;   if ( present(ng)     ) lng     = ng
    if ( lng > mf%ng      ) call bl_error("MULTIFAB_FILL_BOUNDARY_WAITRECV_C: ng too large", lng)
    if ( mf%nc < (c+nc-1) ) call bl_error('MULTIFAB_FILL_BOUNDARY_WAITRECV_C: nc too large', nc)
    if ( present(idim) ) then
       if (idim > 0) lcross = .true. 
    end if
    call build(bpt, "mf_fill_boundary_waitrecv_c")
    call mf_fb_fancy_double_waitrecv(mf, fb_data, c, nc, lng, lcross, idim)
    call destroy(bpt)
  end subroutine multifab_fill_boundary_waitrecv_c

  subroutine multifab_fill_boundary_finish(mf, fb_data, ng, cross, idim)
    type(multifab), intent(inout) :: mf
    type(mf_fb_data), intent(inout) :: fb_data
    integer, intent(in), optional :: ng, idim
    logical, intent(in), optional :: cross
    call multifab_fill_boundary_finish_c(mf, fb_data, 1, mf%nc, ng, cross, idim)
  end subroutine multifab_fill_boundary_finish

  subroutine multifab_fill_boundary_waitrecv(mf, fb_data, ng, cross, idim)
    type(multifab), intent(inout) :: mf
    type(mf_fb_data), intent(inout) :: fb_data
    integer, intent(in), optional :: ng, idim
    logical, intent(in), optional :: cross
    call multifab_fill_boundary_waitrecv_c(mf, fb_data, 1, mf%nc, ng, cross, idim)
  end subroutine multifab_fill_boundary_waitrecv

  subroutine multifab_fill_boundary_test_c(mf, fb_data, c, nc, ng, cross, idim)
    type(multifab), intent(inout) :: mf
    type(mf_fb_data), intent(inout) :: fb_data
    integer, intent(in)           :: c, nc
    integer, intent(in), optional :: ng, idim
    logical, intent(in), optional :: cross
    integer :: lng
    logical :: lcross
    type(bl_prof_timer), save :: bpt
    lcross  = .false.; if ( present(cross)  ) lcross  = cross
    lng     = mf%ng;   if ( present(ng)     ) lng     = ng
    if ( lng > mf%ng      ) call bl_error("MULTIFAB_FILL_BOUNDARY_TEST_C: ng too large", lng)
    if ( mf%nc < (c+nc-1) ) call bl_error('MULTIFAB_FILL_BOUNDARY_TEST_C: nc too large', nc)
    if ( present(idim) ) then
       if (idim > 0) lcross = .true. 
    end if
   
    ! If the boxarray is contained in the domain, then this made sense because nothing will
    !  be done if ng = 0.  However, sometimes fillpatch calls this with a boxarray that is 
    !  not contained in the domain, and we need to use fill_boundary to fill regions of the 
    !  boxarray that are "valid" (i.e. not ghost cells) but that are outside the domain.
    ! if ( lng < 1          ) return

    call build(bpt, "mf_fill_boundary_test_c")
    call mf_fb_fancy_double_test(mf, fb_data, c, nc, lng, lcross, idim)
    call destroy(bpt)
  end subroutine multifab_fill_boundary_test_c

  subroutine multifab_fill_boundary_test(mf, fb_data, ng, cross, idim)
    type(multifab), intent(inout) :: mf
    type(mf_fb_data), intent(inout) :: fb_data
    integer, intent(in), optional :: ng, idim
    logical, intent(in), optional :: cross
    call multifab_fill_boundary_test_c(mf, fb_data, 1, mf%nc, ng, cross, idim)
  end subroutine multifab_fill_boundary_test

  subroutine imultifab_fill_boundary_c(mf, c, nc, ng, cross)
    type(imultifab), intent(inout) :: mf
    integer, intent(in)            :: c, nc
    integer, intent(in), optional  :: ng
    logical, intent(in), optional  :: cross
    integer :: lng
    logical :: lcross
    type(bl_prof_timer), save :: bpt
    lcross  = .false.; if ( present(cross)  ) lcross  = cross
    lng     = mf%ng;   if ( present(ng)     ) lng     = ng
    if ( lng > mf%ng )      call bl_error('IMULTIFAB_FILL_BOUNDARY_C: ng too large', lng)
    if ( mf%nc < (c+nc-1) ) call bl_error('IMULTIFAB_FILL_BOUNDARY_C: nc too large', nc)
    if ( lng < 1 ) return
    call build(bpt, "imf_fill_boundary_c")
    call mf_fb_fancy_integer(mf, c, nc, lng, lcross)
    call destroy(bpt)
  end subroutine imultifab_fill_boundary_c

  subroutine imultifab_fill_boundary(mf, ng, cross)
    type(imultifab), intent(inout) :: mf
    integer, intent(in), optional  :: ng
    logical, intent(in), optional  :: cross
    call imultifab_fill_boundary_c(mf, 1, mf%nc, ng, cross)
  end subroutine imultifab_fill_boundary

  subroutine lmultifab_fill_boundary_c(mf, c, nc, ng, cross)
    type(lmultifab), intent(inout) :: mf
    integer, intent(in)            :: c, nc
    integer, intent(in), optional  :: ng
    logical, intent(in), optional  :: cross
    integer :: lng
    logical :: lcross
    type(bl_prof_timer), save :: bpt
    lcross  = .false.; if ( present(cross)  ) lcross  = cross
    lng     = mf%ng;   if ( present(ng)     ) lng     = ng
    if ( lng > mf%ng )      call bl_error('LMULTIFAB_FILL_BOUNDARY_C: ng too large', lng)
    if ( mf%nc < (c+nc-1) ) call bl_error('LMULTIFAB_FILL_BOUNDARY_C: nc too large', nc)
    if ( lng < 1 ) return
    call build(bpt, "lmf_fill_boundary_c")
    call mf_fb_fancy_logical(mf, c, nc, lng, lcross)
    call destroy(bpt)
  end subroutine lmultifab_fill_boundary_c

  subroutine lmultifab_fill_boundary(mf, ng, cross)
    type(lmultifab), intent(inout) :: mf
    integer, intent(in), optional  :: ng
    logical, intent(in), optional  :: cross
    call lmultifab_fill_boundary_c(mf, 1, mf%nc, ng, cross)
  end subroutine lmultifab_fill_boundary

  subroutine zmultifab_fill_boundary_c(mf, c, nc, ng, cross)
    type(zmultifab), intent(inout) :: mf
    integer, intent(in)            :: c, nc
    integer, intent(in), optional  :: ng
    logical, intent(in), optional  :: cross
    integer :: lng
    logical :: lcross
    type(bl_prof_timer), save :: bpt
    lcross  = .false.; if ( present(cross)  ) lcross  = cross
    lng     = mf%ng;   if ( present(ng)     ) lng     = ng
    if ( lng > mf%ng )      call bl_error('ZMULTIFAB_FILL_BOUNDARY_C: ng too large', lng)
    if ( mf%nc < (c+nc-1) ) call bl_error('ZMULTIFAB_FILL_BOUNDARY_C: nc too large', nc)
    if ( lng < 1 ) return
    call build(bpt, "zmf_fill_boundary_c")
    call mf_fb_fancy_z(mf, c, nc, lng, lcross)
    call destroy(bpt)
  end subroutine zmultifab_fill_boundary_c

  subroutine zmultifab_fill_boundary(mf, ng, cross)
    type(zmultifab), intent(inout) :: mf
    integer, intent(in), optional  :: ng
    logical, intent(in), optional  :: cross
    call zmultifab_fill_boundary_c(mf, 1, mf%nc, ng, cross)
  end subroutine zmultifab_fill_boundary

  subroutine multifab_sum_boundary_c(mf, c, nc, ng)
    type(multifab), intent(inout)        :: mf
    integer,        intent(in)           :: c, nc
    integer,        intent(in), optional :: ng

    real(dp_t), pointer     :: p(:,:,:,:), pdst(:,:,:,:), psrc(:,:,:,:)
    integer, allocatable    :: rst(:)
    integer, parameter      :: tag = 1713
    integer                 :: i, ii, jj, np, sh(MAX_SPACEDIM+1), lng
    type(boxassoc)          :: bxasc
    real(dp_t), allocatable :: g_snd_d(:), g_rcv_d(:)

    lng = mf%ng; if ( present(ng) ) lng = ng

    if ( mf%nc < (c+nc-1) ) call bl_error('MULTIFAB_SUM_BOUNDARY_C: nc too large', nc)
    if ( lng > mf%ng      ) call bl_error('MULTIFAB_SUM_BOUNDARY_C: ng too large', lng)
    if ( lng < 1          ) return

    call boxassoc_build(bxasc, mf%la%lap, lng, mf%nodal, .false., .true.)

    do i = 1, bxasc%l_con%ncpy
       ii   =  local_index(mf,bxasc%l_con%cpy(i)%nd)
       jj   =  local_index(mf,bxasc%l_con%cpy(i)%ns)
       pdst => dataptr(mf%fbs(ii), bxasc%l_con%cpy(i)%dbx, c, nc)
       psrc => dataptr(mf%fbs(jj), bxasc%l_con%cpy(i)%sbx, c, nc)
       call sum_d(pdst,psrc)
    end do

    np = parallel_nprocs()

    if (np == 1) return

    allocate(g_snd_d(nc*bxasc%r_con%svol))
    allocate(g_rcv_d(nc*bxasc%r_con%rvol))

    do i = 1, bxasc%r_con%nsnd
       p => dataptr(mf, local_index(mf,bxasc%r_con%snd(i)%ns), bxasc%r_con%snd(i)%sbx, c, nc)
       call reshape_d_4_1(g_snd_d, 1 + nc*bxasc%r_con%snd(i)%pv, p)
    end do

    allocate(rst(bxasc%r_con%nrp))
    do i = 1, bxasc%r_con%nrp
       rst(i) = parallel_irecv_dv(g_rcv_d(1+nc*bxasc%r_con%rtr(i)%pv:), &
            nc*bxasc%r_con%rtr(i)%sz, bxasc%r_con%rtr(i)%pr, tag)
    end do
    do i = 1, bxasc%r_con%nsp
       call parallel_send_dv(g_snd_d(1+nc*bxasc%r_con%str(i)%pv), &
            nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
    end do
    call parallel_wait(rst)

    do i = 1, bxasc%r_con%nrcv
       sh = bxasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mf, local_index(mf,bxasc%r_con%rcv(i)%nd), bxasc%r_con%rcv(i)%dbx, c, nc)
       call reshape_d_1_4(p, g_rcv_d, 1 + nc*bxasc%r_con%rcv(i)%pv, sh, sum_d)
    end do

  end subroutine multifab_sum_boundary_c

  subroutine multifab_sum_boundary(mf, ng)
    type(multifab), intent(inout)        :: mf
    integer,        intent(in), optional :: ng
    call multifab_sum_boundary_c(mf, 1, mf%nc, ng)
  end subroutine multifab_sum_boundary
  !
  ! This does a logical "or" of ghost cells that overlay valid cells with those valid cells.
  !
  subroutine lmultifab_sum_boundary_c(mf, c, nc, ng)
    type(lmultifab), intent(inout)        :: mf
    integer,         intent(in)           :: c, nc
    integer,         intent(in), optional :: ng

    logical, pointer     :: p(:,:,:,:), pdst(:,:,:,:), psrc(:,:,:,:)
    integer, allocatable :: rst(:)
    integer, parameter   :: tag = 1713
    integer              :: i, ii, jj, np, sh(MAX_SPACEDIM+1), lng
    type(boxassoc)       :: bxasc
    logical, allocatable :: g_snd_l(:), g_rcv_l(:)

    lng = mf%ng; if ( present(ng) ) lng = ng

    if ( mf%nc < (c+nc-1) ) call bl_error('LMULTIFAB_SUM_BOUNDARY_C: nc too large', nc)
    if ( lng > mf%ng      ) call bl_error('LMULTIFAB_SUM_BOUNDARY_C: ng too large', lng)
    if ( lng < 1          ) return

    call boxassoc_build(bxasc, mf%la%lap, lng, mf%nodal, .false., .true.)

    do i = 1, bxasc%l_con%ncpy
       ii   =  local_index(mf,bxasc%l_con%cpy(i)%nd)
       jj   =  local_index(mf,bxasc%l_con%cpy(i)%ns)
       pdst => dataptr(mf%fbs(ii), bxasc%l_con%cpy(i)%dbx, c, nc)
       psrc => dataptr(mf%fbs(jj), bxasc%l_con%cpy(i)%sbx, c, nc)
       call logical_or(pdst,psrc)
    end do

    np = parallel_nprocs()
    
    if (np == 1) then
       call boxassoc_destroy(bxasc)
       return
    end if

    allocate(g_snd_l(nc*bxasc%r_con%svol))
    allocate(g_rcv_l(nc*bxasc%r_con%rvol))

    do i = 1, bxasc%r_con%nsnd
       p => dataptr(mf, local_index(mf,bxasc%r_con%snd(i)%ns), bxasc%r_con%snd(i)%sbx, c, nc)
       call reshape_l_4_1(g_snd_l, 1 + nc*bxasc%r_con%snd(i)%pv, p)
    end do

    allocate(rst(bxasc%r_con%nrp))
    do i = 1, bxasc%r_con%nrp
       rst(i) = parallel_irecv_lv(g_rcv_l(1+nc*bxasc%r_con%rtr(i)%pv:), &
            nc*bxasc%r_con%rtr(i)%sz, bxasc%r_con%rtr(i)%pr, tag)
    end do
    do i = 1, bxasc%r_con%nsp
       call parallel_send_lv(g_snd_l(1+nc*bxasc%r_con%str(i)%pv), &
            nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
    end do
    call parallel_wait(rst)

    do i = 1, bxasc%r_con%nrcv
       sh = bxasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mf, local_index(mf,bxasc%r_con%rcv(i)%nd), bxasc%r_con%rcv(i)%dbx, c, nc)
       call reshape_l_1_4(p, g_rcv_l, 1 + nc*bxasc%r_con%rcv(i)%pv, sh, logical_or)
    end do

    call boxassoc_destroy(bxasc)

  end subroutine lmultifab_sum_boundary_c

  subroutine lmultifab_sum_boundary(mf, ng)
    type(lmultifab), intent(inout)        :: mf
    integer,        intent(in), optional  :: ng
    call lmultifab_sum_boundary_c(mf, 1, mf%nc, ng)
  end subroutine lmultifab_sum_boundary

  subroutine mf_internal_sync_fancy(mf, c, nc, lall, filter)
    type(multifab), intent(inout)           :: mf
    integer, intent(in)                     :: c
    integer, intent(in)                     :: nc
    logical, intent(in)                     :: lall

    real(dp_t), dimension(:,:,:,:), pointer :: pdst, psrc, p
    integer                                 :: i, ii, jj, sh(MAX_SPACEDIM+1), np
    integer, parameter                      :: tag = 1104
    type(syncassoc)                         :: snasc
    integer,    allocatable                 :: rst(:)
    real(dp_t), allocatable                 :: g_snd_d(:), g_rcv_d(:)

    interface
       subroutine filter(out, in)
         use bl_types
         real(dp_t), intent(inout) :: out(:,:,:,:)
         real(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface

    optional filter

    snasc = layout_syncassoc(mf%la, mf%ng, mf%nodal, lall)

    do i = 1, snasc%l_con%ncpy
       ii   =  local_index(mf,snasc%l_con%cpy(i)%nd)
       jj   =  local_index(mf,snasc%l_con%cpy(i)%ns)
       pdst => dataptr(mf%fbs(ii), snasc%l_con%cpy(i)%dbx, c, nc)
       psrc => dataptr(mf%fbs(jj), snasc%l_con%cpy(i)%sbx, c, nc)
       call cpy_d(pdst, psrc, filter)
    end do

    np = parallel_nprocs()

    if (np == 1) return

    allocate(g_snd_d(nc*snasc%r_con%svol))
    allocate(g_rcv_d(nc*snasc%r_con%rvol))

    do i = 1, snasc%r_con%nsnd
       p => dataptr(mf, local_index(mf,snasc%r_con%snd(i)%ns), snasc%r_con%snd(i)%sbx, c, nc)
       call reshape_d_4_1(g_snd_d, 1 + nc*snasc%r_con%snd(i)%pv, p)
    end do

    allocate(rst(snasc%r_con%nrp))
    do i = 1, snasc%r_con%nrp
       rst(i) = parallel_irecv_dv(g_rcv_d(1+nc*snasc%r_con%rtr(i)%pv:), &
            nc*snasc%r_con%rtr(i)%sz, snasc%r_con%rtr(i)%pr, tag)
    end do
    do i = 1, snasc%r_con%nsp
       call parallel_send_dv(g_snd_d(1+nc*snasc%r_con%str(i)%pv), &
            nc*snasc%r_con%str(i)%sz, snasc%r_con%str(i)%pr, tag)
    end do
    call parallel_wait(rst)

    do i = 1, snasc%r_con%nrcv
       sh = snasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mf, local_index(mf,snasc%r_con%rcv(i)%nd), snasc%r_con%rcv(i)%dbx, c, nc)
       call reshape_d_1_4(p, g_rcv_d, 1 + nc*snasc%r_con%rcv(i)%pv, sh, filter)
    end do

  end subroutine mf_internal_sync_fancy
  !!
  !! Internal Sync makes sure that any overlapped values are reconciled
  !! by copying values from the lower index number fabs to the higher index
  !! numbered boxes.  Works cell centered and node centered.  Though in a typical
  !! cell-centered multifab, there are no overlaps to reconcile.
  !! If ALL is true then even ghost cell data is 'reconciled'
  !!
  subroutine multifab_internal_sync_c(mf, c, nc, all, filter)
    type(multifab), intent(inout) :: mf
    integer, intent(in)           :: c
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    integer                       :: lnc
    logical                       :: lall
    type(bl_prof_timer), save     :: bpt
    interface
       subroutine filter(out, in)
         use bl_types
         real(dp_t), intent(inout) :: out(:,:,:,:)
         real(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    lnc  = 1;        if ( present(nc)  ) lnc  = nc
    lall = .false. ; if ( present(all) ) lall = all
    if ( mf%nc < (c+lnc-1) ) call bl_error('MULTIFAB_INTERNAL_SYNC_C: nc too large', lnc)
    call build(bpt, "mf_internal_sync")
    call mf_internal_sync_fancy(mf, c, lnc, lall, filter)
    call destroy(bpt)
  end subroutine multifab_internal_sync_c

  subroutine multifab_internal_sync(mf, all, filter)
    type(multifab), intent(inout) :: mf
    logical, intent(in), optional :: all
    interface
       subroutine filter(out, in)
         use bl_types
         real(dp_t), intent(inout) :: out(:,:,:,:)
         real(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    call multifab_internal_sync_c(mf, 1, mf%nc, all, filter)
  end subroutine multifab_internal_sync

  subroutine lmf_internal_sync_fancy(mf, c, nc, lall, filter)
    type(lmultifab), intent(inout)       :: mf
    integer, intent(in)                  :: c
    integer, intent(in)                  :: nc
    logical, intent(in)                  :: lall

    logical, dimension(:,:,:,:), pointer :: pdst, psrc, p
    integer                              :: i, ii, jj, sh(MAX_SPACEDIM+1), np
    integer, parameter                   :: tag = 1104
    type(syncassoc)                      :: snasc
    integer, allocatable                 :: rst(:)
    logical, allocatable                 :: g_snd_l(:), g_rcv_l(:)

    interface
       subroutine filter(out, in)
         use bl_types
         logical, intent(inout) :: out(:,:,:,:)
         logical, intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface

    optional filter

    snasc = layout_syncassoc(mf%la, mf%ng, mf%nodal, lall)

    do i = 1, snasc%l_con%ncpy
       ii   =  local_index(mf,snasc%l_con%cpy(i)%nd)
       jj   =  local_index(mf,snasc%l_con%cpy(i)%ns)
       pdst => dataptr(mf%fbs(ii), snasc%l_con%cpy(i)%dbx, c, nc)
       psrc => dataptr(mf%fbs(jj), snasc%l_con%cpy(i)%sbx, c, nc)
       call cpy_l(pdst, psrc, filter)
    end do

    np = parallel_nprocs()

    if (np == 1) return

    allocate(g_snd_l(nc*snasc%r_con%svol))
    allocate(g_rcv_l(nc*snasc%r_con%rvol))

    do i = 1, snasc%r_con%nsnd
       p => dataptr(mf, local_index(mf,snasc%r_con%snd(i)%ns), snasc%r_con%snd(i)%sbx, c, nc)
       call reshape_l_4_1(g_snd_l, 1 + nc*snasc%r_con%snd(i)%pv, p)
    end do

    allocate(rst(snasc%r_con%nrp))
    do i = 1, snasc%r_con%nrp
       rst(i) = parallel_irecv_lv(g_rcv_l(1+nc*snasc%r_con%rtr(i)%pv:), &
            nc*snasc%r_con%rtr(i)%sz, snasc%r_con%rtr(i)%pr, tag)
    end do
    do i = 1, snasc%r_con%nsp
       call parallel_send_lv(g_snd_l(1+nc*snasc%r_con%str(i)%pv), &
            nc*snasc%r_con%str(i)%sz, snasc%r_con%str(i)%pr, tag)
    end do
    call parallel_wait(rst)

    do i = 1, snasc%r_con%nrcv
       sh = snasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mf, local_index(mf,snasc%r_con%rcv(i)%nd), snasc%r_con%rcv(i)%dbx, c, nc)
       call reshape_l_1_4(p, g_rcv_l, 1 + nc*snasc%r_con%rcv(i)%pv, sh, filter)
    end do

  end subroutine lmf_internal_sync_fancy
  !!
  !! Internal Sync makes sure that any overlapped values are reconciled
  !! by coping values from the lower index number fabs to the higher index
  !! numbered boxes.  Works cell centered and node centered.  Though in a typical
  !! cell-centered multifab, there are no overlaps to reconcile.
  !!
  subroutine lmultifab_internal_sync_c(mf, c, nc, all, filter)
    type(lmultifab), intent(inout) :: mf
    integer, intent(in)            :: c
    integer, intent(in), optional  :: nc
    logical, intent(in), optional  :: all
    integer                        :: lnc
    logical                        :: lall
    type(bl_prof_timer), save      :: bpt
    interface
       subroutine filter(out, in)
         use bl_types
         logical, intent(inout) :: out(:,:,:,:)
         logical, intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    lnc  = 1;        if ( present(nc)  ) lnc  = nc
    lall = .false. ; if ( present(all) ) lall = all
    if ( mf%nc < (c+lnc-1) ) call bl_error('LMULTIFAB_INTERNAL_SYNC_C: nc too large', lnc)
    call build(bpt, "lmf_internal_sync")
    call lmf_internal_sync_fancy(mf, c, lnc, lall, filter)
    call destroy(bpt)
  end subroutine lmultifab_internal_sync_c

  subroutine multifab_copy_on_shift(mf_out, c_out, mf_in, c_in, nc, len, face)
    type(multifab), intent(in   ) :: mf_in
    type(multifab), intent(inout) :: mf_out
    integer, intent(in)           :: c_in, c_out
    integer, intent(in)           :: len,face
    integer, intent(in), optional :: nc

    type(box)                 :: jbx, abx
    real(dp_t), pointer       :: pdst(:,:,:,:), psrc(:,:,:,:)
    integer                   :: i, j, lnc
    integer, parameter        :: tag = 1204
    type(bl_prof_timer), save :: bpt

    lnc = 1; if ( present(nc)  ) lnc  = nc

    if ( mf_in%nc  < (c_in +lnc-1) ) call bl_error('MULTIFAB_COPY_ON_SHIFT: nc too large', lnc)
    if ( mf_out%nc < (c_out+lnc-1) ) call bl_error('MULTIFAB_COPY_ON_SHIFT: nc too large', lnc)

    call build(bpt, "mf_copy_on_shift")

    do j = 1, nboxes(mf_in%la)
       jbx = shift(box_nodalize(get_box(mf_in%la,j), mf_in%nodal), len, face)
       do i = 1, nboxes(mf_out%la)
          if ( remote(mf_in%la,j) .and. remote(mf_out%la,i) ) cycle
          abx = intersection(box_nodalize(get_box(mf_out%la,i), mf_out%nodal), jbx)
          if ( empty(abx) ) cycle
          if ( local(mf_out%la, i) .and. local(mf_in%la, j) ) then
             pdst => dataptr(mf_out, local_index(mf_out,i), abx,                  c_out, lnc)
             psrc => dataptr(mf_in,  local_index(mf_in, j), shift(abx,-len,face), c_in,  lnc)
             call cpy_d(pdst, psrc)
          else if ( local(mf_in%la,j) ) then ! must send
             psrc => dataptr(mf_in, local_index(mf_in,j), shift(abx,-len,face), c_in, lnc)
             call parallel_send(psrc, get_proc(mf_out%la,i), tag)
          else if ( local(mf_out%la,i) ) then  ! must recv
             pdst => dataptr(mf_out, local_index(mf_out,i), abx, c_out, lnc)
             call parallel_recv(pdst, get_proc(mf_in%la,j), tag)
          end if
       end do
    end do

    call destroy(bpt)

  end subroutine multifab_copy_on_shift

  subroutine lmultifab_internal_sync(mf, all, filter)
    type(lmultifab), intent(inout) :: mf
    logical, intent(in), optional  :: all
    interface
       subroutine filter(out, in)
         logical, intent(inout) :: out(:,:,:,:)
         logical, intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    call lmultifab_internal_sync_c(mf, 1, mf%nc, all, filter)
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
    character(len=6) :: fn
    un = unit_stdout(unit)
    call unit_skip(un, skip)

    if ( parallel_IOProcessor() ) then
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
       write(unit=un, fmt='(" NBOXES  = ",i2)') nboxes(mf%la)
    end if

    do ii = 0, parallel_nprocs()
       if ( ii == parallel_myproc() ) then
          do i = 1, nlocal(mf%la)
             write(unit=fn, fmt='(i6)') global_index(mf%la,i)
             call print(mf%fbs(i), str = fn, unit = unit, all = all, data = data, &
                  skip = unit_get_skip(skip) + 2)
          end do
       end if
       call parallel_barrier()
    end do
  end subroutine multifab_print

  subroutine multifab_print_c(mf, comp, str, unit, all, data, skip)
    use bl_IO_module
    type(multifab), intent(in) :: mf
    integer, intent(in) :: comp
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: all, data
    integer, intent(in), optional :: skip
    integer :: i, ii
    integer :: un
    character(len=6) :: fn
    un = unit_stdout(unit)
    call unit_skip(un, skip)

    if ( parallel_IOProcessor() ) then
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
       write(unit=un, fmt='(" NBOXES  = ",i2)') nboxes(mf%la)
    end if

    do ii = 0, parallel_nprocs()
       if ( ii == parallel_myproc() ) then
          do i = 1, nlocal(mf%la)
             write(unit=fn, fmt='(i6)') global_index(mf%la,i)
             call fab_print(mf%fbs(i), comp, str = fn, unit = unit, all = all, data = data, &
                  skip = unit_get_skip(skip) + 2)
          end do
       end if
       call parallel_barrier()
    end do
  end subroutine multifab_print_c

  subroutine imultifab_print(mf, str, unit, all, data, skip)
    use bl_IO_module
    type(imultifab), intent(in) :: mf
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: all, data
    integer, intent(in), optional :: skip
    integer :: i, ii
    integer :: un
    character(len=6) :: fn
    un = unit_stdout(unit)
    call unit_skip(un, skip)

    if ( parallel_IOProcessor() ) then
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
       write(unit=un, fmt='(" NBOXES  = ",i2)') nboxes(mf%la)
    end if

    do ii = 0, parallel_nprocs()
       if ( ii == parallel_myproc() ) then
          do i = 1, nlocal(mf%la)
             write(unit=fn, fmt='(i6)') global_index(mf%la,i)
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
    character(len=6) :: fn
    un = unit_stdout(unit)
    call unit_skip(un, skip)

    if ( parallel_IOProcessor() ) then
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
       write(unit=un, fmt='(" NBOXES  = ",i2)') nboxes(mf%la)
    end if

    do ii = 0, parallel_nprocs()
       if ( ii == parallel_myproc() ) then
          do i = 1, nlocal(mf%la)
             write(unit=fn, fmt='(i6)') global_index(mf%la,i)
             call print(mf%fbs(i), str = fn, unit = unit, all = all, data = data, &
                  skip = unit_get_skip(skip) + 2)
          end do
       end if
       call parallel_barrier()
    end do
  end subroutine lmultifab_print
  subroutine zmultifab_print(mf, str, unit, all, data, skip)
    use bl_IO_module
    type(zmultifab), intent(in) :: mf
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: all, data
    integer, intent(in), optional :: skip
    integer :: i, ii
    integer :: un
    character(len=6) :: fn
    un = unit_stdout(unit)
    call unit_skip(un, skip)

    if ( parallel_IOProcessor() ) then
       write(unit=un, fmt='("ZMULTIFAB")', advance = 'no')
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
       write(unit=un, fmt='(" NBOXES  = ",i2)') nboxes(mf%la)
    end if

    do ii = 0, parallel_nprocs()
       if ( ii == parallel_myproc() ) then
          do i = 1, nlocal(mf%la)
             write(unit=fn, fmt='(i6)') global_index(mf%la,i)
             call print(mf%fbs(i), str = fn, unit = unit, all = all, data = data, &
                  skip = unit_get_skip(skip) + 2)
          end do
       end if
       call parallel_barrier()
    end do
  end subroutine zmultifab_print

  subroutine mf_copy_fancy_double(mdst, dstcomp, msrc, srccomp, nc, filter)
    type(multifab), intent(inout) :: mdst
    type(multifab), intent(in)    :: msrc
    integer, intent(in)           :: dstcomp, srccomp, nc

    interface
       subroutine filter(out, in)
         use bl_types
         real(dp_t), intent(inout) :: out(:,:,:,:)
         real(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface

    optional filter

    type(copyassoc)         :: cpasc
    real(dp_t), pointer     :: p(:,:,:,:), pdst(:,:,:,:), psrc(:,:,:,:)
    integer, allocatable    :: rst(:)
    integer, parameter      :: tag = 1102
    integer                 :: i, ii, jj, sh(MAX_SPACEDIM+1), np
    real(dp_t), allocatable :: g_snd_d(:), g_rcv_d(:)

    cpasc = layout_copyassoc(mdst%la, msrc%la, mdst%nodal, msrc%nodal)

    do i = 1, cpasc%l_con%ncpy
       ii   =  local_index(mdst,cpasc%l_con%cpy(i)%nd)
       jj   =  local_index(msrc,cpasc%l_con%cpy(i)%ns)
       pdst => dataptr(mdst%fbs(ii), cpasc%l_con%cpy(i)%dbx, dstcomp, nc)
       psrc => dataptr(msrc%fbs(jj), cpasc%l_con%cpy(i)%sbx, srccomp, nc)
       call cpy_d(pdst, psrc, filter)
    end do

    np = parallel_nprocs()

    if (np == 1) return

    allocate(g_snd_d(nc*cpasc%r_con%svol))
    allocate(g_rcv_d(nc*cpasc%r_con%rvol))

    do i = 1, cpasc%r_con%nsnd
       p => dataptr(msrc, local_index(msrc,cpasc%r_con%snd(i)%ns), cpasc%r_con%snd(i)%sbx, srccomp, nc)
       call reshape_d_4_1(g_snd_d, 1 + nc*cpasc%r_con%snd(i)%pv, p)
    end do

    allocate(rst(cpasc%r_con%nrp))
    do i = 1, cpasc%r_con%nrp
       rst(i) = parallel_irecv_dv(g_rcv_d(1+nc*cpasc%r_con%rtr(i)%pv:), &
            nc*cpasc%r_con%rtr(i)%sz, cpasc%r_con%rtr(i)%pr, tag)
    end do
    do i = 1, cpasc%r_con%nsp
       call parallel_send_dv(g_snd_d(1+nc*cpasc%r_con%str(i)%pv), &
            nc*cpasc%r_con%str(i)%sz, cpasc%r_con%str(i)%pr, tag)
    end do
    call parallel_wait(rst)

    do i = 1, cpasc%r_con%nrcv
       sh = cpasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mdst, local_index(mdst,cpasc%r_con%rcv(i)%nd), cpasc%r_con%rcv(i)%dbx, dstcomp, nc)
       call reshape_d_1_4(p, g_rcv_d, 1 + nc*cpasc%r_con%rcv(i)%pv, sh, filter)
    end do

  end subroutine mf_copy_fancy_double

  subroutine mf_copy_fancy_integer(mdst, dstcomp, msrc, srccomp, nc, filter)
    type(imultifab), intent(inout) :: mdst
    type(imultifab), intent(in)    :: msrc
    integer, intent(in)            :: dstcomp, srccomp, nc

    interface
       subroutine filter(out, in)
         use bl_types
         integer, intent(inout) :: out(:,:,:,:)
         integer, intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface

    optional filter

    type(copyassoc)      :: cpasc
    integer, pointer     :: p(:,:,:,:), pdst(:,:,:,:), psrc(:,:,:,:)
    integer, allocatable :: rst(:)
    integer, parameter   :: tag = 1102
    integer              :: i, ii, jj, np, sh(MAX_SPACEDIM+1)
    integer, allocatable :: g_snd_i(:), g_rcv_i(:)

    cpasc = layout_copyassoc(mdst%la, msrc%la, mdst%nodal, msrc%nodal)

    do i = 1, cpasc%l_con%ncpy
       ii   =  local_index(mdst,cpasc%l_con%cpy(i)%nd)
       jj   =  local_index(msrc,cpasc%l_con%cpy(i)%ns)
       pdst => dataptr(mdst%fbs(ii), cpasc%l_con%cpy(i)%dbx, dstcomp, nc)
       psrc => dataptr(msrc%fbs(jj), cpasc%l_con%cpy(i)%sbx, srccomp, nc)
       call cpy_i(pdst, psrc, filter)
    end do

    np = parallel_nprocs()

    if (np == 1) return

    allocate(g_snd_i(nc*cpasc%r_con%svol))
    allocate(g_rcv_i(nc*cpasc%r_con%rvol))

    do i = 1, cpasc%r_con%nsnd
       p => dataptr(msrc, local_index(msrc,cpasc%r_con%snd(i)%ns), cpasc%r_con%snd(i)%sbx, srccomp, nc)
       call reshape_i_4_1(g_snd_i, 1 + nc*cpasc%r_con%snd(i)%pv, p)
    end do

    allocate(rst(cpasc%r_con%nrp))
    do i = 1, cpasc%r_con%nrp
       rst(i) = parallel_irecv_iv(g_rcv_i(1+nc*cpasc%r_con%rtr(i)%pv:), &
            nc*cpasc%r_con%rtr(i)%sz, cpasc%r_con%rtr(i)%pr, tag)
    end do
    do i = 1, cpasc%r_con%nsp
       call parallel_send_iv(g_snd_i(1+nc*cpasc%r_con%str(i)%pv), &
            nc*cpasc%r_con%str(i)%sz, cpasc%r_con%str(i)%pr, tag)
    end do
    call parallel_wait(rst)

    do i = 1, cpasc%r_con%nrcv
       sh = cpasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mdst, local_index(mdst,cpasc%r_con%rcv(i)%nd), cpasc%r_con%rcv(i)%dbx, dstcomp, nc)
       call reshape_i_1_4(p, g_rcv_i, 1 + nc*cpasc%r_con%rcv(i)%pv, sh, filter)
    end do

  end subroutine mf_copy_fancy_integer

  subroutine mf_copy_fancy_logical(mdst, dstcomp, msrc, srccomp, nc, filter)
    type(lmultifab), intent(inout) :: mdst
    type(lmultifab), intent(in)    :: msrc
    integer, intent(in)            :: dstcomp, srccomp, nc

    interface
       subroutine filter(out, in)
         use bl_types
         logical, intent(inout) :: out(:,:,:,:)
         logical, intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface

    optional filter

    type(copyassoc)      :: cpasc
    logical, pointer     :: p(:,:,:,:), pdst(:,:,:,:), psrc(:,:,:,:)
    integer, allocatable :: rst(:)
    integer, parameter   :: tag = 1102
    integer              :: i, ii, jj, np, sh(MAX_SPACEDIM+1)
    logical, allocatable :: g_snd_l(:), g_rcv_l(:)

    cpasc = layout_copyassoc(mdst%la, msrc%la, mdst%nodal, msrc%nodal)

    do i = 1, cpasc%l_con%ncpy
       ii   =  local_index(mdst,cpasc%l_con%cpy(i)%nd)
       jj   =  local_index(msrc,cpasc%l_con%cpy(i)%ns)
       pdst => dataptr(mdst%fbs(ii), cpasc%l_con%cpy(i)%dbx, dstcomp, nc)
       psrc => dataptr(msrc%fbs(jj), cpasc%l_con%cpy(i)%sbx, srccomp, nc)
       call cpy_l(pdst, psrc, filter)
    end do

    np = parallel_nprocs()

    if (np == 1) return

    allocate(g_snd_l(nc*cpasc%r_con%svol))
    allocate(g_rcv_l(nc*cpasc%r_con%rvol))

    do i = 1, cpasc%r_con%nsnd
       p => dataptr(msrc, local_index(msrc,cpasc%r_con%snd(i)%ns), cpasc%r_con%snd(i)%sbx, srccomp, nc)
       call reshape_l_4_1(g_snd_l, 1 + nc*cpasc%r_con%snd(i)%pv, p)
    end do

    allocate(rst(cpasc%r_con%nrp))
    do i = 1, cpasc%r_con%nrp
       rst(i) = parallel_irecv_lv(g_rcv_l(1+nc*cpasc%r_con%rtr(i)%pv:), &
            nc*cpasc%r_con%rtr(i)%sz, cpasc%r_con%rtr(i)%pr, tag)
    end do
    do i = 1, cpasc%r_con%nsp
       call parallel_send_lv(g_snd_l(1+nc*cpasc%r_con%str(i)%pv), &
            nc*cpasc%r_con%str(i)%sz, cpasc%r_con%str(i)%pr, tag)
    end do
    call parallel_wait(rst)

    do i = 1, cpasc%r_con%nrcv
       sh = cpasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mdst, local_index(mdst,cpasc%r_con%rcv(i)%nd), cpasc%r_con%rcv(i)%dbx, dstcomp, nc)
       call reshape_l_1_4(p, g_rcv_l, 1 + nc*cpasc%r_con%rcv(i)%pv, sh, filter)
    end do

  end subroutine mf_copy_fancy_logical

  subroutine mf_copy_fancy_z(mdst, dstcomp, msrc, srccomp, nc, filter)
    type(zmultifab), intent(inout) :: mdst
    type(zmultifab), intent(in)    :: msrc
    integer, intent(in)            :: dstcomp, srccomp, nc

    interface
       subroutine filter(out, in)
         use bl_types
         complex(dp_t), intent(inout) :: out(:,:,:,:)
         complex(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface

    optional filter

    type(copyassoc)            :: cpasc
    complex(dp_t), pointer     :: p(:,:,:,:), pdst(:,:,:,:), psrc(:,:,:,:)
    integer, allocatable       :: rst(:)
    integer, parameter         :: tag = 1102
    integer                    :: i, ii, jj, sh(MAX_SPACEDIM+1)
    complex(dp_t), allocatable :: g_snd_z(:), g_rcv_z(:)

    cpasc = layout_copyassoc(mdst%la, msrc%la, mdst%nodal, msrc%nodal)

    do i = 1, cpasc%l_con%ncpy
       ii   =  local_index(mdst,cpasc%l_con%cpy(i)%nd)
       jj   =  local_index(msrc,cpasc%l_con%cpy(i)%ns)
       pdst => dataptr(mdst%fbs(ii), cpasc%l_con%cpy(i)%dbx, dstcomp, nc)
       psrc => dataptr(msrc%fbs(jj), cpasc%l_con%cpy(i)%sbx, srccomp, nc)
       call cpy_z(pdst, psrc, filter)
    end do

    if (parallel_nprocs() == 1) return

    allocate(g_snd_z(nc*cpasc%r_con%svol))
    allocate(g_rcv_z(nc*cpasc%r_con%rvol))

    do i = 1, cpasc%r_con%nsnd
       p => dataptr(msrc, local_index(msrc,cpasc%r_con%snd(i)%ns), cpasc%r_con%snd(i)%sbx, srccomp, nc)
       call reshape_z_4_1(g_snd_z, 1 + nc*cpasc%r_con%snd(i)%pv, p)
    end do

    allocate(rst(cpasc%r_con%nrp))
    do i = 1, cpasc%r_con%nrp
       rst(i) = parallel_irecv_zv(g_rcv_z(1+nc*cpasc%r_con%rtr(i)%pv:), &
            nc*cpasc%r_con%rtr(i)%sz, cpasc%r_con%rtr(i)%pr, tag)
    end do
    do i = 1, cpasc%r_con%nsp
       call parallel_send_zv(g_snd_z(1+nc*cpasc%r_con%str(i)%pv), &
            nc*cpasc%r_con%str(i)%sz, cpasc%r_con%str(i)%pr, tag)
    end do
    call parallel_wait(rst)

    do i = 1, cpasc%r_con%nrcv
       sh = cpasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mdst, local_index(mdst,cpasc%r_con%rcv(i)%nd), cpasc%r_con%rcv(i)%dbx, dstcomp, nc)
       call reshape_z_1_4(p, g_rcv_z, 1 + nc*cpasc%r_con%rcv(i)%pv, sh, filter)
    end do

  end subroutine mf_copy_fancy_z

  subroutine multifab_copy_c(mdst, dstcomp, msrc, srccomp, nc, ng, filter)
    type(multifab), intent(inout) :: mdst
    type(multifab), intent(in)    :: msrc
    integer, intent(in)           :: dstcomp, srccomp
    integer, intent(in), optional :: nc
    integer, intent(in), optional :: ng
    real(dp_t), pointer           :: pdst(:,:,:,:), psrc(:,:,:,:)
    integer                       :: i, lnc, lng
    interface
       subroutine filter(out, in)
         use bl_types
         real(dp_t), intent(inout) :: out(:,:,:,:)
         real(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    type(bl_prof_timer), save :: bpt
    call build(bpt, "mf_copy_c")
    lnc     = 1;       if ( present(nc)     ) lnc = nc
    lng     = 0;       if ( present(ng)     ) lng = ng
    if ( lnc < 1 )                   call bl_error('MULTIFAB_COPY_C: nc must be >= 1')
    if ( mdst%nc < (dstcomp+lnc-1) ) call bl_error('MULTIFAB_COPY_C: nc too large for dst multifab', lnc)
    if ( msrc%nc < (srccomp+lnc-1) ) call bl_error('MULTIFAB_COPY_C: nc too large for src multifab', lnc)
    if ( lng > 0 )                   call bl_assert(mdst%ng >= ng, msrc%ng >= ng,"not enough ghost cells in multifab_copy_c")
    if ( mdst%la == msrc%la ) then
       do i = 1, nlocal(mdst%la)
          if ( lng > 0 ) then
             pdst => dataptr(mdst, i, grow(get_ibox(mdst, i),lng), dstcomp, lnc)
             psrc => dataptr(msrc, i, grow(get_ibox(msrc, i),lng), srccomp, lnc)
          else
             pdst => dataptr(mdst, i, get_ibox(mdst, i), dstcomp, lnc)
             psrc => dataptr(msrc, i, get_ibox(msrc, i), srccomp, lnc)
          end if
          call cpy_d(pdst, psrc, filter)
       end do
    else
       if ( lng > 0 ) call bl_error('MULTIFAB_COPY_C: copying ghostcells allowed only when layouts are the same')
       call mf_copy_fancy_double(mdst, dstcomp, msrc, srccomp, lnc, filter)
    end if
    call destroy(bpt)
  end subroutine multifab_copy_c

  subroutine multifab_copy(mdst, msrc, ng, filter)
    type(multifab), intent(inout) :: mdst
    type(multifab), intent(in)    :: msrc
    integer, intent(in), optional :: ng
    interface
       subroutine filter(out, in)
         use bl_types
         real(dp_t), intent(inout) :: out(:,:,:,:)
         real(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    if ( mdst%nc .ne. msrc%nc ) call bl_error('MULTIFAB_COPY: multifabs must have same number of components')
    call multifab_copy_c(mdst, 1, msrc, 1, mdst%nc, ng, filter)
  end subroutine multifab_copy

  subroutine imultifab_copy_c(mdst, dstcomp, msrc, srccomp, nc, ng, filter)
    type(imultifab), intent(inout) :: mdst
    type(imultifab), intent(in)    :: msrc
    integer, intent(in)            :: dstcomp, srccomp
    integer, intent(in), optional  :: nc
    integer, intent(in), optional  :: ng
    integer, pointer               :: pdst(:,:,:,:), psrc(:,:,:,:)
    integer                        :: i, lnc, lng
    interface
       subroutine filter(out, in)
         use bl_types
         integer, intent(inout) :: out(:,:,:,:)
         integer, intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    type(bl_prof_timer), save :: bpt
    call build(bpt, "imf_copy_c")
    lnc     = 1;       if ( present(nc)     ) lnc  = nc
    lng     = 0;       if ( present(ng)     ) lng = ng
    if ( lnc < 1 )                   call bl_error('IMULTIFAB_COPY_C: nc must be >= 1')
    if ( mdst%nc < (dstcomp+lnc-1) ) call bl_error('IMULTIFAB_COPY_C: nc too large for dst multifab', lnc)
    if ( msrc%nc < (srccomp+lnc-1) ) call bl_error('IMULTIFAB_COPY_C: nc too large for src multifab', lnc)
    if ( lng > 0 )                   call bl_assert(mdst%ng >= ng, msrc%ng >= ng,"not enough ghost cells in imultifab_copy_c")
    if ( mdst%la == msrc%la ) then
       do i = 1, nlocal(mdst%la)
          if ( lng > 0 ) then
             pdst => dataptr(mdst, i, grow(get_ibox(mdst, i),lng), dstcomp, lnc)
             psrc => dataptr(msrc, i, grow(get_ibox(msrc, i),lng), srccomp, lnc)
          else
             pdst => dataptr(mdst, i, get_ibox(mdst, i), dstcomp, lnc)
             psrc => dataptr(msrc, i, get_ibox(msrc, i), srccomp, lnc)
          end if
          call cpy_i(pdst, psrc, filter)
       end do
    else
       if ( lng > 0 ) call bl_error('IMULTIFAB_COPY_C: copying ghostcells allowed only when layouts are the same')
       call mf_copy_fancy_integer(mdst, dstcomp, msrc, srccomp, lnc, filter)
    end if
    call destroy(bpt)
  end subroutine imultifab_copy_c

  subroutine imultifab_copy(mdst, msrc, ng, filter)
    type(imultifab), intent(inout) :: mdst
    type(imultifab), intent(in)    :: msrc
    integer, intent(in), optional  :: ng 
    interface
       subroutine filter(out, in)
         use bl_types
         integer, intent(inout) :: out(:,:,:,:)
         integer, intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    if ( mdst%nc .ne. msrc%nc ) call bl_error('IMULTIFAB_COPY: multifabs must have same number of components')
    call imultifab_copy_c(mdst, 1, msrc, 1, mdst%nc, ng, filter)
  end subroutine imultifab_copy

  subroutine lmultifab_copy_c(mdst, dstcomp, msrc, srccomp, nc, ng, filter)
    type(lmultifab), intent(inout) :: mdst
    type(lmultifab), intent(in)    :: msrc
    integer, intent(in)            :: dstcomp, srccomp
    integer, intent(in), optional  :: nc
    integer, intent(in), optional  :: ng
    logical, pointer               :: pdst(:,:,:,:), psrc(:,:,:,:)
    integer                        :: i, lnc, lng
    interface
       subroutine filter(out, in)
         use bl_types
         logical, intent(inout) :: out(:,:,:,:)
         logical, intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    type(bl_prof_timer), save :: bpt
    call build(bpt, "lmf_copy_c")
    lnc     = 1;       if ( present(nc)     ) lnc = nc
    lng     = 0;       if ( present(ng)     ) lng = ng
    if ( lnc < 1 )                   call bl_error('LMULTIFAB_COPY_C: nc must be >= 1')
    if ( mdst%nc < (dstcomp+lnc-1) ) call bl_error('LMULTIFAB_COPY_C: nc too large for dst multifab', lnc)
    if ( msrc%nc < (srccomp+lnc-1) ) call bl_error('LMULTIFAB_COPY_C: nc too large for src multifab', lnc)
    if ( lng > 0 )                   call bl_assert(mdst%ng >= ng, msrc%ng >= ng,"not enough ghost cells in lmultifab_copy_c")
    if ( mdst%la == msrc%la ) then
       do i = 1, nlocal(mdst%la)
          if ( lng > 0 ) then
             pdst => dataptr(mdst, i, grow(get_ibox(mdst, i),lng), dstcomp, lnc)
             psrc => dataptr(msrc, i, grow(get_ibox(msrc, i),lng), srccomp, lnc)
          else
             pdst => dataptr(mdst, i, get_ibox(mdst, i), dstcomp, lnc)
             psrc => dataptr(msrc, i, get_ibox(msrc, i), srccomp, lnc)
          end if
          call cpy_l(pdst, psrc, filter)
       end do
    else
       if ( lng > 0 ) call bl_error('LMULTIFAB_COPY_C: copying ghostcells allowed only when layouts are the same')
       call mf_copy_fancy_logical(mdst, dstcomp, msrc, srccomp, lnc, filter)
    end if
    call destroy(bpt)
  end subroutine lmultifab_copy_c

  subroutine lmultifab_copy(mdst, msrc, ng, filter)
    type(lmultifab), intent(inout) :: mdst
    type(lmultifab), intent(in)    :: msrc
    integer, intent(in), optional  :: ng
    interface
       subroutine filter(out, in)
         use bl_types
         logical, intent(inout) :: out(:,:,:,:)
         logical, intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    if ( mdst%nc .ne. msrc%nc ) call bl_error('LMULTIFAB_COPY: multifabs must have same number of components')
    call lmultifab_copy_c(mdst, 1, msrc, 1, mdst%nc, ng, filter)
  end subroutine lmultifab_copy

  subroutine zmultifab_copy_c(mdst, dstcomp, msrc, srccomp, nc, ng, filter)
    type(zmultifab), intent(inout) :: mdst
    type(zmultifab), intent(in)    :: msrc
    integer, intent(in)            :: dstcomp, srccomp
    integer, intent(in), optional  :: nc
    integer, intent(in), optional  :: ng
    complex(dp_t), pointer         :: pdst(:,:,:,:), psrc(:,:,:,:)
    integer                        :: i, lnc, lng
    interface
       subroutine filter(out, in)
         use bl_types
         complex(dp_t), intent(inout) :: out(:,:,:,:)
         complex(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    lnc     = 1;       if ( present(nc)     ) lnc = nc
    lng     = 0;       if ( present(ng)     ) lng = ng
    if ( lnc < 1 )                   call bl_error('ZMULTIFAB_COPY_C: nc must be >= 1')
    if ( mdst%nc < (dstcomp+lnc-1) ) call bl_error('ZMULTIFAB_COPY_C: nc too large for dst multifab', lnc)
    if ( msrc%nc < (srccomp+lnc-1) ) call bl_error('ZMULTIFAB_COPY_C: nc too large for src multifab', lnc)
    if ( lng > 0 )                   call bl_assert(mdst%ng >= ng, msrc%ng >= ng,"not enough ghost cells in zmultifab_copy_c")
    if ( mdst%la == msrc%la ) then
       do i = 1, nlocal(mdst%la)
          if ( lng > 0 ) then
             pdst => dataptr(mdst, i, grow(get_ibox(mdst, i),lng), dstcomp, lnc)
             psrc => dataptr(msrc, i, grow(get_ibox(msrc, i),lng), srccomp, lnc)
          else
             pdst => dataptr(mdst, i, get_ibox(mdst, i), dstcomp, lnc)
             psrc => dataptr(msrc, i, get_ibox(msrc, i), srccomp, lnc)
          end if
          call cpy_z(pdst, psrc, filter)
       end do
    else
       if ( lng > 0 ) call bl_error('ZMULTIFAB_COPY_C: copying ghostcells allowed only when layouts are the same')
       call mf_copy_fancy_z(mdst, dstcomp, msrc, srccomp, lnc, filter)
    end if
  end subroutine zmultifab_copy_c

  subroutine zmultifab_copy(mdst, msrc, ng, filter)
    type(zmultifab), intent(inout) :: mdst
    type(zmultifab), intent(in)    :: msrc
    integer, intent(in), optional  :: ng
    interface
       subroutine filter(out, in)
         use bl_types
         complex(dp_t), intent(inout) :: out(:,:,:,:)
         complex(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    if ( mdst%nc .ne. msrc%nc ) call bl_error('ZMULTIFAB_COPY: multifabs must have same number of components')
    call zmultifab_copy_c(mdst, 1, msrc, 1, mdst%nc, ng, filter)
  end subroutine zmultifab_copy

  subroutine mf_build_nodal_dot_mask(mask, mf)
    type(multifab), intent(in)  :: mf
    type(multifab), intent(out) :: mask
    integer :: i, d
    type(box) :: full_box, shrunk_box, inner_box
    type(bl_prof_timer), save :: bpt

    call build(bpt, "build_nodal_dot_mask")

    if (.not.nodal_q(mf)) then
       call bl_error("In mf_build_nodal_dot_mask with mf not nodal!")
    end if

    call build(mask, mf%la, 1, 0, mf%nodal)

    !$OMP PARALLEL DO PRIVATE(i,d,full_box,shrunk_box,inner_box)
    do i = 1, nlocal(mf%la)

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
    !$OMP END PARALLEL DO

    call destroy(bpt)

  end subroutine mf_build_nodal_dot_mask

  function multifab_dot_cc(mf, comp, mf1, comp1, mask, nodal_mask, local, comm) result(r)
    real(dp_t)                 :: r
    type(multifab), intent(in) :: mf
    type(multifab), intent(in) :: mf1
    integer,        intent(in) :: comp, comp1

    type(lmultifab), intent(in), optional :: mask
    type(multifab),  intent(in), optional :: nodal_mask
    logical,         intent(in), optional :: local
    integer,         intent(in), optional :: comm

    type(multifab)      :: tmask
    real(dp_t), pointer :: mp(:,:,:,:), mp1(:,:,:,:), ma(:,:,:,:)
    logical,    pointer :: lmp(:,:,:,:)
    real(dp_t)          :: r1,r2
    integer             :: i,j,k,n,lo(4),hi(4)
    logical             :: llocal

    if ( present(mask) ) then
       if ( ncomp(mask) /= 1 ) call bl_error('Mask array is multicomponent')
    end if

    llocal = .false.; if ( present(local) ) llocal = local

    r1 = 0.0_dp_t

    if ( cell_centered_q(mf) ) then

       do n = 1, nlocal(mf%la)
          mp  => dataptr(mf,  n, get_ibox(mf, n), comp)
          mp1 => dataptr(mf1, n, get_ibox(mf1,n), comp1)

          lo = lbound(mp); hi = ubound(mp)

          r2 = 0.0_dp_t

          if ( present(mask) )then
             lmp => dataptr(mask, n, get_ibox(mask, n), 1)
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      if ( lmp(i,j,k,1) ) r2 = r2 + mp(i,j,k,1)*mp1(i,j,k,1)
                   end do
                end do
             end do
          else
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      r2 = r2 + mp(i,j,k,1)*mp1(i,j,k,1)
                   end do
                end do
             end do
          endif

          r1 = r1 + r2
       end do

    else if ( nodal_q(mf) ) then

       if ( .not. present(nodal_mask) ) call build_nodal_dot_mask(tmask, mf)

       do n = 1, nlocal(mf%la) 
          mp  => dataptr(mf,  n, get_ibox(mf,  n), comp)
          mp1 => dataptr(mf1, n, get_ibox(mf1, n), comp1)
          if ( present(nodal_mask) ) then
             ma => dataptr(nodal_mask, n, get_ibox(nodal_mask, n))
          else
             ma => dataptr(tmask,      n, get_ibox(tmask,      n))
          endif

          lo = lbound(mp); hi = ubound(mp)

          r2 = 0.0_dp_t

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   r2 = r2 + ma(i,j,k,1)*mp(i,j,k,1)*mp1(i,j,k,1)
                end do
             end do
          end do

          r1 = r1 + r2
       end do

       if ( .not. present(nodal_mask) ) call destroy(tmask)
    else
       call bl_error("MULTIFAB_DOT_CC, fails when not nodal or cell-centered, can be fixed")
    end if

    r = r1

    if ( .not. llocal ) then
       if ( present(comm) ) then
          call parallel_reduce(r, r1, MPI_SUM, comm = comm)
       else
          call parallel_reduce(r, r1, MPI_SUM)
       end if
    end if

  end function multifab_dot_cc

  function multifab_dot_c(mf, mf1, comp, nodal_mask, comm) result(r)
    real(dp_t)                 :: r
    type(multifab), intent(in) :: mf
    type(multifab), intent(in) :: mf1
    integer       , intent(in) :: comp
    type(multifab), intent(in), optional :: nodal_mask
    integer       , intent(in), optional :: comm
    if ( present(comm) ) then
       r = multifab_dot_cc(mf, comp, mf1, comp, nodal_mask = nodal_mask, local = .false., comm = comm);
    else
       r = multifab_dot_cc(mf, comp, mf1, comp, nodal_mask = nodal_mask, local = .false.);
    end if
  end function multifab_dot_c

  function multifab_dot(mf, mf1, nodal_mask, local, comm) result(r)
    real(dp_t) :: r
    type(multifab), intent(in) :: mf
    type(multifab), intent(in) :: mf1
    type(multifab), intent(in), optional :: nodal_mask
    logical, intent(in), optional :: local
    integer, intent(in), optional :: comm
    if ( present(comm) ) then
       r = multifab_dot_cc(mf, 1, mf1, 1, nodal_mask = nodal_mask, local = local, comm = comm);
    else
       r = multifab_dot_cc(mf, 1, mf1, 1, nodal_mask = nodal_mask, local = local);
    end if
  end function multifab_dot

  subroutine multifab_rescale_2(mf, c, min, max, xmin, xmax, clip)
    type(multifab), intent(inout) :: mf
    real(dp_t), intent(in), optional :: min, max
    real(dp_t), intent(in), optional :: xmin, xmax
    logical, intent(in), optional :: clip
    integer, intent(in) :: c
    integer :: i
    real(dp_t), pointer :: mp(:,:,:,:)
    real(dp_t) :: lmin, lmax, lxmin, lxmax
    logical :: lclip
    lclip = .false. ; if ( present(clip) ) lclip = clip
    if ( present(min) ) then
       lmin = min
    else
       lmin = min_val(mf, c)
    end if
    if ( present(max) ) then
       lmax = max
    else
       lmax = max_val(mf, c)
    end if
    lxmin = 0.0_dp_t ; if ( present(xmin) ) lxmin = xmin
    lxmax = 1.0_dp_t ; if ( present(xmax) ) lxmax = xmax
    if ( lclip ) then
       do i = 1, nlocal(mf%la)
          mp => dataptr(mf, i, get_ibox(mf, i), c)
          where ( mp < lmin )
             mp = lxmin
          elsewhere ( mp > lmax )
             mp = lxmax
          elsewhere
             mp = lxmax*(mp-lmin)/(lmax-lmin) + lxmin
          end where
       end do
    else
       do i = 1, nlocal(mf%la)
          mp => dataptr(mf, i, get_ibox(mf, i), c)
          mp = lxmax*(mp-lmin)/(lmax-lmin) + lxmin
       end do
    end if
  end subroutine multifab_rescale_2

  subroutine multifab_rescale_c(mf, c, val, off)
    real(dp_t), intent(in) :: val
    integer, intent(in) :: c
    real(dp_t), intent(in), optional :: off
    type(multifab), intent(inout) :: mf
    real(dp_t), pointer :: mp(:,:,:,:)
    integer :: i
    do i = 1, nlocal(mf%la)
       mp => dataptr(mf, i, get_ibox(mf, i), c)
       if ( present(off) ) then
          mp = mp*val + off
       else
          mp = mp*val
       end if
    end do
  end subroutine multifab_rescale_c
  subroutine multifab_rescale(mf, val, off)
    real(dp_t), intent(in) :: val
    real(dp_t), intent(in), optional :: off
    type(multifab), intent(inout) :: mf
    real(dp_t), pointer :: mp(:,:,:,:)
    integer :: i
    do i = 1, nlocal(mf%la)
       mp => dataptr(mf, i, get_ibox(mf, i))
       if ( present(off) ) then
          mp = mp*val + off
       else
          mp = mp*val
       end if
    end do
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
    do i = 1, nlocal(a%la)
       ap => dataptr(a, i, get_ibox(a, i))
       bp => dataptr(b, i, get_ibox(b, i))
       cp => dataptr(c, i, get_ibox(c, i))
       ap = b1*bp + c1*cp
    end do
    !$OMP END PARALLEL DO
  end subroutine multifab_saxpy_5

  subroutine multifab_saxpy_4(a, b, c1, c)
    real(dp_t),     intent(in)    :: c1
    type(multifab), intent(inout) :: a
    type(multifab), intent(in   ) :: b,c
    real(dp_t), pointer           :: ap(:,:,:,:)
    real(dp_t), pointer           :: bp(:,:,:,:)
    real(dp_t), pointer           :: cp(:,:,:,:)

    integer :: ii, i, j, k, n, lo(4), hi(4)

    do ii = 1, nlocal(a%la)
       ap => dataptr(a, ii, get_ibox(a,ii))
       bp => dataptr(b, ii, get_ibox(b,ii))
       cp => dataptr(c, ii, get_ibox(c,ii))

       lo = lbound(ap)
       hi = ubound(ap)
       
       ! ap = bp + c1*cp

       !$OMP PARALLEL PRIVATE(i,j,k,n) IF((hi(3)-lo(3)).ge.7)
       do n = lo(4), hi(4)
          !$OMP DO
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   ap(i,j,k,n) = bp(i,j,k,n) + c1 * cp(i,j,k,n)
                end do
             end do
          end do
          !$OMP END DO NOWAIT
       end do
       !$OMP END PARALLEL

    end do
  end subroutine multifab_saxpy_4

  subroutine multifab_saxpy_3_doit(ap, b1, bp)

    real(dp_t),     intent(in   ) :: b1
    real(dp_t), pointer           :: ap(:,:,:,:)
    real(dp_t), pointer           :: bp(:,:,:,:)
    integer                       :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! ap = ap + b1*bp

    !$OMP PARALLEL PRIVATE(i,j,k,n) IF((hi(3)-lo(3)).ge.7)
    do n = lo(4), hi(4)
       !$OMP DO
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) + b1 * bp(i,j,k,n)
             end do
          end do
       end do
       !$OMP END DO NOWAIT
    end do
    !$OMP END PARALLEL

  end subroutine multifab_saxpy_3_doit

  subroutine multifab_saxpy_3(a, b1, b, all)

    real(dp_t),     intent(in   ) :: b1
    type(multifab), intent(inout) :: a
    type(multifab), intent(in   ) :: b
    real(dp_t), pointer           :: ap(:,:,:,:)
    real(dp_t), pointer           :: bp(:,:,:,:)
    logical, intent(in), optional :: all
    integer :: ii
    logical :: lall

    lall = .false.; if ( present(all) ) lall = all

    do ii = 1, nlocal(a%la)
       if ( lall ) then
          ap => dataptr(a,ii)
          bp => dataptr(b,ii)
       else
          ap => dataptr(a, ii, get_ibox(a,ii))
          bp => dataptr(b, ii, get_ibox(b,ii))
       end if
       call multifab_saxpy_3_doit(ap,b1,bp)
    end do

  end subroutine multifab_saxpy_3

  subroutine multifab_saxpy_3_c(a, ia, b1, b, all)
    real(dp_t), intent(in) :: b1
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    integer, intent(in) :: ia
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    logical, intent(in), optional :: all
    integer :: ii
    logical :: lall

    lall = .false.; if ( present(all) ) lall = all

    do ii = 1, nlocal(a%la)
       if ( lall ) then
          ap => dataptr(a,ii,ia)
          bp => dataptr(b,ii)
       else
          ap => dataptr(a, ii, get_ibox(a,ii), ia)
          bp => dataptr(b, ii, get_ibox(b,ii))
       end if
       call multifab_saxpy_3_doit(ap,b1,bp)
    end do

  end subroutine multifab_saxpy_3_c

  subroutine multifab_saxpy_3_cc(a, ia, b1, b, ib, nc, all)
    real(dp_t), intent(in) :: b1
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    integer, intent(in) :: ia, ib
    integer, intent(in), optional :: nc
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    logical, intent(in), optional :: all
    integer :: ii
    logical :: lall

    lall = .false.; if ( present(all) ) lall = all

    do ii = 1, nlocal(a%la)
       if ( lall ) then
          ap => dataptr(a, ii, ia, nc)
          bp => dataptr(b, ii, ib, nc)
       else
          ap => dataptr(a, ii, get_ibox(a,ii), ia, nc)
          bp => dataptr(b, ii, get_ibox(b,ii), ib, nc)
       end if
       call multifab_saxpy_3_doit(ap,b1,bp)
    end do

  end subroutine multifab_saxpy_3_cc

  function multifab_norm_l1_c(mf, comp, nc, mask, all) result(r)
    real(dp_t) :: r
    logical, intent(in), optional :: all
    integer, intent(in) :: comp
    integer, intent(in), optional :: nc
    type(multifab), intent(in) :: mf
    logical, pointer :: lp(:,:,:,:)
    type(lmultifab), intent(in), optional :: mask
    real(dp_t), pointer :: mp(:,:,:,:)
    integer :: i, n
    real(dp_t) :: r1
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0.0_dp_t
    if ( present(mask) ) then
       do i = 1, nlocal(mf%la)
          if ( lall ) then
             lp => dataptr(mask, i, get_pbox(mask, i))
          else
             lp => dataptr(mask, i, get_ibox(mask, i))
          end if
          do n = comp, comp+nc-1
             if ( lall ) then
                mp => dataptr(mf, i, get_pbox(mf, i), n)
             else
                mp => dataptr(mf, i, get_ibox(mf, i), n)
             end if
             r1 = r1 + sum(abs(mp), mask = lp)
          end do
       end do
    else
       do i = 1, nlocal(mf%la)
          if ( lall ) then
             mp => dataptr(mf, i, get_pbox(mf, i), comp, nc)
          else
             mp => dataptr(mf, i, get_ibox(mf, i), comp, nc)
          end if
          r1 = r1 + sum(abs(mp))
       end do
    end if
    call parallel_reduce(r, r1, MPI_SUM)

  end function multifab_norm_l1_c
  function multifab_norm_l1(mf, all) result(r)
    real(dp_t)                    :: r
    type(multifab), intent(in)    :: mf
    logical, intent(in), optional :: all
    r = multifab_norm_l1_c(mf, 1, mf%nc, all = all)
  end function multifab_norm_l1

  function multifab_sum_c(mf, comp, nc, mask, all) result(r)
    real(dp_t) :: r
    integer, intent(in) :: comp
    logical, intent(in), optional :: all
    integer, intent(in), optional :: nc
    type(multifab), intent(in) :: mf
    type(lmultifab), intent(in), optional :: mask
    real(dp_t), pointer :: mp(:,:,:,:)
    logical, pointer :: lp(:,:,:,:)
    integer :: i, n
    real(dp_t) :: r1
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0.0_dp_t
    if ( present(mask) ) then
       do i = 1, nlocal(mf%la)
          if ( lall ) then
             lp => dataptr(mask, i, get_pbox(mask, i))
          else
             lp => dataptr(mask, i, get_ibox(mask, i))
          end if
          do n = comp, comp+nc-1
             if ( lall ) then
                mp => dataptr(mf, i, get_pbox(mf, i), n)
             else
                mp => dataptr(mf, i, get_ibox(mf, i), n)
             end if
             r1 = r1 + sum(mp, mask=lp)
          end do
       end do
    else
       do i = 1, nlocal(mf%la)
          if ( lall ) then
             mp => dataptr(mf, i, get_pbox(mf, i), comp, nc)
          else
             mp => dataptr(mf, i, get_ibox(mf, i), comp, nc)
          end if
          r1 = r1 + sum(mp)
       end do
    end if
    call parallel_reduce(r, r1, MPI_SUM)
  end function multifab_sum_c
  function multifab_sum(mf, mask, all) result(r)
    real(dp_t)                            :: r
    type(multifab), intent(in)            :: mf
    type(lmultifab), intent(in), optional :: mask
    logical, intent(in), optional         :: all
    r = multifab_sum_c(mf, 1, mf%nc, mask, all)
  end function multifab_sum

  function multifab_norm_l2_doit(ap, lp) result(r)

    real(dp_t), pointer        :: ap(:,:,:,:)
    logical, pointer, optional :: lp(:,:,:,:)
    real(dp_t)                 :: r,r1
    integer                    :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    r1 = 0.0_dp_t

    !$OMP PARALLEL PRIVATE(i,j,k,n) REDUCTION(+:r1) IF((hi(3)-lo(3)).ge.7)
    if ( present(lp) ) then
       do n = lo(4), hi(4)
          !$OMP DO
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   if (lp(i,j,k,n)) r1 = r1 + ap(i,j,k,n)*ap(i,j,k,n)
                end do
             end do
          end do
          !$OMP END DO NOWAIT
       end do
    else
       do n = lo(4), hi(4)
          !$OMP DO
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   r1 = r1 + ap(i,j,k,n)*ap(i,j,k,n)
                end do
             end do
          end do
          !$OMP END DO NOWAIT
       end do
    end if
    !$OMP END PARALLEL

    r = r1

  end function multifab_norm_l2_doit

  function multifab_norm_l2_c(mf, comp, nc, mask, all) result(r)
    real(dp_t) :: r
    integer, intent(in) :: comp
    logical, intent(in), optional :: all
    integer, intent(in), optional :: nc
    type(multifab), intent(in) :: mf
    type(lmultifab), intent(in), optional :: mask
    real(dp_t), pointer :: mp(:,:,:,:)
    logical, pointer :: lp(:,:,:,:)
    integer :: i, n
    real(dp_t) :: r1
    logical :: lall
    integer :: lnc
    lnc  = 1; if ( present(nc) ) lnc = nc
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0.0_dp_t

    if ( present(mask) ) then
       do i = 1, nlocal(mf%la)
          if ( lall ) then
             lp => dataptr(mask, i, get_pbox(mask, i))
          else
             lp => dataptr(mask, i, get_ibox(mask, i))
          end if
          do n = comp, comp + lnc - 1
             if ( lall ) then
                mp => dataptr(mf, i, get_pbox(mf, i), n)
             else
                mp => dataptr(mf, i, get_ibox(mf, i), n)
             end if
             r1 = r1 + multifab_norm_l2_doit(mp,lp)
          end do
       end do
    else
       do i = 1, nlocal(mf%la)
          if ( lall ) then
             mp => dataptr(mf, i, get_pbox(mf, i), comp, nc)
          else
             mp => dataptr(mf, i, get_ibox(mf, i), comp, nc)
          end if
          r1 = r1 + multifab_norm_l2_doit(mp)
       end do
    end if
    call parallel_reduce(r, r1, MPI_SUM)
    r = sqrt(r)
  end function multifab_norm_l2_c
  function multifab_norm_l2(mf, mask, all) result(r)
    real(dp_t)                            :: r
    logical, intent(in), optional         :: all
    type(multifab), intent(in)            :: mf
    type(lmultifab), intent(in), optional :: mask
    r = multifab_norm_l2_c(mf, 1, mf%nc, mask, all)
  end function multifab_norm_l2

  function multifab_norm_inf_doit(ap, lp) result(r)

    real(dp_t), pointer        :: ap(:,:,:,:)
    logical, pointer, optional :: lp(:,:,:,:)
    real(dp_t)                 :: r,r1
    integer                    :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    r1 = 0.0_dp_t

    !$OMP PARALLEL PRIVATE(i,j,k,n) REDUCTION(MAX : r1) IF((hi(3)-lo(3)).ge.7)
    if ( present(lp) ) then
       do n = lo(4), hi(4)
          !$OMP DO
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   if (lp(i,j,k,n)) r1 = max(r1,abs(ap(i,j,k,n)))
                end do
             end do
          end do
          !$OMP END DO NOWAIT
       end do
    else
       do n = lo(4), hi(4)
          !$OMP DO
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   r1 = max(r1,abs(ap(i,j,k,n)))
                end do
             end do
          end do
          !$OMP END DO NOWAIT
       end do
    end if
    !$OMP END PARALLEL

    r = r1

  end function multifab_norm_inf_doit

  function multifab_norm_inf_c(mf, comp, nc, mask, all, local, comm) result(r)
    real(dp_t) :: r
    logical, intent(in), optional :: all, local
    integer, intent(in) :: comp
    integer, intent(in), optional :: nc
    type(lmultifab), intent(in), optional :: mask
    type(multifab), intent(in) :: mf
    integer, intent(in), optional :: comm

    logical,    pointer :: lp(:,:,:,:)
    real(dp_t), pointer :: mp(:,:,:,:)
    integer             :: i, n
    real(dp_t)          :: r1
    logical             :: lall, llocal

    lall   = .false.; if ( present(all)   ) lall   = all
    llocal = .false.; if ( present(local) ) llocal = local

    r1 = 0.0_dp_t

    if ( present(mask) ) then
       do i = 1, nlocal(mf%la)
          if ( lall ) then
             lp => dataptr(mask, i, get_pbox(mask, i))
          else
             lp => dataptr(mask, i, get_ibox(mask, i))
          end if
          do n = comp, comp+nc-1
             if ( lall ) then
                mp => dataptr(mf, i, get_pbox(mf, i), n)
             else
                mp => dataptr(mf, i, get_ibox(mf, i), n)
             end if
             r1 = max(r1, multifab_norm_inf_doit(mp,lp))
          end do
       end do
    else
       do i = 1, nlocal(mf%la)
          if ( lall ) then
             mp => dataptr(mf, i, get_pbox(mf, i), comp, nc)
          else
             mp => dataptr(mf, i, get_ibox(mf, i), comp, nc)
          end if
          r1 = max(r1, multifab_norm_inf_doit(mp))
       end do
    end if

    r = r1

    if ( .not. llocal ) then
       if ( present(comm) ) then
          call parallel_reduce(r, r1, MPI_MAX, comm = comm)
       else
          call parallel_reduce(r, r1, MPI_MAX)
       end if
    end if

  end function multifab_norm_inf_c

  function multifab_norm_inf(mf, mask, all, local, comm) result(r)
    real(dp_t)                            :: r
    type(multifab), intent(in)            :: mf
    type(lmultifab), intent(in), optional :: mask
    logical, intent(in), optional         :: all, local
    integer, intent(in), optional         :: comm
    r = multifab_norm_inf_c(mf, 1, mf%nc, mask, all, local, comm)
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
    do i = 1, nlocal(mf%la)
       if ( lall ) then
          mp => dataptr(mf, i, get_pbox(mf, i), comp, nc)
       else
          mp => dataptr(mf, i, get_ibox(mf, i), comp, nc)
       end if
       r1 = max(r1, maxval(abs(mp)))
    end do
    call parallel_reduce(r, r1, MPI_MAX)
  end function imultifab_norm_inf_c
  function imultifab_norm_inf(mf, all) result(r)
    integer                       :: r
    logical, intent(in), optional :: all
    type(imultifab), intent(in)   :: mf
    r = imultifab_norm_inf_c(mf, 1, mf%nc, all)
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
    do i = 1, nlocal(mf%la)
       if ( lall ) then
          mp => dataptr(mf, i, get_pbox(mf, i), comp, nc)
       else
          mp => dataptr(mf, i, get_ibox(mf, i), comp, nc)
       end if
       r1 = r1 + sum(mp)
    end do
    call parallel_reduce(r, r1, MPI_SUM)
  end function imultifab_sum_c
  function imultifab_sum(mf, all) result(r)
    integer :: r
    logical, intent(in), optional :: all
    type(imultifab), intent(in) :: mf
    r = imultifab_sum_c(mf, 1, mf%nc, all)
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
    do i = 1, nlocal(mf%la)
       if ( lall ) then
          mp => dataptr(mf, i, get_pbox(mf, i))
       else
          mp => dataptr(mf, i, get_ibox(mf, i))
       end if
       r1 = r1 + count(mp)
    end do
    call parallel_reduce(r, r1, MPI_SUM)
  end function lmultifab_count

  subroutine multifab_div_div_c_doit(ap, bp)

    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)

    integer :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! ap = ap/bp

    !$OMP PARALLEL PRIVATE(i,j,k,n) IF((hi(3)-lo(3)).ge.7)
    do n = lo(4), hi(4)
       !$OMP DO
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) / bp(i,j,k,n)
             end do
          end do
       end do
       !$OMP END DO NOWAIT
    end do
    !$OMP END PARALLEL

  end subroutine multifab_div_div_c_doit

  subroutine multifab_div_div_s_doit(ap, b)

    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t)          :: b

    integer :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! ap = ap/b

    !$OMP PARALLEL PRIVATE(i,j,k,n) IF((hi(3)-lo(3)).ge.7)
    do n = lo(4), hi(4)
       !$OMP DO
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) / b
             end do
          end do
       end do
       !$OMP END DO NOWAIT
    end do
    !$OMP END PARALLEL

  end subroutine multifab_div_div_s_doit

  subroutine multifab_div_div(a, b, ng)
    type(multifab), intent(inout) :: a
    type(multifab), intent(in   ) :: b
    integer, intent(in), optional :: ng
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i,lng

    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= ng, b%ng >= ng,"not enough ghost cells in multifab_div_div")

    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a,i),lng) )
          bp => dataptr(b, i, grow(get_ibox(b,i),lng) )
       else
          ap => dataptr(a, i, get_ibox(a, i))
          bp => dataptr(b, i, get_ibox(b, i))
       end if
       if ( any(bp == 0.0_dp_t) ) then
          call bl_error("MULTIFAB_DIV_DIV: divide by zero")
       end if
       call multifab_div_div_c_doit(ap, bp)
    end do

  end subroutine multifab_div_div

  subroutine multifab_div_div_s(a, b, ng)
    type(multifab), intent(inout) :: a
    real(dp_t), intent(in)  :: b
    integer, intent(in), optional :: ng
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i,lng

    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= ng,"not enough ghost cells in multifab_div_div_s")
    if ( b == 0.0_dp_t ) then
       call bl_error("MULTIFAB_DIV_DIV_S: divide by zero")
    end if

    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a, i),lng))
       else
          ap => dataptr(a, i, get_ibox(a, i))
       end if
       call multifab_div_div_s_doit(ap, b)
    end do

  end subroutine multifab_div_div_s

  subroutine multifab_div_div_c(a, targ, b, src, nc, ng)
    integer, intent(in) :: targ, src
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i,lng

    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= ng,"not enough ghost cells in multifab_div_div_c")

    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a,i),lng), targ, nc)
          bp => dataptr(b, i, grow(get_ibox(b,i),lng), src, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
          bp => dataptr(b, i, get_ibox(b, i), src, nc)
       end if
       if ( any(bp == 0.0_dp_t) ) then
          call bl_error("MULTIFAB_DIV_DIV: divide by zero")
       end if
       call multifab_div_div_c_doit(ap, bp)
    end do

  end subroutine multifab_div_div_c

  subroutine multifab_div_div_s_c(a, targ, b, nc, ng)
    integer, intent(in) :: targ
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    real(dp_t), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i,lng

    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= ng,"not enough ghost cells in multifab_div_div_s_c")
    if ( b == 0.0_dp_t ) then
       call bl_error("MULTIFAB_DIV_DIV_S_C: divide by zero")
    end if

    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a,i),lng), targ, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
       end if
       call multifab_div_div_s_doit(ap, b)
    end do

  end subroutine multifab_div_div_s_c

  subroutine multifab_div_s_c(a, ia, b, ib, val, nc, ng)
    integer, intent(in) :: ia, ib
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)    :: b
    real(dp_t), intent(in)  :: val
    real(dp_t), pointer :: ap(:,:,:,:), bp(:,:,:,:)
    integer :: i,lng

    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= ng,"not enough ghost cells in multifab_div_s_c")
    if ( val == 0.0_dp_t ) then
       call bl_error("MULTIFAB_DIV_DIV_S_C: divide by zero")
    end if
    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a, i), lng), ia, nc)
          bp => dataptr(a, i, grow(get_ibox(b, i), lng), ib, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), ia, nc)
          bp => dataptr(a, i, get_ibox(b, i), ib, nc)
       end if
       ap = bp/val
    end do
  end subroutine multifab_div_s_c

  subroutine multifab_mult_mult_c_doit(ap, bp)
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer             :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! ap = ap*bp

    !$OMP PARALLEL PRIVATE(i,j,k,n) IF((hi(3)-lo(3)).ge.7)
    do n = lo(4), hi(4)
       !$OMP DO
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) * bp(i,j,k,n)
             end do
          end do
       end do
       !$OMP END DO NOWAIT
    end do
    !$OMP END PARALLEL

  end subroutine multifab_mult_mult_c_doit

  subroutine multifab_mult_mult_s_doit(ap, b)
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t)          :: b
    integer             :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! ap = ap*b

    !$OMP PARALLEL PRIVATE(i,j,k,n) IF((hi(3)-lo(3)).ge.7)
    do n = lo(4), hi(4)
       !$OMP DO
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) * b
             end do
          end do
       end do
       !$OMP END DO NOWAIT
    end do
    !$OMP END PARALLEL

  end subroutine multifab_mult_mult_s_doit

  subroutine multifab_mult_mult(a, b, ng)
    type(multifab), intent(inout) :: a
    type(multifab), intent(in   ) :: b
    integer, intent(in), optional :: ng
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i,lng
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= ng,"not enough ghost cells in multifab_mult_mult")

    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a, i),lng))
          bp => dataptr(b, i, grow(get_ibox(b, i),lng))
       else
          ap => dataptr(a, i, get_ibox(a, i))
          bp => dataptr(b, i, get_ibox(b, i))
       end if
       call multifab_mult_mult_c_doit(ap, bp)
    end do

  end subroutine multifab_mult_mult
  subroutine multifab_mult_mult_s(a, b, ng)
    type(multifab), intent(inout) :: a
    real(dp_t)    , intent(in   ) :: b
    integer, intent(in), optional :: ng
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i,lng
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= ng,"not enough ghost cells in multifab_mult_mult_s")

    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a, i),lng))
       else
          ap => dataptr(a, i, get_ibox(a, i))
       end if
       call multifab_mult_mult_s_doit(ap, b)
    end do

  end subroutine multifab_mult_mult_s

  subroutine multifab_mult_mult_c(a, targ, b, src, nc, ng)
    integer, intent(in) :: targ, src
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i,lng
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= ng,"not enough ghost cells in multifab_mult_mult_c")

    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a, i),lng), targ, nc)
          bp => dataptr(b, i, grow(get_ibox(b, i),lng), src, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
          bp => dataptr(b, i, get_ibox(b, i), src, nc)
       end if
       call multifab_mult_mult_c_doit(ap, bp)
    end do

  end subroutine multifab_mult_mult_c

  subroutine multifab_mult_mult_s_c(a, targ, b, nc, ng)
    integer, intent(in) :: targ
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    real(dp_t), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i,lng
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= ng,"not enough ghost cells in multifab_mult_mult_s_c")

    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a, i),lng), targ, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
       end if
       call multifab_mult_mult_s_doit(ap, b)
    end do

  end subroutine multifab_mult_mult_s_c

  subroutine multifab_mult_s_c(a, ia, b, ib, val, nc, ng)
    integer, intent(in) :: ia, ib
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)    :: b
    real(dp_t), intent(in)  :: val
    real(dp_t), pointer :: ap(:,:,:,:), bp(:,:,:,:)
    integer :: i,lng
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= ng,"not enough ghost cells in multifab_mult_s_c")

    !$OMP PARALLEL DO PRIVATE(i,ap,bp)
    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a, i), lng), ia, nc)
          bp => dataptr(a, i, grow(get_ibox(b, i), lng), ib, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), ia, nc)
          bp => dataptr(a, i, get_ibox(b, i), ib, nc)
       end if
       ap = bp * val
    end do
    !$OMP END PARALLEL DO

  end subroutine multifab_mult_s_c

  subroutine multifab_sub_sub_c_doit(ap, bp)

    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)

    integer :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! ap = ap - bp

    !$OMP PARALLEL PRIVATE(i,j,k,n) IF((hi(3)-lo(3)).ge.7)
    do n = lo(4), hi(4)
       !$OMP DO
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) - bp(i,j,k,n)
             end do
          end do
       end do
       !$OMP END DO NOWAIT
    end do
    !$OMP END PARALLEL

  end subroutine multifab_sub_sub_c_doit

  subroutine multifab_sub_sub_s_doit(ap, b)

    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t)          :: b

    integer :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! ap = ap - b

    !$OMP PARALLEL PRIVATE(i,j,k,n) IF((hi(3)-lo(3)).ge.7)
    do n = lo(4), hi(4)
       !$OMP DO
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) - b
             end do
          end do
       end do
       !$OMP END DO NOWAIT
    end do
    !$OMP END PARALLEL

  end subroutine multifab_sub_sub_s_doit

  subroutine multifab_sub_sub(a, b, ng)
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    integer, intent(in), optional :: ng
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i,lng
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= ng, b%ng >= ng, "not enough ghost cells in multifab_sub_sub")

    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a, i),ng))
          bp => dataptr(b, i, grow(get_ibox(b, i),ng))
       else
          ap => dataptr(a, i, get_ibox(a, i))
          bp => dataptr(b, i, get_ibox(b, i))
       end if
       call multifab_sub_sub_c_doit(ap, bp)
    end do

  end subroutine multifab_sub_sub

  subroutine multifab_sub_sub_s(a, b, ng)
    type(multifab), intent(inout) :: a
    real(dp_t), intent(in)  :: b
    integer, intent(in), optional :: ng
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i,lng
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= ng, "not enough ghost cells in multifab_sub_sub_s")

    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a, i),ng))
       else
          ap => dataptr(a, i, get_ibox(a, i))
       end if
       call multifab_sub_sub_s_doit(ap, b)
    end do

  end subroutine multifab_sub_sub_s

  subroutine multifab_sub_sub_c(a, targ, b, src, nc, ng)
    integer, intent(in) :: targ, src
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i,lng
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= ng, b%ng >= ng, "not enough ghost cells in multifab_sub_sub_c")

    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a, i),ng), targ, nc)
          bp => dataptr(b, i, grow(get_ibox(b, i),ng), src, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
          bp => dataptr(b, i, get_ibox(b, i), src, nc)
       end if
       call multifab_sub_sub_c_doit(ap, bp)
    end do

  end subroutine multifab_sub_sub_c

  subroutine multifab_sub_sub_s_c(a, targ, b, nc, ng)
    integer, intent(in) :: targ
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    real(dp_t), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i,lng
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= ng, "not enough ghost cells in multifab_sub_sub_s_c")

    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a, i),ng), targ, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
       end if
       call multifab_sub_sub_s_doit(ap, b)
    end do

  end subroutine multifab_sub_sub_s_c

  subroutine multifab_plus_plus_c_doit(ap, bp)

    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)

    integer :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! ap = ap + bp

    !$OMP PARALLEL PRIVATE(i,j,k,n) IF((hi(3)-lo(3)).ge.7)
    do n = lo(4), hi(4)
       !$OMP DO
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) + bp(i,j,k,n)
             end do
          end do
       end do
       !$OMP END DO NOWAIT
    end do
    !$OMP END PARALLEL

  end subroutine multifab_plus_plus_c_doit

  subroutine multifab_plus_plus_c(a, dst, b, src, nc, ng)
    integer, intent(in) :: dst, src
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i,lng

    lng = 0; if ( present(ng) ) lng = ng

    if ( lng > 0 ) call bl_assert(a%ng >= ng, b%ng >= ng,"not enough ghost cells in multifab_plus_plus_c")

    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a,i),lng), dst, nc)
          bp => dataptr(b, i, grow(get_ibox(b,i),lng), src, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), dst, nc)
          bp => dataptr(b, i, get_ibox(b, i), src, nc)
       end if
       call multifab_plus_plus_c_doit(ap, bp)
    end do

  end subroutine multifab_plus_plus_c

  subroutine multifab_plus_plus(a, b, ng)
    type(multifab), intent(inout) :: a
    type(multifab), intent(in   ) :: b
    integer, intent(in), optional :: ng
    if ( present(ng) ) then
      call multifab_plus_plus_c(a, 1, b, 1, a%nc, ng)
    else
      call multifab_plus_plus_c(a, 1, b, 1, a%nc)
    end if
  end subroutine multifab_plus_plus

  subroutine multifab_plus_plus_s_doit(ap, b)

    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t)          :: b

    integer :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! ap = ap + b

    !$OMP PARALLEL PRIVATE(i,j,k,n) IF((hi(3)-lo(3)).ge.7)
    do n = lo(4), hi(4)
       !$OMP DO
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) + b
             end do
          end do
       end do
       !$OMP END DO NOWAIT
    end do
    !$OMP END PARALLEL

  end subroutine multifab_plus_plus_s_doit

  subroutine multifab_plus_plus_s_c(a, dst, b, nc, ng)
    integer, intent(in) :: dst
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    real(dp_t), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i,lng

    lng = 0; if ( present(ng) ) lng = ng

    do i = 1, nlocal(a%la)
       if ( lng > 0 ) then
          ap => dataptr(a, i, grow(get_ibox(a,i),lng), dst, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), dst, nc)
       end if
       call multifab_plus_plus_s_doit(ap, b)
    end do

  end subroutine multifab_plus_plus_s_c

  subroutine multifab_plus_plus_s(a, b, ng)
    type(multifab), intent(inout) :: a
    real(dp_t), intent(in)  :: b
    integer, intent(in), optional :: ng
    if ( present(ng) ) then
      call multifab_plus_plus_s_c(a, 1, b, a%nc, ng)
    else
      call multifab_plus_plus_s_c(a, 1, b, a%nc)
    end if
  end subroutine multifab_plus_plus_s

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
    real(dp_t), pointer :: mp(:,:,:,:)
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

    do i = 1, nlocal(mf%la)
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
      real(dp_t), intent(inout)  :: f(lo(1):,:)
      real(dp_t), intent(in) :: x(xo(1):,:)
      integer :: i
      do i = max(lbound(f,1),lbound(x,1)), min(ubound(f,1),ubound(x,1))
         f(i,:) = x(i,:)
      end do
    end subroutine c_1d
    subroutine c_2d(f, lo, x, xo)
      integer, intent(in) :: lo(:), xo(:)
      real(dp_t), intent(inout)  :: f(lo(1):,lo(2):,:)
      real(dp_t), intent(in) :: x(xo(1):,xo(2):,:)
      integer :: i, j
      do j = max(lbound(f,2),lbound(x,2)), min(ubound(f,2),ubound(x,2))
         do i = max(lbound(f,1),lbound(x,1)), min(ubound(f,1),ubound(x,1))
            f(i,j,:) = x(i,j,:)
         end do
      end do
    end subroutine c_2d
    subroutine c_3d(f, lo, x, xo)
      integer, intent(in) :: lo(:), xo(:)
      real(dp_t), intent(inout) :: f(lo(1):,lo(2):,lo(3):,:)
      real(dp_t), intent(in) :: x(xo(1):,xo(2):,xo(3):,:)
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

  function multifab_min(mf, all) result(r)
    real(kind=dp_t) :: r
    type(multifab), intent(in) :: mf
    logical, intent(in), optional :: all
    integer :: i
    real(kind=dp_t) :: r1
    r1 = +Huge(r1)
    do i = 1, nlocal(mf%la)
       r1 = min(r1,min_val(mf%fbs(i), all))
    end do
    call parallel_reduce(r, r1, MPI_MIN)
  end function multifab_min
  function multifab_min_c(mf, c, nc, all) result(r)
    real(kind=dp_t) :: r
    type(multifab), intent(in) :: mf
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    real(kind=dp_t) :: r1
    integer :: i
    r1 = +Huge(r1)
    do i = 1, nlocal(mf%la)
       r1 = min(r1, min_val(mf%fbs(i), c, nc, all))
    end do
    call parallel_reduce(r, r1, MPI_MIN)
  end function multifab_min_c
  
  function multifab_max(mf, all, allow_empty) result(r)
    real(kind=dp_t) :: r
    type(multifab), intent(in) :: mf
    logical, intent(in), optional :: all, allow_empty
    logical :: lallow
    integer :: i
    real(kind=dp_t) :: r1
    r1 = -Huge(r1)
    lallow = .false.; if (present(allow_empty)) lallow=allow_empty
    if (lallow) then
       do i = 1, nlocal(mf%la)
          if (empty(get_box(mf%fbs(i)))) cycle
          r1 = max(r1, max_val(mf%fbs(i), all))
       end do
    else
       do i = 1, nlocal(mf%la)
          r1 = max(r1, max_val(mf%fbs(i), all))
       end do
    endif
    call parallel_reduce(r, r1, MPI_MAX)
  end function multifab_max
  function multifab_max_c(mf, c, nc, all) result(r)
    real(kind=dp_t) :: r
    type(multifab), intent(in) :: mf
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    integer :: i
    real(kind=dp_t) :: r1
    r1 = -Huge(r1)
    do i = 1, nlocal(mf%la)
       r1 = max(r1, max_val(mf%fbs(i), c, nc, all))
    end do
    call parallel_reduce(r, r1, MPI_MAX)
  end function multifab_max_c
  
end module multifab_module
