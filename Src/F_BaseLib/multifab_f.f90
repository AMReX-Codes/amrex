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
     logical :: sent = .false.
     logical :: rcvd = .false.
     integer, pointer :: send_request(:) => Null()
     integer, pointer :: recv_request(:) => Null()
     real(dp_t), pointer :: send_buffer(:) => Null()
     real(dp_t), pointer :: recv_buffer(:) => Null()
  end type mf_fb_data

  type mfiter
     integer :: dim      = 0
     integer :: ng       = 0
     logical :: nodal(3) = .false.
     integer :: it       = 0    ! current tile index local to this thread
     integer :: ntiles   = 0
     type(tilearray)  :: ta
     type(layout_rep), pointer :: lap => Null()
     logical :: built    = .false.
  end type mfiter

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

  interface mfiter_build
     module procedure  multifab_iter_build
     module procedure zmultifab_iter_build
     module procedure lmultifab_iter_build
     module procedure imultifab_iter_build
  end interface mfiter_build

  interface is_equal
     module procedure multifab_equal
  end interface is_equal

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

  private :: multifab_iter_build, zmultifab_iter_build, lmultifab_iter_build, imultifab_iter_build
  private :: iter_build_doit

  private :: multifab_equal

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

  subroutine multifab_build(mf, la, nc, ng, nodal, stencil, fab_alloc)
    type(multifab), intent(out)   :: mf
    type(layout),   intent(in )   :: la
    integer, intent(in), optional :: nc, ng
    logical, intent(in), optional :: nodal(:), stencil, fab_alloc
    integer :: i, lnc, lng
    logical :: lfab_alloc
    if ( built_q(mf) ) call bl_error("MULTIFAB_BUILD: already built")
    lng = 0; if ( present(ng) ) lng = ng
    lnc = 1; if ( present(nc) ) lnc = nc
    mf%dim = get_dim(la)
    mf%la  = la
    mf%nc  = lnc
    mf%ng  = lng
    allocate(mf%nodal(mf%dim))
    mf%nodal = .False.; if ( present(nodal) ) mf%nodal = nodal(1:mf%dim)

    lfab_alloc = .true.;  if ( present(fab_alloc) ) lfab_alloc = fab_alloc

    if (lfab_alloc) then
       allocate(mf%fbs(nlocal(mf%la)))

       do i = 1, nlocal(mf%la)
          call fab_build( &
               mf%fbs(i), get_box(mf%la, global_index(mf%la,i)), &
               mf%nc, mf%ng, mf%nodal,  &
               alloc = .true., stencil = stencil)
       end do
       call mem_stats_alloc(multifab_ms, volume(mf, all = .TRUE.))
    end if
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
       call ifab_build(mf%fbs(i), get_box(mf%la, global_index(mf%la,i)), &
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
       call lfab_build(mf%fbs(i), get_box(mf%la, global_index(mf%la,i)), &
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
       call zfab_build(mf%fbs(i), get_box(mf%la, global_index(mf%la,i)), &
            mf%nc, mf%ng, mf%nodal, alloc = .true.)
    end do
    call mem_stats_alloc(zmultifab_ms, volume(mf, all = .TRUE.))
  end subroutine zmultifab_build

  subroutine multifab_build_copy(m1, m2)
    type(multifab), intent(inout) :: m1
    type(multifab), intent(in) :: m2
    integer :: i
    if ( built_q(m1) ) call bl_error("MULTIFAB_BUILD_COPY: already built")
    if ( built_q(m1) ) call multifab_destroy(m1)
    m1%dim = m2%dim
    m1%la  = m2%la
    m1%nc  = m2%nc
    m1%ng  = m2%ng
    allocate(m1%nodal(m1%dim))
    m1%nodal = m2%nodal
    allocate(m1%fbs(nlocal(m1%la)))
    do i = 1, nlocal(m1%la)
       call fab_build(m1%fbs(i), get_box(m2%fbs(i)), m1%nc, m1%ng, m1%nodal)
    end do
    call multifab_copy(m1, m2, ng=m1%ng)
    call mem_stats_alloc(multifab_ms, volume(m1, all = .TRUE.))
  end subroutine multifab_build_copy
  subroutine imultifab_build_copy(m1, m2)
    type(imultifab), intent(inout) :: m1
    type(imultifab), intent(in) :: m2
    integer :: i
    if ( built_q(m1) ) call bl_error("IMULTIFAB_BUILD_COPY: already built")
    if ( built_q(m1) ) call imultifab_destroy(m1)
    m1%dim = m2%dim
    m1%la  = m2%la
    m1%nc  = m2%nc
    m1%ng  = m2%ng
    allocate(m1%nodal(m1%dim))
    m1%nodal = m2%nodal
    allocate(m1%fbs(nlocal(m1%la)))
    do i = 1, nlocal(m1%la)
       call ifab_build(m1%fbs(i), get_box(m2%fbs(i)), m1%nc, m1%ng, m1%nodal)
    end do
    call imultifab_copy(m1, m2, ng=m1%ng)
    call mem_stats_alloc(imultifab_ms, volume(m1, all = .TRUE.))
  end subroutine imultifab_build_copy
  subroutine lmultifab_build_copy(m1, m2)
    type(lmultifab), intent(inout) :: m1
    type(lmultifab), intent(in) :: m2
    integer :: i
    if ( built_q(m1) ) call bl_error("LMULTIFAB_BUILD_COPY: already built")
    if ( built_q(m1) ) call lmultifab_destroy(m1)
    m1%dim = m2%dim
    m1%la  = m2%la
    m1%nc  = m2%nc
    m1%ng  = m2%ng
    allocate(m1%nodal(m1%dim))
    m1%nodal = m2%nodal
    allocate(m1%fbs(nlocal(m1%la)))
    do i = 1, nlocal(m1%la)
       call lfab_build(m1%fbs(i), get_box(m2%fbs(i)), m1%nc, m1%ng, m1%nodal)
    end do
    call lmultifab_copy(m1, m2, ng=m1%ng)
    call mem_stats_alloc(lmultifab_ms, volume(m1, all = .TRUE.))
  end subroutine lmultifab_build_copy
  subroutine zmultifab_build_copy(m1, m2)
    type(zmultifab), intent(inout) :: m1
    type(zmultifab), intent(in) :: m2
    integer :: i
    if ( built_q(m1) ) call bl_error("ZMULTIFAB_BUILD_COPY: already built")
    if ( built_q(m1) ) call zmultifab_destroy(m1)
    m1%dim = m2%dim
    m1%la  = m2%la
    m1%nc  = m2%nc
    m1%ng  = m2%ng
    allocate(m1%nodal(m1%dim))
    m1%nodal = m2%nodal
    allocate(m1%fbs(nlocal(m1%la)))
    do i = 1, nlocal(m1%la)
       call zfab_build(m1%fbs(i), get_box(m2%fbs(i)), m1%nc, m1%ng, m1%nodal)
    end do
    call zmultifab_copy(m1, m2, m1%ng)
    call mem_stats_alloc(zmultifab_ms, volume(m1, all = .TRUE.))
  end subroutine zmultifab_build_copy

  subroutine multifab_destroy(mf)
    type(multifab), intent(inout) :: mf
    integer :: i
    if (associated(mf%fbs)) then
       call mem_stats_dealloc(multifab_ms, volume(mf, all = .TRUE.))
       do i = 1, nlocal(mf%la)
          call fab_destroy(mf%fbs(i))
       end do
       deallocate(mf%fbs)
    end if
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
       call ifab_destroy(mf%fbs(i))
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
       call lfab_destroy(mf%fbs(i))
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
       call zfab_destroy(mf%fbs(i))
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
       !$omp parallel do reduction(+:r)
       do i = 1, nboxes(mf%la)
          r = r + volume(grow(box_nodalize(get_box(mf%la,i),mf%nodal),mf%ng))
       end do
       !$omp end parallel do
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
       !$omp parallel do reduction(+:r)
       do i = 1, nboxes(mf%la)
          r = r + volume(grow(box_nodalize(get_box(mf%la,i),mf%nodal),mf%ng))
       end do
       !$omp end parallel do
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
       !$omp parallel do reduction(+:r)
       do i = 1, nboxes(mf%la)
          r = r + volume(grow(box_nodalize(get_box(mf%la,i),mf%nodal),mf%ng))
       end do
       !$omp end parallel do
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
       !$omp parallel do reduction(+:r)
       do i = 1, nboxes(mf%la)
          r = r + volume(grow(box_nodalize(get_box(mf%la,i),mf%nodal),mf%ng))
       end do
       !$omp end parallel do
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

  subroutine multifab_setval(mf, val, all)
    type(multifab), intent(inout) :: mf
    real(dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    call multifab_setval_c(mf, val, 1, mf%nc, all)
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
    type(mfiter) :: mfi
    lall = .FALSE.; if ( present(all) ) lall = all
    !$omp parallel private(i,bx1,mfi)
    call mfiter_build(mfi, mf, .true.)
    do while (next_tile(mfi,i))
       if ( lall ) then
          bx1 = intersection(bx, get_growntilebox(mfi))
       else
          bx1 = intersection(bx, get_tilebox(mfi))
       end if
       if ( .not. empty(bx1) ) then
          call setval(mf%fbs(i), val, bx1, c, nc)
       end if
    end do
    !$omp end parallel
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
    type(mfiter) :: mfi
    lall = .FALSE.; if ( present(all) ) lall = all
    !$omp parallel private(i,bx1,mfi)
    call mfiter_build(mfi, mf, .true.)
    do while (next_tile(mfi,i))
       if ( lall ) then
          bx1 = intersection(bx, get_growntilebox(mfi))
       else
          bx1 = intersection(bx, get_tilebox(mfi))
       end if
       if ( .not. empty(bx1) ) then
          call setval(mf%fbs(i), val, bx1, c, nc)
       end if
    end do
    !$omp end parallel
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
    type(mfiter) :: mfi
    lall = .FALSE.; if ( present(all) ) lall = all
    !$omp parallel private(i,bx1,mfi)
    call mfiter_build(mfi, mf, .true.)
    do while (next_tile(mfi,i))
       if ( lall ) then
          bx1 = intersection(bx, get_growntilebox(mfi))
       else
          bx1 = intersection(bx, get_tilebox(mfi))
       end if
       if ( .not. empty(bx1) ) then
          call setval(mf%fbs(i), val, bx1, c, nc)
       end if
    end do
    !$omp end parallel
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
    type(mfiter) :: mfi
    lall = .FALSE.; if ( present(all) ) lall = all
    !$omp parallel private(i,bx1,mfi)
    call mfiter_build(mfi, mf, .true.)
    do while (next_tile(mfi,i))
       if ( lall ) then
          bx1 = intersection(bx, get_growntilebox(mfi))
       else
          bx1 = intersection(bx, get_tilebox(mfi))
       end if
       if ( .not. empty(bx1) ) then
          call setval(mf%fbs(i), val, bx1, c, nc)
       end if
    end do
    !$omp end parallel
  end subroutine zmultifab_setval_bx_c

  subroutine multifab_setval_c(mf, val, c, nc, all)
    type(multifab), intent(inout) :: mf
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    real(dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    integer :: i
    logical lall
    type(mfiter) :: mfi
    type(bl_prof_timer), save :: bpt
    call build(bpt, "mf_setval_c")
    lall = .FALSE.; if ( present(all) ) lall = all
    !$omp parallel private(mfi,i)
    call mfiter_build(mfi,mf,.true.)
    do while (next_tile(mfi,i))
       if ( lall ) then
          call setval(mf%fbs(i), val, get_growntilebox(mfi), c, nc)
       else
          call setval(mf%fbs(i), val, get_tilebox(mfi), c, nc)
       end if
    end do
    !$omp end parallel
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
    type(mfiter) :: mfi
    type(bl_prof_timer), save :: bpt
    call build(bpt, "imf_setval_c")
    lall = .FALSE.; if ( present(all) ) lall = all
    !$omp parallel private(mfi,i)
    call mfiter_build(mfi,mf,.true.)
    do while (next_tile(mfi,i))
       if ( lall ) then
          call setval(mf%fbs(i), val, get_growntilebox(mfi), c, nc)
       else
          call setval(mf%fbs(i), val, get_tilebox(mfi), c, nc)
       end if
    end do
    !$omp end parallel
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
    type(mfiter) :: mfi
    type(bl_prof_timer), save :: bpt
    call build(bpt, "lmf_setval_c")
    lall = .FALSE.; if ( present(all) ) lall = all
    !$omp parallel private(mfi,i)
    call mfiter_build(mfi,mf,.true.)
    do while (next_tile(mfi,i))
       if ( lall ) then
          call setval(mf%fbs(i), val, get_growntilebox(mfi), c, nc)
       else
          call setval(mf%fbs(i), val, get_tilebox(mfi), c, nc)
       end if
    end do
    !$omp end parallel
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
    type(mfiter) :: mfi
    lall = .FALSE.; if ( present(all) ) lall = all
    !$omp parallel private(mfi,i)
    call mfiter_build(mfi,mf,.true.)
    do while (next_tile(mfi,i))
       if ( lall ) then
          call setval(mf%fbs(i), val, get_growntilebox(mfi), c, nc)
       else
          call setval(mf%fbs(i), val, get_tilebox(mfi), c, nc)
       end if
    end do
    !$omp end parallel
  end subroutine zmultifab_setval_c

  subroutine multifab_set_corner(mf, val)
    type(multifab), intent(inout) :: mf
    real(dp_t), intent(in) :: val
    integer :: i
    type(mfiter) :: mfi
    if (mf%ng > 0) then
       !$omp parallel private(mfi,i)
       call mfiter_build(mfi,mf,.true.)
       do while (next_tile(mfi,i))
          call fab_set_corner(mf%fbs(i), val, get_growntilebox(mfi), 1, mf%nc)
       end do
       !$omp end parallel
    end if
  end subroutine multifab_set_corner

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
                dst(c) = src(i,j,k,n)
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
    integer                 :: tag
    integer                 :: i, ii, jj, np, sh(MAX_SPACEDIM+1)
    type(boxassoc)          :: bxasc
    real(dp_t), allocatable :: g_snd_d(:), g_rcv_d(:)

    bxasc = layout_boxassoc(mf%la, ng, mf%nodal, lcross, idim)

    !$OMP PARALLEL DO PRIVATE(i,ii,jj,p1,p2) if (bxasc%l_con%threadsafe)
    do i = 1, bxasc%l_con%ncpy
       ii  =  bxasc%l_con%cpy(i)%lnd
       jj  =  bxasc%l_con%cpy(i)%lns
       p1  => dataptr(mf%fbs(ii), bxasc%l_con%cpy(i)%dbx, c, nc)
       p2  => dataptr(mf%fbs(jj), bxasc%l_con%cpy(i)%sbx, c, nc)
       call cpy_d(p1,p2)
    end do
    !$OMP END PARALLEL DO

    np = parallel_nprocs()

    if (np == 1) return

    tag = parallel_tag()

    allocate(g_snd_d(nc*bxasc%r_con%svol))
    allocate(g_rcv_d(nc*bxasc%r_con%rvol))

    !$OMP PARALLEL DO PRIVATE(i,p)
    do i = 1, bxasc%r_con%nsnd
       p => dataptr(mf, bxasc%r_con%snd(i)%lns, bxasc%r_con%snd(i)%sbx, c, nc)
       call reshape_d_4_1(g_snd_d, 1 + nc*bxasc%r_con%snd(i)%pv, p)
    end do
    !$OMP END PARALLEL DO

    if (bxasc%r_con%nrp+bxasc%r_con%nsp .gt. 0) then
       allocate(rst(bxasc%r_con%nrp+bxasc%r_con%nsp))
       do i = 1, bxasc%r_con%nrp
          rst(i) = parallel_irecv_dv(g_rcv_d(1+nc*bxasc%r_con%rtr(i)%pv), &
               nc*bxasc%r_con%rtr(i)%sz, bxasc%r_con%rtr(i)%pr, tag)
       end do
       do i = 1, bxasc%r_con%nsp
          rst(i+bxasc%r_con%nrp) = parallel_isend_dv(g_snd_d(1+nc*bxasc%r_con%str(i)%pv), &
               nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
       end do
       call parallel_wait(rst)
    end if

    !$OMP PARALLEL DO PRIVATE(i,sh,p) if (bxasc%r_con%threadsafe)
    do i = 1, bxasc%r_con%nrcv
       sh = bxasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mf, bxasc%r_con%rcv(i)%lnd, bxasc%r_con%rcv(i)%dbx, c, nc)
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
    integer                 :: i, ii, jj, np, istart, iend, nsize, tag
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

    !$OMP PARALLEL DO PRIVATE(i,ii,jj,p1,p2) if (bxasc%l_con%threadsafe)
    do i = 1, bxasc%l_con%ncpy
       ii  =  bxasc%l_con%cpy(i)%lnd
       jj  =  bxasc%l_con%cpy(i)%lns
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

    tag = parallel_tag()

    allocate(fb_data%send_buffer(nc*bxasc%r_con%svol))
    allocate(fb_data%recv_buffer(nc*bxasc%r_con%rvol))

    !$OMP PARALLEL DO PRIVATE(i,p)
    do i = 1, bxasc%r_con%nsnd
       p => dataptr(mf, bxasc%r_con%snd(i)%lns, bxasc%r_con%snd(i)%sbx, c, nc)
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
            nsize, bxasc%r_con%rtr(i)%pr, tag)
    end do

    do i = 1, bxasc%r_con%nsp
       istart = nc*bxasc%r_con%str(i)%pv + 1
       nsize = nc*bxasc%r_con%str(i)%sz
       iend = istart + nsize - 1
       fb_data%send_request(i) = parallel_isend_dv(fb_data%send_buffer(istart:iend), &
            nsize, bxasc%r_con%str(i)%pr, tag)
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

       !$omp parallel do private(i,sh,p) if (bxasc%r_con%threadsafe)
       do i = 1, bxasc%r_con%nrcv
          sh = bxasc%r_con%rcv(i)%sh
          sh(4) = nc
          p => dataptr(mf, bxasc%r_con%rcv(i)%lnd, bxasc%r_con%rcv(i)%dbx, c, nc)
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

       !$omp parallel do private(i,sh,p) if (bxasc%r_con%threadsafe)
       do i = 1, bxasc%r_con%nrcv
          sh = bxasc%r_con%rcv(i)%sh
          sh(4) = nc
          p => dataptr(mf, bxasc%r_con%rcv(i)%lnd, bxasc%r_con%rcv(i)%dbx, c, nc)
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

       !$omp parallel do private(i,sh,p) if (bxasc%r_con%threadsafe)
       do i = 1, bxasc%r_con%nrcv
          sh = bxasc%r_con%rcv(i)%sh
          sh(4) = nc
          p => dataptr(mf, bxasc%r_con%rcv(i)%lnd, bxasc%r_con%rcv(i)%dbx, c, nc)
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
    integer              :: tag
    integer              :: i, ii, jj, np, sh(MAX_SPACEDIM+1)
    type(boxassoc)       :: bxasc
    integer, allocatable :: g_snd_i(:), g_rcv_i(:)

    bxasc = layout_boxassoc(mf%la, ng, mf%nodal, lcross)

    !$OMP PARALLEL DO PRIVATE(i,ii,jj,p1,p2)
    do i = 1, bxasc%l_con%ncpy
       ii  =  bxasc%l_con%cpy(i)%lnd
       jj  =  bxasc%l_con%cpy(i)%lns
       p1  => dataptr(mf%fbs(ii), bxasc%l_con%cpy(i)%dbx, c, nc)
       p2  => dataptr(mf%fbs(jj), bxasc%l_con%cpy(i)%sbx, c, nc)
       call cpy_i(p1,p2)
    end do
    !$OMP END PARALLEL DO

    np = parallel_nprocs()

    if (np == 1) return

    tag = parallel_tag()

    allocate(g_snd_i(nc*bxasc%r_con%svol))
    allocate(g_rcv_i(nc*bxasc%r_con%rvol))

    !$OMP PARALLEL DO PRIVATE(i,p)
    do i = 1, bxasc%r_con%nsnd
       p => dataptr(mf, bxasc%r_con%snd(i)%lns, bxasc%r_con%snd(i)%sbx, c, nc)
       call reshape_i_4_1(g_snd_i, 1 + nc*bxasc%r_con%snd(i)%pv, p)
    end do
    !$OMP END PARALLEL DO

    if (bxasc%r_con%nrp+bxasc%r_con%nsp .gt. 0) then
       allocate(rst(bxasc%r_con%nrp+bxasc%r_con%nsp))
       do i = 1, bxasc%r_con%nrp
          rst(i) = parallel_irecv_iv(g_rcv_i(1+nc*bxasc%r_con%rtr(i)%pv), &
               nc*bxasc%r_con%rtr(i)%sz, bxasc%r_con%rtr(i)%pr, tag)
       end do
       do i = 1, bxasc%r_con%nsp
          rst(i+bxasc%r_con%nrp) = parallel_isend_iv(g_snd_i(1+nc*bxasc%r_con%str(i)%pv), &
               nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
       end do
       call parallel_wait(rst)
    end if

    !$OMP PARALLEL DO PRIVATE(i,sh,p)
    do i = 1, bxasc%r_con%nrcv
       sh = bxasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mf, bxasc%r_con%rcv(i)%lnd, bxasc%r_con%rcv(i)%dbx, c, nc)
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
    integer              :: tag
    integer              :: i, ii, jj, np, sh(MAX_SPACEDIM+1)
    type(boxassoc)       :: bxasc
    logical, allocatable :: g_snd_l(:), g_rcv_l(:)

    bxasc = layout_boxassoc(mf%la, ng, mf%nodal, lcross)

    !$OMP PARALLEL DO PRIVATE(i,ii,jj,p1,p2)
    do i = 1, bxasc%l_con%ncpy
       ii  =  bxasc%l_con%cpy(i)%lnd
       jj  =  bxasc%l_con%cpy(i)%lns
       p1  => dataptr(mf%fbs(ii), bxasc%l_con%cpy(i)%dbx, c, nc)
       p2  => dataptr(mf%fbs(jj), bxasc%l_con%cpy(i)%sbx, c, nc)
       call cpy_l(p1,p2)
    end do
    !$OMP END PARALLEL DO

    np = parallel_nprocs()

    if (np == 1) return

    tag = parallel_tag()

    allocate(g_snd_l(nc*bxasc%r_con%svol))
    allocate(g_rcv_l(nc*bxasc%r_con%rvol))

    !$OMP PARALLEL DO PRIVATE(i,p)
    do i = 1, bxasc%r_con%nsnd
       p => dataptr(mf, bxasc%r_con%snd(i)%lns, bxasc%r_con%snd(i)%sbx, c, nc)
       call reshape_l_4_1(g_snd_l, 1 + nc*bxasc%r_con%snd(i)%pv, p)
    end do
    !$OMP END PARALLEL DO

    if (bxasc%r_con%nrp+bxasc%r_con%nsp .gt. 0) then
       allocate(rst(bxasc%r_con%nrp+bxasc%r_con%nsp))
       do i = 1, bxasc%r_con%nrp
          rst(i) = parallel_irecv_lv(g_rcv_l(1+nc*bxasc%r_con%rtr(i)%pv), &
               nc*bxasc%r_con%rtr(i)%sz, bxasc%r_con%rtr(i)%pr, tag)
       end do
       do i = 1, bxasc%r_con%nsp
          rst(i+bxasc%r_con%nrp) =  parallel_isend_lv(g_snd_l(1+nc*bxasc%r_con%str(i)%pv), &
               nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
       end do
       call parallel_wait(rst)
    end if

    !$OMP PARALLEL DO PRIVATE(i,sh,p)
    do i = 1, bxasc%r_con%nrcv
       sh = bxasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mf, bxasc%r_con%rcv(i)%lnd, bxasc%r_con%rcv(i)%dbx, c, nc)
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
    integer                    :: tag
    integer                    :: i, ii, jj, sh(MAX_SPACEDIM+1)
    type(boxassoc)             :: bxasc
    complex(dp_t), allocatable :: g_snd_z(:), g_rcv_z(:)

    bxasc = layout_boxassoc(mf%la, ng, mf%nodal, lcross)

    do i = 1, bxasc%l_con%ncpy
       ii  =  bxasc%l_con%cpy(i)%lnd
       jj  =  bxasc%l_con%cpy(i)%lns
       p1  => dataptr(mf%fbs(ii), bxasc%l_con%cpy(i)%dbx, c, nc)
       p2  => dataptr(mf%fbs(jj), bxasc%l_con%cpy(i)%sbx, c, nc)
       call cpy_z(p1,p2)
    end do

    allocate(g_snd_z(nc*bxasc%r_con%svol))
    allocate(g_rcv_z(nc*bxasc%r_con%rvol))

    do i = 1, bxasc%r_con%nsnd
       p => dataptr(mf, bxasc%r_con%snd(i)%lns, bxasc%r_con%snd(i)%sbx, c, nc)
       call reshape_z_4_1(g_snd_z, 1 + nc*bxasc%r_con%snd(i)%pv, p)
    end do

    tag = parallel_tag()

    if (bxasc%r_con%nrp+bxasc%r_con%nsp .gt. 0) then
       allocate(rst(bxasc%r_con%nrp+bxasc%r_con%nsp))
       do i = 1, bxasc%r_con%nrp
          rst(i) = parallel_irecv_zv(g_rcv_z(1+nc*bxasc%r_con%rtr(i)%pv), &
               nc*bxasc%r_con%rtr(i)%sz, bxasc%r_con%rtr(i)%pr, tag)
       end do
       do i = 1, bxasc%r_con%nsp
          rst(i+bxasc%r_con%nrp) = parallel_isend_zv(g_snd_z(1+nc*bxasc%r_con%str(i)%pv), &
               nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
       end do
       call parallel_wait(rst)
    end if

    do i = 1, bxasc%r_con%nrcv
       sh = bxasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mf, bxasc%r_con%rcv(i)%lnd, bxasc%r_con%rcv(i)%dbx, c, nc)
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

    ! multifab_fill_ghost_cells calls fillpatch with a multifab built
    ! from a boxarray of ghost boxes of the fine multifab to be filled. 
    ! And the multifab has ng=0.  So there are "valid" cells that are outside the domain.
    ! For non-periodic boundaries, these cells outside the domain will be filled
    ! by multifab_physbc later in multifab_fill_ghost_cells.  However, for
    ! periodic boundaries, mf_fb_fancy_double must be called for these to be filled.
    if ( lng < 1 .and. .not.any(mf%la%lap%pmask) ) return

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

    ! multifab_fill_ghost_cells calls fillpatch with a multifab built
    ! from a boxarray of ghost boxes of the fine multifab to be filled. 
    ! And the multifab has ng=0.  So there are "valid" cells that are outside the domain.
    ! For non-periodic boundaries, these cells outside the domain will be filled
    ! by multifab_physbc later in multifab_fill_ghost_cells.  However, for
    ! periodic boundaries, mf_fb_fancy_double must be called for these to be filled.
    if ( lng < 1 .and. .not.any(mf%la%lap%pmask) ) return

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

    ! multifab_fill_ghost_cells calls fillpatch with a multifab built
    ! from a boxarray of ghost boxes of the fine multifab to be filled. 
    ! And the multifab has ng=0.  So there are "valid" cells that are outside the domain.
    ! For non-periodic boundaries, these cells outside the domain will be filled
    ! by multifab_physbc later in multifab_fill_ghost_cells.  However, for
    ! periodic boundaries, mf_fb_fancy_double must be called for these to be filled.
    if ( lng < 1 .and. .not.any(mf%la%lap%pmask) ) return

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
    integer                 :: tag
    integer                 :: i, ii, jj, np, sh(MAX_SPACEDIM+1), lng
    type(boxassoc)          :: bxasc
    real(dp_t), allocatable :: g_snd_d(:), g_rcv_d(:)

    lng = mf%ng; if ( present(ng) ) lng = ng

    if ( mf%nc < (c+nc-1) ) call bl_error('MULTIFAB_SUM_BOUNDARY_C: nc too large', nc)
    if ( lng > mf%ng      ) call bl_error('MULTIFAB_SUM_BOUNDARY_C: ng too large', lng)
    if ( lng < 1          ) return

    call boxassoc_build(bxasc, mf%la%lap, lng, mf%nodal, .false., .true.)

    ! unsafe to do OMP
    do i = 1, bxasc%l_con%ncpy
       ii   =  bxasc%l_con%cpy(i)%lnd
       jj   =  bxasc%l_con%cpy(i)%lns
       pdst => dataptr(mf%fbs(ii), bxasc%l_con%cpy(i)%dbx, c, nc)
       psrc => dataptr(mf%fbs(jj), bxasc%l_con%cpy(i)%sbx, c, nc)
       call sum_d(pdst,psrc)
    end do

    np = parallel_nprocs()
    
    if (np == 1) then
       call boxassoc_destroy(bxasc)
       return
    end if

    tag = parallel_tag()

    allocate(g_snd_d(nc*bxasc%r_con%svol))
    allocate(g_rcv_d(nc*bxasc%r_con%rvol))

    !$OMP PARALLEL DO PRIVATE(i,p)
    do i = 1, bxasc%r_con%nsnd
       p => dataptr(mf, bxasc%r_con%snd(i)%lns, bxasc%r_con%snd(i)%sbx, c, nc)
       call reshape_d_4_1(g_snd_d, 1 + nc*bxasc%r_con%snd(i)%pv, p)
    end do
    !$OMP END PARALLEL DO

    if (bxasc%r_con%nrp+bxasc%r_con%nsp .gt. 0) then
       allocate(rst(bxasc%r_con%nrp+bxasc%r_con%nsp))
       do i = 1, bxasc%r_con%nrp
          rst(i) = parallel_irecv_dv(g_rcv_d(1+nc*bxasc%r_con%rtr(i)%pv), &
               nc*bxasc%r_con%rtr(i)%sz, bxasc%r_con%rtr(i)%pr, tag)
       end do
       do i = 1, bxasc%r_con%nsp
          rst(i+bxasc%r_con%nrp) = parallel_isend_dv(g_snd_d(1+nc*bxasc%r_con%str(i)%pv), &
               nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
       end do
       call parallel_wait(rst)
    end if

    ! unsafe to do OMP
    do i = 1, bxasc%r_con%nrcv
       sh = bxasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mf, bxasc%r_con%rcv(i)%lnd, bxasc%r_con%rcv(i)%dbx, c, nc)
       call reshape_d_1_4(p, g_rcv_d, 1 + nc*bxasc%r_con%rcv(i)%pv, sh, sum_d)
    end do

    call boxassoc_destroy(bxasc)

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
    integer              :: tag
    integer              :: i, ii, jj, np, sh(MAX_SPACEDIM+1), lng
    type(boxassoc)       :: bxasc
    logical, allocatable :: g_snd_l(:), g_rcv_l(:)

    lng = mf%ng; if ( present(ng) ) lng = ng

    if ( mf%nc < (c+nc-1) ) call bl_error('LMULTIFAB_SUM_BOUNDARY_C: nc too large', nc)
    if ( lng > mf%ng      ) call bl_error('LMULTIFAB_SUM_BOUNDARY_C: ng too large', lng)
    if ( lng < 1          ) return

    call boxassoc_build(bxasc, mf%la%lap, lng, mf%nodal, .false., .true.)

    do i = 1, bxasc%l_con%ncpy
       ii   =  bxasc%l_con%cpy(i)%lnd
       jj   =  bxasc%l_con%cpy(i)%lns
       pdst => dataptr(mf%fbs(ii), bxasc%l_con%cpy(i)%dbx, c, nc)
       psrc => dataptr(mf%fbs(jj), bxasc%l_con%cpy(i)%sbx, c, nc)
       call logical_or(pdst,psrc)
    end do

    np = parallel_nprocs()
    
    if (np == 1) then
       call boxassoc_destroy(bxasc)
       return
    end if

    tag = parallel_tag()

    allocate(g_snd_l(nc*bxasc%r_con%svol))
    allocate(g_rcv_l(nc*bxasc%r_con%rvol))

    do i = 1, bxasc%r_con%nsnd
       p => dataptr(mf, bxasc%r_con%snd(i)%lns, bxasc%r_con%snd(i)%sbx, c, nc)
       call reshape_l_4_1(g_snd_l, 1 + nc*bxasc%r_con%snd(i)%pv, p)
    end do

    if (bxasc%r_con%nrp+bxasc%r_con%nsp .gt. 0) then
       allocate(rst(bxasc%r_con%nrp+bxasc%r_con%nsp))
       do i = 1, bxasc%r_con%nrp
          rst(i) = parallel_irecv_lv(g_rcv_l(1+nc*bxasc%r_con%rtr(i)%pv), &
               nc*bxasc%r_con%rtr(i)%sz, bxasc%r_con%rtr(i)%pr, tag)
       end do
       do i = 1, bxasc%r_con%nsp
          rst(i+bxasc%r_con%nrp) = parallel_isend_lv(g_snd_l(1+nc*bxasc%r_con%str(i)%pv), &
               nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
       end do
       call parallel_wait(rst)
    end if

    do i = 1, bxasc%r_con%nrcv
       sh = bxasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mf, bxasc%r_con%rcv(i)%lnd, bxasc%r_con%rcv(i)%dbx, c, nc)
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
    integer                                 :: tag
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
       ii   =  snasc%l_con%cpy(i)%lnd
       jj   =  snasc%l_con%cpy(i)%lns
       pdst => dataptr(mf%fbs(ii), snasc%l_con%cpy(i)%dbx, c, nc)
       psrc => dataptr(mf%fbs(jj), snasc%l_con%cpy(i)%sbx, c, nc)
       call cpy_d(pdst, psrc, filter)
    end do

    np = parallel_nprocs()

    if (np == 1) return

    tag = parallel_tag()

    allocate(g_snd_d(nc*snasc%r_con%svol))
    allocate(g_rcv_d(nc*snasc%r_con%rvol))

    do i = 1, snasc%r_con%nsnd
       p => dataptr(mf, snasc%r_con%snd(i)%lns, snasc%r_con%snd(i)%sbx, c, nc)
       call reshape_d_4_1(g_snd_d, 1 + nc*snasc%r_con%snd(i)%pv, p)
    end do

    if (snasc%r_con%nrp+snasc%r_con%nsp .gt. 0) then
       allocate(rst(snasc%r_con%nrp+snasc%r_con%nsp))
       do i = 1, snasc%r_con%nrp
          rst(i) = parallel_irecv_dv(g_rcv_d(1+nc*snasc%r_con%rtr(i)%pv), &
               nc*snasc%r_con%rtr(i)%sz, snasc%r_con%rtr(i)%pr, tag)
       end do
       do i = 1, snasc%r_con%nsp
          rst(i+snasc%r_con%nrp) = parallel_isend_dv(g_snd_d(1+nc*snasc%r_con%str(i)%pv), &
               nc*snasc%r_con%str(i)%sz, snasc%r_con%str(i)%pr, tag)
       end do
       call parallel_wait(rst)
    end if

    do i = 1, snasc%r_con%nrcv
       sh = snasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mf, snasc%r_con%rcv(i)%lnd, snasc%r_con%rcv(i)%dbx, c, nc)
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
    integer                              :: tag
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
       ii   =  snasc%l_con%cpy(i)%lnd
       jj   =  snasc%l_con%cpy(i)%lns
       pdst => dataptr(mf%fbs(ii), snasc%l_con%cpy(i)%dbx, c, nc)
       psrc => dataptr(mf%fbs(jj), snasc%l_con%cpy(i)%sbx, c, nc)
       call cpy_l(pdst, psrc, filter)
    end do

    np = parallel_nprocs()

    if (np == 1) return

    tag = parallel_tag()

    allocate(g_snd_l(nc*snasc%r_con%svol))
    allocate(g_rcv_l(nc*snasc%r_con%rvol))

    do i = 1, snasc%r_con%nsnd
       p => dataptr(mf, snasc%r_con%snd(i)%lns, snasc%r_con%snd(i)%sbx, c, nc)
       call reshape_l_4_1(g_snd_l, 1 + nc*snasc%r_con%snd(i)%pv, p)
    end do

    if (snasc%r_con%nrp+snasc%r_con%nsp .gt. 0) then
       allocate(rst(snasc%r_con%nrp+snasc%r_con%nsp))
       do i = 1, snasc%r_con%nrp
          rst(i) = parallel_irecv_lv(g_rcv_l(1+nc*snasc%r_con%rtr(i)%pv), &
               nc*snasc%r_con%rtr(i)%sz, snasc%r_con%rtr(i)%pr, tag)
       end do
       do i = 1, snasc%r_con%nsp
          rst(i+snasc%r_con%nrp) = parallel_isend_lv(g_snd_l(1+nc*snasc%r_con%str(i)%pv), &
               nc*snasc%r_con%str(i)%sz, snasc%r_con%str(i)%pr, tag)
       end do
       call parallel_wait(rst)
    end if

    do i = 1, snasc%r_con%nrcv
       sh = snasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mf, snasc%r_con%rcv(i)%lnd, snasc%r_con%rcv(i)%dbx, c, nc)
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
    integer                   :: tag
    type(bl_prof_timer), save :: bpt

    lnc = 1; if ( present(nc)  ) lnc  = nc

    if ( mf_in%nc  < (c_in +lnc-1) ) call bl_error('MULTIFAB_COPY_ON_SHIFT: nc too large', lnc)
    if ( mf_out%nc < (c_out+lnc-1) ) call bl_error('MULTIFAB_COPY_ON_SHIFT: nc too large', lnc)

    call build(bpt, "mf_copy_on_shift")

    tag = parallel_tag()

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

  subroutine mf_copy_fancy_double(mdst, dstcomp, msrc, srccomp, nc, filter, bndry_reg_to_other)
    type(multifab), intent(inout) :: mdst
    type(multifab), intent(in)    :: msrc
    integer, intent(in)           :: dstcomp, srccomp, nc
    logical, intent(in), optional :: bndry_reg_to_other

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
    integer                 :: tag
    integer                 :: i, ii, jj, sh(MAX_SPACEDIM+1), np
    real(dp_t), allocatable :: g_snd_d(:), g_rcv_d(:)
    logical                 :: br_to_other
    
    br_to_other = .false.  
    if (present(bndry_reg_to_other)) br_to_other = bndry_reg_to_other

    if (br_to_other) then
       call copyassoc_build_br_to_other(cpasc, mdst%la, msrc%la, mdst%nodal, msrc%nodal)
    else
       cpasc = layout_copyassoc(mdst%la, msrc%la, mdst%nodal, msrc%nodal)
    end if

    !$OMP PARALLEL DO PRIVATE(i,ii,jj,pdst,psrc) if (cpasc%l_con%threadsafe)
    do i = 1, cpasc%l_con%ncpy
       ii   =  cpasc%l_con%cpy(i)%lnd
       jj   =  cpasc%l_con%cpy(i)%lns
       pdst => dataptr(mdst%fbs(ii), cpasc%l_con%cpy(i)%dbx, dstcomp, nc)
       psrc => dataptr(msrc%fbs(jj), cpasc%l_con%cpy(i)%sbx, srccomp, nc)
       call cpy_d(pdst, psrc, filter)
    end do
    !$OMP END PARALLEL DO

    np = parallel_nprocs()

    if (np == 1) then
       if (br_to_other) call copyassoc_destroy(cpasc)
       return
    end if

    tag = parallel_tag()

    allocate(g_snd_d(nc*cpasc%r_con%svol))
    allocate(g_rcv_d(nc*cpasc%r_con%rvol))

    !$OMP PARALLEL DO PRIVATE(i,p)
    do i = 1, cpasc%r_con%nsnd
       p => dataptr(msrc, cpasc%r_con%snd(i)%lns, cpasc%r_con%snd(i)%sbx, srccomp, nc)
       call reshape_d_4_1(g_snd_d, 1 + nc*cpasc%r_con%snd(i)%pv, p)
    end do
    !$OMP END PARALLEL DO

    if (cpasc%r_con%nrp+cpasc%r_con%nsp .gt. 0) then
       allocate(rst(cpasc%r_con%nrp+cpasc%r_con%nsp))
       do i = 1, cpasc%r_con%nrp
          rst(i) = parallel_irecv_dv(g_rcv_d(1+nc*cpasc%r_con%rtr(i)%pv), &
               nc*cpasc%r_con%rtr(i)%sz, cpasc%r_con%rtr(i)%pr, tag)
       end do
       do i = 1, cpasc%r_con%nsp
          rst(i+cpasc%r_con%nrp) = parallel_isend_dv(g_snd_d(1+nc*cpasc%r_con%str(i)%pv), &
               nc*cpasc%r_con%str(i)%sz, cpasc%r_con%str(i)%pr, tag)
       end do
       call parallel_wait(rst)
    end if

    !$OMP PARALLEL DO PRIVATE(i,sh,p) if (cpasc%r_con%threadsafe)
    do i = 1, cpasc%r_con%nrcv
       sh = cpasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mdst, cpasc%r_con%rcv(i)%lnd, cpasc%r_con%rcv(i)%dbx, dstcomp, nc)
       call reshape_d_1_4(p, g_rcv_d, 1 + nc*cpasc%r_con%rcv(i)%pv, sh, filter)
    end do
    !$OMP END PARALLEL DO    

    if (br_to_other) call copyassoc_destroy(cpasc)

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
    integer              :: tag
    integer              :: i, ii, jj, np, sh(MAX_SPACEDIM+1)
    integer, allocatable :: g_snd_i(:), g_rcv_i(:)

    cpasc = layout_copyassoc(mdst%la, msrc%la, mdst%nodal, msrc%nodal)

    do i = 1, cpasc%l_con%ncpy
       ii   =  cpasc%l_con%cpy(i)%lnd
       jj   =  cpasc%l_con%cpy(i)%lns
       pdst => dataptr(mdst%fbs(ii), cpasc%l_con%cpy(i)%dbx, dstcomp, nc)
       psrc => dataptr(msrc%fbs(jj), cpasc%l_con%cpy(i)%sbx, srccomp, nc)
       call cpy_i(pdst, psrc, filter)
    end do

    np = parallel_nprocs()

    if (np == 1) return

    tag = parallel_tag()

    allocate(g_snd_i(nc*cpasc%r_con%svol))
    allocate(g_rcv_i(nc*cpasc%r_con%rvol))

    do i = 1, cpasc%r_con%nsnd
       p => dataptr(msrc, cpasc%r_con%snd(i)%lns, cpasc%r_con%snd(i)%sbx, srccomp, nc)
       call reshape_i_4_1(g_snd_i, 1 + nc*cpasc%r_con%snd(i)%pv, p)
    end do

    if (cpasc%r_con%nrp+cpasc%r_con%nsp .gt. 0) then
       allocate(rst(cpasc%r_con%nrp+cpasc%r_con%nsp))
       do i = 1, cpasc%r_con%nrp
          rst(i) = parallel_irecv_iv(g_rcv_i(1+nc*cpasc%r_con%rtr(i)%pv), &
               nc*cpasc%r_con%rtr(i)%sz, cpasc%r_con%rtr(i)%pr, tag)
       end do
       do i = 1, cpasc%r_con%nsp
          rst(i+cpasc%r_con%nrp) = parallel_isend_iv(g_snd_i(1+nc*cpasc%r_con%str(i)%pv), &
               nc*cpasc%r_con%str(i)%sz, cpasc%r_con%str(i)%pr, tag)
       end do
       call parallel_wait(rst)
    end if

    do i = 1, cpasc%r_con%nrcv
       sh = cpasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mdst, cpasc%r_con%rcv(i)%lnd, cpasc%r_con%rcv(i)%dbx, dstcomp, nc)
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
    integer              :: tag 
    integer              :: i, ii, jj, np, sh(MAX_SPACEDIM+1)
    logical, allocatable :: g_snd_l(:), g_rcv_l(:)

    cpasc = layout_copyassoc(mdst%la, msrc%la, mdst%nodal, msrc%nodal)

    do i = 1, cpasc%l_con%ncpy
       ii   =  cpasc%l_con%cpy(i)%lnd
       jj   =  cpasc%l_con%cpy(i)%lns
       pdst => dataptr(mdst%fbs(ii), cpasc%l_con%cpy(i)%dbx, dstcomp, nc)
       psrc => dataptr(msrc%fbs(jj), cpasc%l_con%cpy(i)%sbx, srccomp, nc)
       call cpy_l(pdst, psrc, filter)
    end do

    np = parallel_nprocs()

    if (np == 1) return

    tag = parallel_tag()

    allocate(g_snd_l(nc*cpasc%r_con%svol))
    allocate(g_rcv_l(nc*cpasc%r_con%rvol))

    do i = 1, cpasc%r_con%nsnd
       p => dataptr(msrc, cpasc%r_con%snd(i)%lns, cpasc%r_con%snd(i)%sbx, srccomp, nc)
       call reshape_l_4_1(g_snd_l, 1 + nc*cpasc%r_con%snd(i)%pv, p)
    end do

    if (cpasc%r_con%nrp+cpasc%r_con%nsp .gt. 0) then
       allocate(rst(cpasc%r_con%nrp+cpasc%r_con%nsp))
       do i = 1, cpasc%r_con%nrp
          rst(i) = parallel_irecv_lv(g_rcv_l(1+nc*cpasc%r_con%rtr(i)%pv), &
               nc*cpasc%r_con%rtr(i)%sz, cpasc%r_con%rtr(i)%pr, tag)
       end do
       do i = 1, cpasc%r_con%nsp
          rst(i+cpasc%r_con%nrp) = parallel_isend_lv(g_snd_l(1+nc*cpasc%r_con%str(i)%pv), &
               nc*cpasc%r_con%str(i)%sz, cpasc%r_con%str(i)%pr, tag)
       end do
       call parallel_wait(rst)
    end if

    do i = 1, cpasc%r_con%nrcv
       sh = cpasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mdst, cpasc%r_con%rcv(i)%lnd, cpasc%r_con%rcv(i)%dbx, dstcomp, nc)
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
    integer                    :: tag
    integer                    :: i, ii, jj, sh(MAX_SPACEDIM+1)
    complex(dp_t), allocatable :: g_snd_z(:), g_rcv_z(:)

    cpasc = layout_copyassoc(mdst%la, msrc%la, mdst%nodal, msrc%nodal)

    do i = 1, cpasc%l_con%ncpy
       ii   =  cpasc%l_con%cpy(i)%lnd
       jj   =  cpasc%l_con%cpy(i)%lns
       pdst => dataptr(mdst%fbs(ii), cpasc%l_con%cpy(i)%dbx, dstcomp, nc)
       psrc => dataptr(msrc%fbs(jj), cpasc%l_con%cpy(i)%sbx, srccomp, nc)
       call cpy_z(pdst, psrc, filter)
    end do

    if (parallel_nprocs() == 1) return

    tag = parallel_tag()

    allocate(g_snd_z(nc*cpasc%r_con%svol))
    allocate(g_rcv_z(nc*cpasc%r_con%rvol))

    do i = 1, cpasc%r_con%nsnd
       p => dataptr(msrc, cpasc%r_con%snd(i)%lns, cpasc%r_con%snd(i)%sbx, srccomp, nc)
       call reshape_z_4_1(g_snd_z, 1 + nc*cpasc%r_con%snd(i)%pv, p)
    end do

    if (cpasc%r_con%nrp+cpasc%r_con%nsp .gt. 0) then
       allocate(rst(cpasc%r_con%nrp+cpasc%r_con%nsp))
       do i = 1, cpasc%r_con%nrp
          rst(i) = parallel_irecv_zv(g_rcv_z(1+nc*cpasc%r_con%rtr(i)%pv), &
               nc*cpasc%r_con%rtr(i)%sz, cpasc%r_con%rtr(i)%pr, tag)
       end do
       do i = 1, cpasc%r_con%nsp
          rst(i+cpasc%r_con%nrp) = parallel_isend_zv(g_snd_z(1+nc*cpasc%r_con%str(i)%pv), &
               nc*cpasc%r_con%str(i)%sz, cpasc%r_con%str(i)%pr, tag)
       end do
       call parallel_wait(rst)
    end if

    do i = 1, cpasc%r_con%nrcv
       sh = cpasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mdst, cpasc%r_con%rcv(i)%lnd, cpasc%r_con%rcv(i)%dbx, dstcomp, nc)
       call reshape_z_1_4(p, g_rcv_z, 1 + nc*cpasc%r_con%rcv(i)%pv, sh, filter)
    end do

  end subroutine mf_copy_fancy_z

  subroutine multifab_copy_c(mdst, dstcomp, msrc, srccomp, nc, ng, filter, bndry_reg_to_other, ngsrc)
    type(multifab), intent(inout)        :: mdst
    type(multifab), intent(in)   ,target :: msrc
    integer, intent(in)           :: dstcomp, srccomp
    integer, intent(in), optional :: nc
    integer, intent(in), optional :: ng, ngsrc
    logical, intent(in), optional :: bndry_reg_to_other
    real(dp_t), pointer           :: pdst(:,:,:,:), psrc(:,:,:,:)
    type(multifab), pointer       :: pmfsrc
    type(multifab), target        :: msrctmp
    type(layout)                  :: lasrctmp
    type(boxarray)                :: batmp
    type(list_box)                :: bl
    integer                       :: i, lnc, lng, lngsrc
    type(mfiter)                  :: mfi
    type(box)                     :: bx
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
    lngsrc  = 0;       if ( present(ngsrc)  ) lngsrc = ngsrc
    if ( lnc < 1 )                   call bl_error('MULTIFAB_COPY_C: nc must be >= 1')
    if ( mdst%nc < (dstcomp+lnc-1) ) call bl_error('MULTIFAB_COPY_C: nc too large for dst multifab', lnc)
    if ( msrc%nc < (srccomp+lnc-1) ) call bl_error('MULTIFAB_COPY_C: nc too large for src multifab', lnc)
    if ( mdst%la == msrc%la ) then
       if ( lng > 0 ) &
            call bl_assert(mdst%ng >= lng, msrc%ng >= lng,"not enough ghost cells in multifab_copy_c")

       !$omp parallel private(mfi,i,bx,pdst,psrc)
       call mfiter_build(mfi,mdst,.true.)
       do while (next_tile(mfi,i))
          bx = get_growntilebox(mfi,lng)
          pdst => dataptr(mdst, i, bx, dstcomp, lnc)
          psrc => dataptr(msrc, i, bx, srccomp, lnc)
          call cpy_d(pdst, psrc, filter)
       end do
       !$omp end parallel
    else
       if (lng    >       0) call bl_error('MULTIFAB_COPY_C: ng > 0 not supported in parallel copy')
       if (lngsrc > msrc%ng) call bl_error('MULTIFAB_COPY_C: ngsrc > msrc%ng')

       if (lngsrc > 0) then
          do i = 1, nboxes(msrc%la)
             call push_back(bl, grow(box_nodalize(get_box(msrc%la,i),msrc%nodal),lngsrc))
          end do
          call boxarray_build_l(batmp, bl, sort = .false.)
          call destroy(bl)
          call layout_build_ba(lasrctmp, batmp, boxarray_bbox(batmp), explicit_mapping = get_proc(msrc%la))
          call boxarray_destroy(batmp)
          call multifab_build(msrctmp, lasrctmp, nc = msrc%nc, ng = msrc%ng-lngsrc, fab_alloc=.false.)
          msrctmp%fbs => msrc%fbs

          pmfsrc => msrctmp
       else
          pmfsrc => msrc
       end if

       call mf_copy_fancy_double(mdst, dstcomp, pmfsrc, srccomp, lnc, filter, bndry_reg_to_other)

       if (lngsrc > 0) then
          nullify(msrctmp%fbs)
          call multifab_destroy(msrctmp)
          call layout_destroy(lasrctmp)
       end if
    end if
    call destroy(bpt)
  end subroutine multifab_copy_c

  subroutine multifab_copy(mdst, msrc, ng, filter, bndry_reg_to_other, ngsrc)
    type(multifab), intent(inout) :: mdst
    type(multifab), intent(in)    :: msrc
    integer, intent(in), optional :: ng, ngsrc
    logical, intent(in), optional :: bndry_reg_to_other
    interface
       subroutine filter(out, in)
         use bl_types
         real(dp_t), intent(inout) :: out(:,:,:,:)
         real(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    if ( mdst%nc .ne. msrc%nc ) call bl_error('MULTIFAB_COPY: multifabs must have same number of components')
    call multifab_copy_c(mdst, 1, msrc, 1, mdst%nc, ng, filter, bndry_reg_to_other, ngsrc)
  end subroutine multifab_copy

  subroutine imultifab_copy_c(mdst, dstcomp, msrc, srccomp, nc, ng, filter)
    type(imultifab), intent(inout) :: mdst
    type(imultifab), intent(in)    :: msrc
    integer, intent(in)            :: dstcomp, srccomp
    integer, intent(in), optional  :: nc
    integer, intent(in), optional  :: ng
    integer, pointer               :: pdst(:,:,:,:), psrc(:,:,:,:)
    integer                        :: i, lnc, lng
    type(mfiter)                   :: mfi
    type(box)                      :: bx
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
    if ( lng > 0 )                   call bl_assert(mdst%ng >= lng, msrc%ng >= lng,"not enough ghost cells in imultifab_copy_c")
    if ( mdst%la == msrc%la ) then
       !$omp parallel private(mfi,i,bx,pdst,psrc)
       call mfiter_build(mfi,mdst,.true.)
       do while(next_tile(mfi,i))
          bx = get_growntilebox(mfi, lng)
          pdst => dataptr(mdst, i, bx, dstcomp, lnc)
          psrc => dataptr(msrc, i, bx, srccomp, lnc)
          call cpy_i(pdst, psrc, filter)
       end do
       !$omp end parallel
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
    type(mfiter)                   :: mfi
    type(box)                      :: bx
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
    if ( lng > 0 )                   call bl_assert(mdst%ng >= lng, msrc%ng >= lng,"not enough ghost cells in lmultifab_copy_c")
    if ( mdst%la == msrc%la ) then
       !$omp parallel private(mfi,i,bx,pdst,psrc)
       call mfiter_build(mfi,mdst,.true.)
       do while (next_tile(mfi,i))
          bx = get_growntilebox(mfi,lng)
          pdst => dataptr(mdst, i, bx, dstcomp, lnc)
          psrc => dataptr(msrc, i, bx, srccomp, lnc)
          call cpy_l(pdst, psrc, filter)
       end do
       !$omp end parallel
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
    type(mfiter)                   :: mfi
    type(box)                      :: bx
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
    if ( lng > 0 )                   call bl_assert(mdst%ng >= lng, msrc%ng >= lng,"not enough ghost cells in zmultifab_copy_c")
    if ( mdst%la == msrc%la ) then
       !$omp parallel private(mfi,i,bx,pdst,psrc)
       call mfiter_build(mfi,mdst,.true.)
       do while (next_tile(mfi,i))
          bx = get_growntilebox(mfi,lng)
          pdst => dataptr(mdst, i, bx, dstcomp, lnc)
          psrc => dataptr(msrc, i, bx, srccomp, lnc)
          call cpy_z(pdst, psrc, filter)
       end do
       !$omp end parallel
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

    call multifab_build(mask, mf%la, 1, 0, mf%nodal)

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
    use omp_module
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
    integer             :: i,j,k,n,lo(4),hi(4)
    logical             :: llocal
    type(mfiter)        :: mfi
    type(box)           :: bx
    integer             :: tid, nthreads
    real(dp_t)          :: r1, r2
    real(dp_t), allocatable :: rt(:,:)

    if ( present(mask) ) then
       if ( ncomp(mask) /= 1 ) call bl_error('Mask array is multicomponent')
    end if
    llocal = .false.; if ( present(local) ) llocal = local

    ! a bit of hack to get consistent answer with OMP
    nthreads = omp_get_max_threads()
    allocate(rt(16,0:nthreads-1)) ! extra padding to avoid false sharing

    rt = 0.0_dp_t

    if ( cell_centered_q(mf) ) then
       !$omp parallel private(mp,mp1,lmp,i,j,k,n,lo,hi,mfi,bx,tid,r2)
       tid = omp_get_thread_num()
       call mfiter_build(mfi,mf,.true.)
       do while(next_tile(mfi,n))
          bx = get_tilebox(mfi)

          mp  => dataptr(mf %fbs(n), bx, comp )
          mp1 => dataptr(mf1%fbs(n), bx, comp1)

          lo = lbound(mp); hi = ubound(mp)

          r2 = 0.0_dp_t

          if ( present(mask) )then
             lmp => dataptr(mask%fbs(n), bx)
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

          rt(1,tid) = rt(1,tid) + r2

       end do
       !$omp end parallel
    else if ( nodal_q(mf) ) then

       if ( .not. present(nodal_mask) ) call build_nodal_dot_mask(tmask, mf)

       !$omp parallel private(mp,mp1,ma,lmp,i,j,k,n,lo,hi,mfi,bx,tid,r2)
       tid = omp_get_thread_num()
       call mfiter_build(mfi,mf,.true.)
       do while(next_tile(mfi,n))
          bx = get_tilebox(mfi)

          mp  => dataptr(mf %fbs(n), bx, comp )
          mp1 => dataptr(mf1%fbs(n), bx, comp1)
          if ( present(nodal_mask) ) then
             ma => dataptr(nodal_mask%fbs(n), bx)
          else
             ma => dataptr(     tmask%fbs(n), bx)
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

          rt(1,tid) = rt(1,tid) + r2

       end do
       !$omp end parallel

       if ( .not. present(nodal_mask) ) call multifab_destroy(tmask)
    else
       call bl_error("MULTIFAB_DOT_CC, fails when not nodal or cell-centered, can be fixed")
    end if

    r1 = sum(rt(1,:))
    r = r1

    if ( .not. llocal ) then
       if ( present(comm) ) then
          call parallel_reduce(r, r1, MPI_SUM, comm = comm)
       else
          call parallel_reduce(r, r1, MPI_SUM)
       end if
    end if

  end function multifab_dot_cc

  function multifab_dot_c(mf, mf1, comp, nodal_mask, local, comm) result(r)
    real(dp_t)                 :: r
    type(multifab), intent(in) :: mf
    type(multifab), intent(in) :: mf1
    integer       , intent(in) :: comp
    type(multifab), intent(in), optional :: nodal_mask
    logical       , intent(in), optional :: local
    integer       , intent(in), optional :: comm
    r = multifab_dot_cc(mf, comp, mf1, comp, nodal_mask = nodal_mask, local = local, comm = comm)
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
    type(mfiter) :: mfi
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
    !$omp parallel private(i,mp,mfi)
    call mfiter_build(mfi,mf,.true.)
    do while(next_tile(mfi,i))
       mp => dataptr(mf%fbs(i), get_tilebox(mfi), c)
       if (lclip) then
          where ( mp < lmin )
             mp = lxmin
          elsewhere ( mp > lmax )
             mp = lxmax
          elsewhere
             mp = lxmax*(mp-lmin)/(lmax-lmin) + lxmin
          end where
       else
          mp = lxmax*(mp-lmin)/(lmax-lmin) + lxmin
       end if
    end do
    !$omp end parallel
  end subroutine multifab_rescale_2

  subroutine multifab_rescale_c(mf, c, val, off)
    real(dp_t), intent(in) :: val
    integer, intent(in) :: c
    real(dp_t), intent(in), optional :: off
    type(multifab), intent(inout) :: mf
    real(dp_t), pointer :: mp(:,:,:,:)
    integer :: i
    type(mfiter) :: mfi
    !$omp parallel private(mp,i,mfi)
    call mfiter_build(mfi,mf,.true.)
    do while(next_tile(mfi,i))
       mp => dataptr(mf%fbs(i), get_tilebox(mfi), c)
       if ( present(off) ) then
          mp = mp*val + off
       else
          mp = mp*val
       end if
    end do
    !$omp end parallel
  end subroutine multifab_rescale_c
  subroutine multifab_rescale(mf, val, off)
    real(dp_t), intent(in) :: val
    real(dp_t), intent(in), optional :: off
    type(multifab), intent(inout) :: mf
    real(dp_t), pointer :: mp(:,:,:,:)
    integer :: i
    type(mfiter) :: mfi
    !$omp parallel private(mp,i,mfi)
    call mfiter_build(mfi,mf,.true.)
    do while(next_tile(mfi,i))
       mp => dataptr(mf%fbs(i), get_tilebox(mfi))
       if ( present(off) ) then
          mp = mp*val + off
       else
          mp = mp*val
       end if
    end do
    !$omp end parallel
  end subroutine multifab_rescale

  subroutine multifab_saxpy_5(a, b1, b, c1, c, all)
    real(dp_t), intent(in) :: b1, c1
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b, c
    logical, intent(in), optional :: all
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    real(dp_t), pointer :: cp(:,:,:,:)
    integer :: i
    logical :: lall
    type(mfiter) :: mfi
    type(box) :: bx
    lall = .false.; if ( present(all) ) lall = all
    !$omp parallel private(ap,bp,cp,i,mfi,bx)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       if ( lall ) then
          bx = get_growntilebox(mfi)
       else
          bx = get_tilebox(mfi)
       end if
       ap => dataptr(a%fbs(i), bx)
       bp => dataptr(b%fbs(i), bx)
       cp => dataptr(c%fbs(i), bx)
       ap = b1*bp + c1*cp
    end do
    !$omp end parallel
  end subroutine multifab_saxpy_5

  subroutine multifab_saxpy_4(a, b, c1, c, all)
    real(dp_t),     intent(in)    :: c1
    type(multifab), intent(inout) :: a
    type(multifab), intent(in   ) :: b,c
    logical, intent(in), optional :: all
    real(dp_t), pointer           :: ap(:,:,:,:)
    real(dp_t), pointer           :: bp(:,:,:,:)
    real(dp_t), pointer           :: cp(:,:,:,:)
    
    integer :: ii, i, j, k, n, lo(4), hi(4)
    logical :: lall
    type(mfiter) :: mfi
    type(box) :: bx

    lall = .false.; if ( present(all) ) lall = all

    !$omp parallel private(ap,bp,cp,ii,i,j,k,n,lo,hi,mfi,bx)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,ii))
       if ( lall ) then
          bx = get_growntilebox(mfi)
       else
          bx = get_tilebox(mfi)
       end if

       ap => dataptr(a%fbs(ii), bx)
       bp => dataptr(b%fbs(ii), bx)
       cp => dataptr(c%fbs(ii), bx)

       lo = lbound(ap)
       hi = ubound(ap)
       
       ! ap = bp + c1*cp

       do n = lo(4), hi(4)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   ap(i,j,k,n) = bp(i,j,k,n) + c1 * cp(i,j,k,n)
                end do
             end do
          end do
       end do
    end do
    !$omp end parallel
  end subroutine multifab_saxpy_4

  subroutine multifab_saxpy_3_doit(ap, b1, bp)

    real(dp_t),     intent(in   ) :: b1
    real(dp_t), pointer           :: ap(:,:,:,:)
    real(dp_t), pointer           :: bp(:,:,:,:)
    integer                       :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! ap = ap + b1*bp

    do n = lo(4), hi(4)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) + b1 * bp(i,j,k,n)
             end do
          end do
       end do
    end do

  end subroutine multifab_saxpy_3_doit

  subroutine multifab_saxpy_3(a, b1, b, all)

    real(dp_t),     intent(in   ) :: b1
    type(multifab), intent(inout) :: a
    type(multifab), intent(in   ) :: b
    real(dp_t), pointer           :: ap(:,:,:,:)
    real(dp_t), pointer           :: bp(:,:,:,:)
    logical, intent(in), optional :: all
    integer :: i
    logical :: lall
    type(mfiter) :: mfi
    type(box) :: bx
    lall = .false.; if ( present(all) ) lall = all
    !$omp parallel private(ap,bp,i,mfi,bx)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       if ( lall ) then
          bx = get_growntilebox(mfi)
       else
          bx = get_tilebox(mfi)
       end if
       ap => dataptr(a%fbs(i), bx)
       bp => dataptr(b%fbs(i), bx)
       call multifab_saxpy_3_doit(ap,b1,bp)
    end do
    !$omp end parallel
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
    type(mfiter) :: mfi
    type(box) :: bx
    lall = .false.; if ( present(all) ) lall = all
    !$omp parallel private(ap,bp,i,mfi,bx)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       if ( lall ) then
          bx = get_growntilebox(mfi)
       else
          bx = get_tilebox(mfi)
       end if
       ap => dataptr(a%fbs(i), bx, ia)
       bp => dataptr(b%fbs(i), bx)
       call multifab_saxpy_3_doit(ap,b1,bp)
    end do
    !$omp end parallel
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
    integer :: i
    logical :: lall
    type(mfiter) :: mfi
    type(box) :: bx 
    lall = .false.; if ( present(all) ) lall = all
    !$omp parallel private(ap,bp,i,mfi,bx)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       if ( lall ) then
          bx = get_growntilebox(mfi)
       else
          bx = get_tilebox(mfi)
       end if
       ap => dataptr(a%fbs(i), bx, ia, nc)
       bp => dataptr(b%fbs(i), bx, ib, nc)
       call multifab_saxpy_3_doit(ap,b1,bp)
    end do
    !$omp end parallel
  end subroutine multifab_saxpy_3_cc

  function multifab_norm_l1_c(mf, comp, nc, mask, all) result(r)
    real(dp_t) :: r
    logical, intent(in), optional :: all
    integer, intent(in) :: comp
    integer, intent(in), optional :: nc
    type(multifab), intent(in) :: mf
    type(lmultifab), intent(in), optional :: mask
    logical, pointer :: lp(:,:,:,:)
    real(dp_t), pointer :: mp(:,:,:,:)
    integer :: i, n
    real(dp_t) :: r1
    logical :: lall
    integer :: lnc
    type(mfiter) :: mfi
    type(box) :: bx
    lnc  = 1; if ( present(nc) ) lnc = nc
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0.0_dp_t
    !$omp parallel private(mp,lp,i,n,mfi,bx) reduction(+:r1)
    call mfiter_build(mfi,mf,.true.)
    do while(next_tile(mfi,i))
       if ( lall ) then
          bx = get_growntilebox(mfi)
       else
          bx = get_tilebox(mfi)
       end if
       if ( present(mask) ) then
          lp => dataptr(mask%fbs(i), bx)
          do n = comp, comp + lnc - 1
             mp => dataptr(mf%fbs(i), bx, n)
             r1 = r1 + sum(abs(mp), mask = lp)
          end do
       else
          mp => dataptr(mf%fbs(i), bx, comp, lnc)
          r1 = r1 + sum(abs(mp))
       end if
    end do
    !$omp end parallel
    call parallel_reduce(r, r1, MPI_SUM)
  end function multifab_norm_l1_c
  function multifab_norm_l1(mf, all) result(r)
    real(dp_t)                    :: r
    type(multifab), intent(in)    :: mf
    logical, intent(in), optional :: all
    r = multifab_norm_l1_c(mf, 1, mf%nc, all = all)
  end function multifab_norm_l1

  function multifab_sum_c(mf, comp, nc, mask, all, local) result(r)
    real(dp_t) :: r
    integer, intent(in) :: comp
    logical, intent(in), optional :: all
    logical, intent(in), optional :: local
    integer, intent(in), optional :: nc
    type(multifab), intent(in) :: mf
    type(lmultifab), intent(in), optional :: mask
    real(dp_t), pointer :: mp(:,:,:,:)
    logical, pointer :: lp(:,:,:,:)
    integer :: i, n
    real(dp_t) :: r1
    logical :: lall, llocal
    integer :: lnc
    type(mfiter) :: mfi
    type(box) :: bx
    lnc  = 1; if ( present(nc) ) lnc = nc
    lall = .false.; if ( present(all) ) lall = all
    llocal = .false.; if ( present(local) ) llocal = local
    r1 = 0.0_dp_t
    !$omp parallel private(mp,lp,i,n,mfi,bx) reduction(+:r1)
    call mfiter_build(mfi,mf,.true.)
    do while(next_tile(mfi,i))
       if ( lall ) then
          bx = get_growntilebox(mfi)
       else
          bx = get_tilebox(mfi)
       end if
       if ( present(mask) ) then
          lp => dataptr(mask%fbs(i), bx)
          do n = comp, comp + lnc -1
             mp => dataptr(mf%fbs(i), bx, n)
             r1 = r1 + sum(mp, mask=lp)
          end do
       else
          mp => dataptr(mf%fbs(i), bx, comp, lnc)
          r1 = r1 + sum(mp)
       end if
    end do
    !$omp end parallel
    if (llocal) then
       r = r1
    else
       call parallel_reduce(r, r1, MPI_SUM)
    end if
  end function multifab_sum_c
  function multifab_sum(mf, mask, all, local) result(r)
    real(dp_t)                            :: r
    type(multifab), intent(in)            :: mf
    type(lmultifab), intent(in), optional :: mask
    logical, intent(in), optional         :: all
    logical, intent(in), optional         :: local
    r = multifab_sum_c(mf, 1, mf%nc, mask, all, local)
  end function multifab_sum

  function multifab_norm_l2_doit(ap, lp) result(r)

    real(dp_t), pointer        :: ap(:,:,:,:)
    logical, pointer, optional :: lp(:,:,:,:)
    real(dp_t)                 :: r,r1
    integer                    :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    r1 = 0.0_dp_t

    if ( present(lp) ) then
       do n = lo(4), hi(4)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   if (lp(i,j,k,n)) r1 = r1 + ap(i,j,k,n)*ap(i,j,k,n)
                end do
             end do
          end do
       end do
    else
       do n = lo(4), hi(4)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   r1 = r1 + ap(i,j,k,n)*ap(i,j,k,n)
                end do
             end do
          end do
       end do
    end if

    r = r1

  end function multifab_norm_l2_doit

  function multifab_norm_l2_c(mf, comp, nc, mask, all, local) result(r)
    real(dp_t) :: r
    integer, intent(in) :: comp
    logical, intent(in), optional :: all
    logical, intent(in), optional :: local
    integer, intent(in), optional :: nc
    type(multifab), intent(in) :: mf
    type(lmultifab), intent(in), optional :: mask
    real(dp_t), pointer :: mp(:,:,:,:)
    logical, pointer :: lp(:,:,:,:)
    integer :: i, n
    real(dp_t) :: r1
    logical :: lall, llocal
    integer :: lnc
    type(mfiter) :: mfi
    type(box) :: bx
    lnc  = 1; if ( present(nc) ) lnc = nc
    lall = .false.; if ( present(all) ) lall = all
    llocal = .false.; if ( present(local) ) llocal = local
    r1 = 0.0_dp_t
    !$omp parallel private(mp,lp,i,n,mfi,bx) reduction(+:r1)
    call mfiter_build(mfi,mf,.true.)
    do while(next_tile(mfi,i))
       if ( lall ) then
          bx = get_growntilebox(mfi)
       else
          bx = get_tilebox(mfi)
       end if
       if ( present(mask) ) then
          lp => dataptr(mask%fbs(i), bx)
          do n = comp, comp + lnc - 1
             mp => dataptr(mf%fbs(i), bx, n)
             r1 = r1 + multifab_norm_l2_doit(mp,lp)
          end do
       else
          mp => dataptr(mf%fbs(i), bx, comp, lnc)
          r1 = r1 + multifab_norm_l2_doit(mp)
       end if
    end do
    !$omp end parallel
    if (llocal) then
       r = r1
    else
       call parallel_reduce(r, r1, MPI_SUM)
    end if
    r = sqrt(r)
  end function multifab_norm_l2_c
  function multifab_norm_l2(mf, mask, all, local) result(r)
    real(dp_t)                            :: r
    logical, intent(in), optional         :: all
    logical, intent(in), optional         :: local
    type(multifab), intent(in)            :: mf
    type(lmultifab), intent(in), optional :: mask
    r = multifab_norm_l2_c(mf, 1, mf%nc, mask, all, local)
  end function multifab_norm_l2

  function multifab_norm_inf_doit(ap, lp) result(r)

    real(dp_t), pointer        :: ap(:,:,:,:)
    logical, pointer, optional :: lp(:,:,:,:)
    real(dp_t)                 :: r,r1
    integer                    :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    r1 = 0.0_dp_t

    if ( present(lp) ) then
       do n = lo(4), hi(4)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   if (lp(i,j,k,n)) r1 = max(r1,abs(ap(i,j,k,n)))
                end do
             end do
          end do
       end do
    else
       do n = lo(4), hi(4)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   r1 = max(r1,abs(ap(i,j,k,n)))
                end do
             end do
          end do
       end do
    end if

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
    type(mfiter)        :: mfi
    type(box)           :: bx

    lall   = .false.; if ( present(all)   ) lall   = all
    llocal = .false.; if ( present(local) ) llocal = local

    r1 = 0.0_dp_t

    if ( present(mask) ) then
       !$omp parallel private(lp,mp,i,n,mfi,bx) reduction(max:r1)
       call mfiter_build(mfi,mf,.true.)
       do while(next_tile(mfi,i))
          if ( lall ) then
             bx = get_growntilebox(mfi)
          else
             bx = get_tilebox(mfi)
          end if
          lp => dataptr(mask%fbs(i), bx)
          do n = comp, comp+nc-1
             mp => dataptr(mf%fbs(i), bx, n)
             r1 = max(r1, multifab_norm_inf_doit(mp,lp))
          end do
       end do
       !$omp end parallel
    else
       !$omp parallel private(mp,i,mfi) reduction(max:r1)
       call mfiter_build(mfi,mf,.true.)
       do while(next_tile(mfi,i))
          if ( lall ) then
             mp => dataptr(mf%fbs(i), get_growntilebox(mfi), comp, nc)
          else
             mp => dataptr(mf%fbs(i), get_tilebox(mfi), comp, nc)
          end if
          r1 = max(r1, multifab_norm_inf_doit(mp))
       end do
       !$omp end parallel
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
    type(mfiter) :: mfi
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0
    !$omp parallel private(mp,i,mfi) reduction(max:r1)
    call mfiter_build(mfi,mf,.true.)
    do while(next_tile(mfi,i))
       if ( lall ) then
          mp => dataptr(mf%fbs(i), get_growntilebox(mfi), comp, nc)
       else
          mp => dataptr(mf%fbs(i), get_tilebox(mfi), comp, nc)
       end if
       r1 = max(r1, maxval(abs(mp)))
    end do
    !$omp end parallel
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
    type(mfiter) :: mfi
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0
    !$omp parallel private(mp,i,mfi) reduction(+:r1)
    call mfiter_build(mfi,mf,.true.)
    do while(next_tile(mfi,i))
       if ( lall ) then
          mp => dataptr(mf%fbs(i), get_growntilebox(mfi), comp, nc)
       else
          mp => dataptr(mf%fbs(i), get_tilebox(mfi), comp, nc)
       end if
       r1 = r1 + sum(mp)
    end do
    !$omp end parallel
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
    type(mfiter) :: mfi
    lall = .false.; if ( present(all) ) lall = all
    r1 = 0
    !$omp parallel private(mp,i,mfi) reduction(+:r1)
    call mfiter_build(mfi,mf,.true.)
    do while(next_tile(mfi,i))
       if ( lall ) then
          mp => dataptr(mf%fbs(i), get_growntilebox(mfi))
       else
          mp => dataptr(mf%fbs(i), get_tilebox(mfi))
       end if
       r1 = r1 + count(mp)
    end do
    !$omp end parallel
    call parallel_reduce(r, r1, MPI_SUM)
  end function lmultifab_count

  subroutine multifab_div_div_c_doit(ap, bp)

    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)

    integer :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! ap = ap/bp

    do n = lo(4), hi(4)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) / bp(i,j,k,n)
             end do
          end do
       end do
    end do

  end subroutine multifab_div_div_c_doit

  subroutine multifab_div_div_s_doit(ap, b)

    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t)          :: b

    integer :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! ap = ap/b

    do n = lo(4), hi(4)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) * (1.0_dp_t / b)
             end do
          end do
       end do
    end do

  end subroutine multifab_div_div_s_doit

  subroutine multifab_div_div(a, b, ng)
    type(multifab), intent(inout) :: a
    type(multifab), intent(in   ) :: b
    integer, intent(in), optional :: ng
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i,lng
    type(mfiter) :: mfi
    type(box) :: bx
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= lng, b%ng >= lng,"not enough ghost cells in multifab_div_div")
    !$omp parallel private(ap,bp,i,mfi,bx)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       bx = get_growntilebox(mfi,lng)
       ap => dataptr(a%fbs(i), bx)
       bp => dataptr(b%fbs(i), bx)
       if ( any(bp == 0.0_dp_t) ) then
          call bl_error("MULTIFAB_DIV_DIV: divide by zero")
       end if
       call multifab_div_div_c_doit(ap, bp)
    end do
    !$omp end parallel
  end subroutine multifab_div_div

  subroutine multifab_div_div_s(a, b, ng)
    type(multifab), intent(inout) :: a
    real(dp_t), intent(in)  :: b
    integer, intent(in), optional :: ng
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i,lng
    type(mfiter) :: mfi
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= lng,"not enough ghost cells in multifab_div_div_s")
    if ( b == 0.0_dp_t ) then
       call bl_error("MULTIFAB_DIV_DIV_S: divide by zero")
    end if
    !$omp parallel private(ap,i,mfi)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       ap => dataptr(a%fbs(i), get_growntilebox(mfi,lng))
       call multifab_div_div_s_doit(ap, b)
    end do
    !$omp end parallel
  end subroutine multifab_div_div_s

  subroutine multifab_div_div_c(a, ia, b, ib, nc, ng)
    integer, intent(in) :: ia, ib
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i,lng
    type(mfiter) :: mfi
    type(box) :: bx
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= lng,"not enough ghost cells in multifab_div_div_c")
    !$omp parallel private(ap,bp,i,mfi,bx)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       bx = get_growntilebox(mfi, lng)
       ap => dataptr(a%fbs(i), bx, ia, nc)
       bp => dataptr(b%fbs(i), bx, ib, nc)
       if ( any(bp == 0.0_dp_t) ) then
          call bl_error("MULTIFAB_DIV_DIV: divide by zero")
       end if
       call multifab_div_div_c_doit(ap, bp)
    end do
    !$omp end parallel
  end subroutine multifab_div_div_c

  subroutine multifab_div_div_s_c(a, ia, b, nc, ng)
    integer, intent(in) :: ia
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    real(dp_t), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i,lng
    type(mfiter) :: mfi
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= lng,"not enough ghost cells in multifab_div_div_s_c")
    if ( b == 0.0_dp_t ) then
       call bl_error("MULTIFAB_DIV_DIV_S_C: divide by zero")
    end if
    !$omp parallel private(ap,i,mfi)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       ap => dataptr(a%fbs(i), get_growntilebox(mfi,lng), ia, nc)
       call multifab_div_div_s_doit(ap, b)
    end do
    !$omp end parallel
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
    type(mfiter) :: mfi
    type(box) :: bx
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= lng,"not enough ghost cells in multifab_div_s_c")
    if ( val == 0.0_dp_t ) then
       call bl_error("MULTIFAB_DIV_DIV_S_C: divide by zero")
    end if
    !$omp parallel private(ap,bp,i,mfi,bx)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       bx = get_growntilebox(mfi, lng)
       ap => dataptr(a%fbs(i), bx, ia, nc)
       bp => dataptr(b%fbs(i), bx, ib, nc)
       ap = bp*(1.0_dp_t/val)
    end do
    !$omp end parallel
  end subroutine multifab_div_s_c

  subroutine multifab_mult_mult_c_doit(ap, bp)
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer             :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! ap = ap*bp

    do n = lo(4), hi(4)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) * bp(i,j,k,n)
             end do
          end do
       end do
    end do

  end subroutine multifab_mult_mult_c_doit

  subroutine multifab_mult_mult_s_doit(ap, b)
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t)          :: b
    integer             :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! ap = ap*b

    do n = lo(4), hi(4)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) * b
             end do
          end do
       end do
    end do

  end subroutine multifab_mult_mult_s_doit

  subroutine multifab_mult_mult(a, b, ng)
    type(multifab), intent(inout) :: a
    type(multifab), intent(in   ) :: b
    integer, intent(in), optional :: ng
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i,lng
    type(mfiter) :: mfi
    type(box) :: bx
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= lng,"not enough ghost cells in multifab_mult_mult")
    !$omp parallel private(ap,bp,i,mfi,bx)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       bx = get_growntilebox(mfi, lng)
       ap => dataptr(a%fbs(i), bx)
       bp => dataptr(b%fbs(i), bx)
       call multifab_mult_mult_c_doit(ap, bp)
    end do
    !$omp end parallel
  end subroutine multifab_mult_mult
  subroutine multifab_mult_mult_s(a, b, ng)
    type(multifab), intent(inout) :: a
    real(dp_t)    , intent(in   ) :: b
    integer, intent(in), optional :: ng
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i,lng
    type(mfiter) :: mfi
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= lng,"not enough ghost cells in multifab_mult_mult_s")
    !$omp parallel private(ap,i,mfi)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       ap => dataptr(a%fbs(i), get_growntilebox(mfi, lng))
       call multifab_mult_mult_s_doit(ap, b)
    end do
    !$omp end parallel
  end subroutine multifab_mult_mult_s

  subroutine multifab_mult_mult_c(a, ia, b, ib, nc, ng)
    integer, intent(in) :: ia, ib
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i,lng
    type(mfiter) :: mfi
    type(box) :: bx
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= lng,"not enough ghost cells in multifab_mult_mult_c")
    !$omp parallel private(ap,bp,i,mfi,bx)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       bx = get_growntilebox(mfi, lng)
       ap => dataptr(a%fbs(i), bx, ia, nc)
       bp => dataptr(b%fbs(i), bx, ib, nc)
       call multifab_mult_mult_c_doit(ap, bp)
    end do
    !$omp end parallel
  end subroutine multifab_mult_mult_c

  subroutine multifab_mult_mult_s_c(a, ia, b, nc, ng)
    integer, intent(in) :: ia
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    real(dp_t), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i,lng
    type(mfiter) :: mfi
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= lng,"not enough ghost cells in multifab_mult_mult_s_c")
    !$omp parallel private(ap,i,mfi)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       ap => dataptr(a%fbs(i), get_growntilebox(mfi, lng), ia, nc)
       call multifab_mult_mult_s_doit(ap, b)
    end do
    !$omp end parallel
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
    type(mfiter) :: mfi
    type(box) :: bx
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= lng,"not enough ghost cells in multifab_mult_s_c")
    !$omp parallel private(ap,bp,i,mfi,bx)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       bx = get_growntilebox(mfi, lng)
       ap => dataptr(a%fbs(i), bx, ia, nc)
       bp => dataptr(b%fbs(i), bx, ib, nc)
       ap = bp * val
    end do
    !$omp end parallel
  end subroutine multifab_mult_s_c

  subroutine multifab_sub_sub_c_doit(ap, bp)

    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)

    integer :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! ap = ap - bp

    do n = lo(4), hi(4)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) - bp(i,j,k,n)
             end do
          end do
       end do
    end do

  end subroutine multifab_sub_sub_c_doit

  subroutine multifab_sub_sub_s_doit(ap, b)

    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t)          :: b

    integer :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! ap = ap - b

    do n = lo(4), hi(4)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) - b
             end do
          end do
       end do
    end do

  end subroutine multifab_sub_sub_s_doit

  subroutine multifab_sub_sub(a, b, ng)
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    integer, intent(in), optional :: ng
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i,lng
    type(mfiter) :: mfi
    type(box) :: bx
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= lng, b%ng >= lng, "not enough ghost cells in multifab_sub_sub")
    !$omp parallel private(ap,bp,i,mfi,bx)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       bx = get_growntilebox(mfi, lng)
       ap => dataptr(a%fbs(i), bx)
       bp => dataptr(b%fbs(i), bx)
       call multifab_sub_sub_c_doit(ap, bp)
    end do
    !$omp end parallel
  end subroutine multifab_sub_sub

  subroutine multifab_sub_sub_s(a, b, ng)
    type(multifab), intent(inout) :: a
    real(dp_t), intent(in)  :: b
    integer, intent(in), optional :: ng
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i,lng
    type(mfiter) :: mfi
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= lng, "not enough ghost cells in multifab_sub_sub_s")
    !$omp parallel private(ap,i,mfi)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       ap => dataptr(a%fbs(i), get_growntilebox(mfi,lng))
       call multifab_sub_sub_s_doit(ap, b)
    end do
    !$omp end parallel
  end subroutine multifab_sub_sub_s

  subroutine multifab_sub_sub_c(a, ia, b, ib, nc, ng)
    integer, intent(in) :: ia, ib
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i,lng
    type(mfiter) :: mfi
    type(box) :: bx
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= lng, b%ng >= lng, "not enough ghost cells in multifab_sub_sub_c")
    !$omp parallel private(i,mfi,ap,bp,bx)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       bx = get_growntilebox(mfi, lng)
       ap => dataptr(a%fbs(i), bx, ia, nc)
       bp => dataptr(b%fbs(i), bx, ib, nc)
       call multifab_sub_sub_c_doit(ap, bp)
    end do
    !$omp end parallel
  end subroutine multifab_sub_sub_c

  subroutine multifab_sub_sub_s_c(a, ia, b, nc, ng)
    integer, intent(in) :: ia
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    real(dp_t), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i,lng
    type(mfiter) :: mfi
    type(box) :: bx
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= lng, "not enough ghost cells in multifab_sub_sub_s_c")
    !$omp parallel private(i,mfi,ap,bx)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       bx = get_growntilebox(mfi,lng)
       ap => dataptr(a%fbs(i), bx, ia, nc)
       call multifab_sub_sub_s_doit(ap, b)
    end do
    !$omp end parallel
  end subroutine multifab_sub_sub_s_c

  subroutine multifab_plus_plus_c_doit(ap, bp)

    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)

    integer :: i, j, k, n, lo(4), hi(4)

    lo = lbound(ap)
    hi = ubound(ap)

    ! ap = ap + bp

    do n = lo(4), hi(4)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) + bp(i,j,k,n)
             end do
          end do
       end do
    end do

  end subroutine multifab_plus_plus_c_doit

  subroutine multifab_plus_plus_c(a, ia, b, ib, nc, ng)
    integer, intent(in) :: ia, ib
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    type(multifab), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    integer :: i,lng
    type(mfiter) :: mfi
    type(box) :: bx
    lng = 0; if ( present(ng) ) lng = ng
    if ( lng > 0 ) call bl_assert(a%ng >= lng, b%ng >= lng,"not enough ghost cells in multifab_plus_plus_c")
    !$omp parallel private(i,mfi,ap,bp,bx)
    call mfiter_build(mfi,a,.true.)
    do while (next_tile(mfi,i))
       bx = get_growntilebox(mfi,lng)
       ap => dataptr(a%fbs(i), bx, ia, nc)
       bp => dataptr(b%fbs(i), bx, ib, nc)
       call multifab_plus_plus_c_doit(ap, bp)
    end do
    !$omp end parallel
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
    do n = lo(4), hi(4)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ap(i,j,k,n) = ap(i,j,k,n) + b
             end do
          end do
       end do
    end do

  end subroutine multifab_plus_plus_s_doit

  subroutine multifab_plus_plus_s_c(a, ia, b, nc, ng)
    integer, intent(in) :: ia
    integer, intent(in)           :: nc
    integer, intent(in), optional :: ng
    type(multifab), intent(inout) :: a
    real(dp_t), intent(in)  :: b
    real(dp_t), pointer :: ap(:,:,:,:)
    integer :: i,lng
    type(mfiter) :: mfi
    lng = 0; if ( present(ng) ) lng = ng
    !$omp parallel private(i,mfi,ap)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,i))
       ap => dataptr(a%fbs(i), get_growntilebox(mfi,lng), ia, nc)
       call multifab_plus_plus_s_doit(ap, b)
    end do
    !$omp end parallel
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

  function multifab_min(mf, all, local) result(r)
    real(kind=dp_t) :: r
    type(multifab), intent(in) :: mf
    logical, intent(in), optional :: all, local
    integer :: i
    real(kind=dp_t) :: r1
    logical :: llocal, lall
    type(mfiter) :: mfi
    lall   = .false.; if (present(all)  ) lall  =  all
    llocal = .false.; if (present(local)) llocal=local
    r1 = +Huge(r1)
    !$omp parallel private(i,mfi) reduction(min:r1)
    call mfiter_build(mfi,mf,.true.)
    do while(next_tile(mfi,i))
       if (lall) then
          r1 = min(r1, min_val(mf%fbs(i), get_growntilebox(mfi)))
       else
          r1 = min(r1, min_val(mf%fbs(i), get_tilebox(mfi)))
       end if
    end do
    !$omp end parallel
    if (llocal) then
       r = r1
    else
       call parallel_reduce(r, r1, MPI_MIN)
    end if
  end function multifab_min
  function multifab_min_c(mf, c, nc, all, local) result(r)
    real(kind=dp_t) :: r
    type(multifab), intent(in) :: mf
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all, local
    real(kind=dp_t) :: r1
    integer :: i
    logical :: llocal, lall
    type(mfiter) :: mfi
    lall   = .false.; if (present(all)  ) lall  =  all
    llocal = .false.; if (present(local)) llocal=local
    r1 = +Huge(r1)
    !$omp parallel private(i,mfi) reduction(min:r1)
    call mfiter_build(mfi,mf,.true.)
    do while(next_tile(mfi,i))
       if (lall) then
          r1 = min(r1, min_val(mf%fbs(i), get_growntilebox(mfi), c, nc))
       else
          r1 = min(r1, min_val(mf%fbs(i), get_tilebox(mfi), c, nc))
       end if
    end do
    !$omp end parallel
    if (llocal) then
       r = r1
    else
       call parallel_reduce(r, r1, MPI_MIN)
    end if
  end function multifab_min_c

  function multifab_max(mf, all, local) result(r)
    real(kind=dp_t) :: r
    type(multifab), intent(in) :: mf
    logical, intent(in), optional :: all, local
    integer :: i
    real(kind=dp_t) :: r1
    logical :: llocal, lall
    type(mfiter) :: mfi
    lall   = .false.; if (present(all)  ) lall  =  all
    llocal = .false.; if (present(local)) llocal=local
    r1 = -Huge(r1)
    !$omp parallel private(i,mfi) reduction(max:r1)
    call mfiter_build(mfi,mf,.true.)
    do while(next_tile(mfi,i))
       if (lall) then
          r1 = max(r1, max_val(mf%fbs(i), get_growntilebox(mfi)))
       else
          r1 = max(r1, max_val(mf%fbs(i), get_tilebox(mfi)))
       end if
    end do
    !$omp end parallel
    if (llocal) then
       r = r1
    else
       call parallel_reduce(r, r1, MPI_MAX)
    end if
  end function multifab_max
  function multifab_max_c(mf, c, nc, all, local) result(r)
    real(kind=dp_t) :: r
    type(multifab), intent(in) :: mf
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all, local
    integer :: i
    real(kind=dp_t) :: r1
    logical :: llocal, lall
    type(mfiter) :: mfi
    lall   = .false.; if (present(all)  ) lall  =  all
    llocal = .false.; if (present(local)) llocal=local
    r1 = -Huge(r1)
    !$omp parallel private(i,mfi) reduction(max:r1)
    call mfiter_build(mfi,mf,.true.)
    do while(next_tile(mfi,i))
       if (lall) then
          r1 = max(r1, max_val(mf%fbs(i), get_growntilebox(mfi), c, nc))
       else
          r1 = max(r1, max_val(mf%fbs(i), get_tilebox(mfi), c, nc))
       end if
    end do
    !$omp end parallel
    if (llocal) then
       r = r1
    else
       call parallel_reduce(r, r1, MPI_MAX)
    end if
  end function multifab_max_c

  function multifab_equal(a, b) result (global_is_equal)
    type(multifab), intent(in) :: a
    type(multifab), intent(in) :: b
    real(dp_t), pointer           :: ap(:,:,:,:)
    real(dp_t), pointer           :: bp(:,:,:,:)
    
    integer :: ii, i, j, k, n, lo(4), hi(4)
    type(mfiter) :: mfi
    type(box) :: bx

    logical :: is_equal, global_is_equal

    is_equal = .true.

    if (.not. layout_equal(a%la, b%la)) then
       global_is_equal = .false.
       return
    endif

    !$omp parallel private(ap,bp,ii,i,j,k,n,lo,hi,mfi,bx) reduction (.and.: is_equal)
    call mfiter_build(mfi,a,.true.)
    do while(next_tile(mfi,ii))
       bx = get_tilebox(mfi)

       ap => dataptr(a%fbs(ii), bx)
       bp => dataptr(b%fbs(ii), bx)

       lo = lbound(ap)
       hi = ubound(ap)
       
       do n = lo(4), hi(4)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   if (.not. ap(i,j,k,n) == bp(i,j,k,n)) is_equal = .false.
                end do
             end do
          end do
       end do
    end do
    !$omp end parallel

    call parallel_reduce(global_is_equal, is_equal, MPI_LAND)

  end function multifab_equal


  subroutine multifab_iter_build(mfi, mf, tiling, tilesize)
    type(mfiter),   intent(inout) :: mfi
    type(multifab), intent(in)    :: mf
    logical, intent(in), optional :: tiling
    integer, intent(in), optional :: tilesize(:) 
    type(layout) :: la
    la = mf%la
    call iter_build_doit(mfi, la, mf%dim, mf%ng, mf%nodal, tiling, tilesize)
  end subroutine multifab_iter_build

  subroutine zmultifab_iter_build(mfi, mf, tiling, tilesize)
    type(mfiter),   intent(inout) :: mfi
    type(zmultifab), intent(in)   :: mf
    logical, intent(in), optional :: tiling
    integer, intent(in), optional :: tilesize(:)    
    type(layout) :: la
    la = mf%la
    call iter_build_doit(mfi, la, mf%dim, mf%ng, mf%nodal, tiling, tilesize)
  end subroutine zmultifab_iter_build

  subroutine imultifab_iter_build(mfi, mf, tiling, tilesize)
    type(mfiter),   intent(inout) :: mfi
    type(imultifab), intent(in)   :: mf
    logical, intent(in), optional :: tiling
    integer, intent(in), optional :: tilesize(:)    
    type(layout) :: la
    la = mf%la
    call iter_build_doit(mfi, la, mf%dim, mf%ng, mf%nodal, tiling, tilesize)
  end subroutine imultifab_iter_build

  subroutine lmultifab_iter_build(mfi, mf, tiling, tilesize)
    type(mfiter),   intent(inout) :: mfi
    type(lmultifab), intent(in)   :: mf
    logical, intent(in), optional :: tiling
    integer, intent(in), optional :: tilesize(:)    
    type(layout) :: la
    la = mf%la
    call iter_build_doit(mfi, la, mf%dim, mf%ng, mf%nodal, tiling, tilesize)
  end subroutine lmultifab_iter_build

  subroutine iter_build_doit(mfi, la, dim, ng, nodal, tiling, tilesize)
    use omp_module
    type(mfiter), intent(inout) :: mfi
    type(layout), intent(inout) :: la
    integer, intent(in) :: dim, ng
    logical, intent(in) :: nodal(:)
    logical, intent(in), optional :: tiling
    integer, intent(in), optional :: tilesize(:)

    integer :: ltilesize(3), tid, nthreads

    if (present(tilesize)) then
       ltilesize(1:size(tilesize)) = tilesize
    else if (present(tiling)) then
       if (tiling .eqv. .true.) then
          ltilesize = layout_get_tilesize()
       else
          ltilesize = (/ 1024000, 1024000, 1024000 /) ! large tile size turn off tiling
       end if
    else
       ltilesize = (/ 1024000, 1024000, 1024000 /) ! large tile size turn off tiling
    end if

    tid = omp_get_thread_num()
    nthreads = omp_get_num_threads()

    call init_layout_tilearray(mfi%ta, la, ltilesize, tid, nthreads)
    
    if (mfi%ta%dim .gt. 0) then
       mfi%dim = dim
       mfi%ng  = ng
       mfi%nodal(1:mfi%dim) = nodal(1:mfi%dim)

       mfi%it = 0
       mfi%ntiles = mfi%ta%ntiles
       mfi%lap => la%lap
    else
       mfi%it = 0
       mfi%ntiles = 0
    end if

    mfi%built = .true.  ! even when mfi%ntiles is 0.
  end subroutine iter_build_doit

  subroutine mfiter_reset(mfi)
    type(mfiter), intent(inout) :: mfi
    mfi%it = 0
  end subroutine mfiter_reset

  function next_tile(mfi,fi) result(r)
    logical :: r
    type(mfiter), intent(inout) :: mfi
    integer, intent(out) :: fi
    call bl_assert(mfi%built, "mfiter is not built")
    mfi%it = mfi%it + 1
    if (mfi%it .le. mfi%ntiles) then
       r = .true.
       fi = mfi%ta%lidx(mfi%it) ! current fab index
    else
       mfi%it = 0
       r = .false.
       fi = 0
    end if
  end function next_tile

  ! deprecated
  function more_tile(mfi) result(r)
    logical :: r
    type(mfiter), intent(inout) :: mfi
    call bl_assert(mfi%built, "mfiter is not built")
    mfi%it = mfi%it + 1
    if (mfi%it .le. mfi%ntiles) then
       r = .true.
    else
       mfi%it = 0
       r = .false.
    end if
  end function more_tile

  function get_fab_index(mfi) result(r)
    integer :: r
    type(mfiter), intent(in) :: mfi
    r = mfi%ta%lidx(mfi%it)
  end function get_fab_index

  function get_tilebox(mfi) result(r)
    type(box) :: r
    type(mfiter), intent(in) :: mfi
    integer :: i, gridhi(mfi%dim)
    r = make_box(mfi%ta%tilelo(:,mfi%it),mfi%ta%tilehi(:,mfi%it))
    if (any(mfi%nodal(1:mfi%dim))) then
       gridhi = upb(get_gridbox(mfi))
       do i = 1, mfi%dim
          if (mfi%nodal(i) .and. gridhi(i).eq.upb(r,i)) r%hi(i) = r%hi(i)+1
       end do
    end if
  end function get_tilebox

  function get_nodaltilebox(mfi, dir) result(r)
    type(box) :: r
    type(mfiter), intent(in) :: mfi
    integer, intent(in) :: dir
    integer :: gridhi, d, d0, d1
    r = get_tilebox(mfi)
    if (dir >=1 .and. dir <= mfi%dim) then
       d0 = dir
       d1 = dir
    else
       d0 = 1
       d1 = mfi%dim
    end if
    do d = d0, d1
       if (.not. mfi%nodal(d)) then
          gridhi = upb(get_gridbox(mfi), d)
          if (gridhi .eq. upb(r,d)) r%hi(d) = r%hi(d)+1
       end if
    end do
  end function get_nodaltilebox

  function get_allnodaltilebox(mfi) result(r)
    type(box) :: r
    type(mfiter), intent(in) :: mfi
    integer :: gridhi,dir
    r = get_tilebox(mfi)
    do dir=1,mfi%dim
       if (.not. mfi%nodal(dir)) then
          gridhi = upb(get_gridbox(mfi), dir)
          if (gridhi .eq. upb(r,dir)) r%hi(dir) = r%hi(dir)+1
       end if
    end do
  end function get_allnodaltilebox

  function get_growntilebox(mfi, ng_in) result(r)
    type(box) :: r
    type(mfiter), intent(in) :: mfi
    integer, intent(in), optional :: ng_in
    integer :: gridlo(mfi%dim), gridhi(mfi%dim), ng, i
    type(box) :: gridbox
    r = get_tilebox(mfi)
    ng = mfi%ng
    if (present(ng_in)) ng = ng_in
    r = get_tilebox(mfi)
    if (ng .gt. 0) then
       gridbox = get_gridbox(mfi)
       gridlo = lwb(gridbox)
       gridhi = upb(gridbox)
       do i=1,mfi%dim
          if (gridlo(i) .eq. lwb(r,i)) r%lo(i) = r%lo(i)-ng
          if (gridhi(i) .le. upb(r,i)) r%hi(i) = r%hi(i)+ng ! .le., not .eq., because of nodal
       end do
    end if
  end function get_growntilebox

  function get_grownnodaltilebox(mfi, dir, ng_in) result(r)
    type(box) :: r
    type(mfiter), intent(in) :: mfi
    integer, intent(in) :: dir
    integer, intent(in), optional :: ng_in
    integer :: gridhi, tilehi, d, d0, d1
    r = get_growntilebox(mfi, ng_in)
    if (dir >=1 .and. dir <= mfi%dim) then
       d0 = dir
       d1 = dir
    else
       d0 = 1
       d1 = mfi%dim
    end if
    do d = d0, d1
       if (.not. mfi%nodal(d)) then
          tilehi = upb(get_tilebox(mfi), d)
          gridhi = upb(get_gridbox(mfi), d)
          if (gridhi .eq. tilehi) r%hi(d) = r%hi(d)+1       
       end if
    end do
  end function get_grownnodaltilebox

  function get_gridbox(mfi) result (r)
    type(box) :: r
    type(mfiter), intent(in) :: mfi
    r = get_box(mfi%lap%bxa, mfi%ta%gidx(mfi%it))
  end function get_gridbox

end module multifab_module
