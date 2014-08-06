module layout_module

  use parallel
  use boxarray_module
  use bl_mem_stat_module

  implicit none

  integer, private, parameter :: LA_UNDF = 0
  integer, private, parameter :: LA_BASE = 1
  integer, private, parameter :: LA_CRSN = 2
  integer, private, parameter :: LA_PCHD = 3
  integer, private, parameter :: LA_PPRT = 4

  integer, parameter :: LA_KNAPSACK   = 101
  integer, parameter :: LA_ROUNDROBIN = 102
  integer, parameter :: LA_LOCAL      = 103
  integer, parameter :: LA_EXPLICIT   = 104

  integer, private :: verbose       = 0
  integer, private :: sfc_threshold = 0
  integer, private :: def_mapping   = LA_KNAPSACK

  type comm_dsc
     integer   :: nd = 0                 ! dst box number
     integer   :: ns = 0                 ! src box number
     type(box) :: dbx                    ! dst sub-box
     type(box) :: sbx                    ! src sub-box
     integer   :: pv = 0                 ! number of points in buf prior to this
     integer   :: sh(MAX_SPACEDIM+1) = 0 ! shape for data from rcvbuf
     integer   :: pr                     ! Processors number of src or dest
  end type comm_dsc

  type trns_dsc
     integer :: sz = 0              ! Size of chunk 
     integer :: pv = -1             ! Number points in buf prior to this
     integer :: pr = MPI_ANY_SOURCE ! src or destination processor
  end type trns_dsc

  type remote_conn
     integer                 :: svol = 0          ! numpts in snd volume
     integer                 :: rvol = 0          ! numpts in rcv volume
     integer                 :: nsnd = 0          ! Number of snd chunks
     integer                 :: nrcv = 0          ! Number of rcv chunks
     type(comm_dsc), pointer :: snd(:) => Null()
     type(comm_dsc), pointer :: rcv(:) => Null()
     integer                 :: nrp  = 0          ! Number of processes receiving from
     integer                 :: nsp  = 0          ! Number of processes sending to
     type(trns_dsc), pointer :: str(:) => Null()
     type(trns_dsc), pointer :: rtr(:) => Null()
  end type remote_conn

  type local_copy_desc
     integer   :: ns = 0    ! Source box in layout
     integer   :: nd = 0    ! Destination box in layout
     type(box) :: sbx       ! Sub-box for this copy
     type(box) :: dbx       ! Sub-box for this copy
  end type local_copy_desc

  type local_conn
     integer                        :: ncpy   ! Number of cpy chunks
     type(local_copy_desc), pointer :: cpy(:) => Null()
  end type local_conn

  type boxassoc
     integer                 :: dim      = 0       ! spatial dimension 1, 2, or 3
     integer                 :: nboxes   = 0       ! number of boxes
     integer                 :: grwth    = 0       ! growth factor
     integer                 :: idim     = 0       ! if 1, 2, or 3, grow only in that direction
     logical, pointer        :: nodal(:) => Null() ! nodal flag
     logical                 :: cross    = .false. ! cross/full stencil?
     type(local_conn)        :: l_con
     type(remote_conn)       :: r_con
     type(boxassoc), pointer :: next     => Null()
  end type boxassoc
  !
  ! Used to cache the boxarray used by multifab_fill_ghost_cells().
  !
  type fgassoc
     integer                :: dim   = 0       ! spatial dimension 1, 2, or 3
     integer                :: grwth = 0       ! growth factor
     type(boxarray)         :: ba
     type(fgassoc), pointer :: next  => Null()
  end type fgassoc

  type syncassoc
     integer                  :: dim      = 0       ! spatial dimension 1, 2, or 3
     integer                  :: nboxes   = 0       ! number of boxes
     integer                  :: grwth    = 0       ! growth factor
     logical                  :: lall     = .false. ! use valid region or everything
     logical, pointer         :: nodal(:) => Null() ! nodal flag
     type(local_conn)         :: l_con
     type(remote_conn)        :: r_con
     type(syncassoc), pointer :: next     => Null()
  end type syncassoc

  type copyassoc
     integer                   :: dim        = 0       ! spatial dimension 1, 2, or 3
     integer                   :: reused     = 0
     integer                   :: hash
     logical, pointer          :: nd_dst(:)  => Null() ! dst nodal flag
     logical, pointer          :: nd_src(:)  => Null() ! src nodal flag
     type(local_conn)          :: l_con
     type(remote_conn)         :: r_con
     type(copyassoc),  pointer :: next       => Null()
     type(boxarray)            :: ba_src
     type(boxarray)            :: ba_dst
     integer, pointer          :: prc_src(:) => Null()
     integer, pointer          :: prc_dst(:) => Null()
  end type copyassoc
  !
  ! Note that a fluxassoc contains two copyassocs.
  !
  type fluxassoc
     integer                   :: dim = 0             ! spatial dimension 1, 2, or 3
     integer                   :: side
     integer                   :: ir(3)
     type(box)                 :: crse_domain
     logical,          pointer :: nd_dst(:) => Null() ! dst nodal flag
     logical,          pointer :: nd_src(:) => Null() ! src nodal flag
     type(copyassoc)           :: flux
     type(copyassoc)           :: mask
     type(box),        pointer :: fbxs(:) => Null()
     type(fluxassoc),  pointer :: next    => Null()
     type(layout_rep), pointer :: lap_dst => Null()
     type(layout_rep), pointer :: lap_src => Null()
  end type fluxassoc
  !
  ! Used by layout_get_box_intersector().
  !
  type box_intersector
     integer   :: i
     type(box) :: bx
  end type box_intersector

  type box_hash_bin
     integer, pointer :: iv(:) => Null()
  end type box_hash_bin
  !
  ! Global list of copyassoc's used by multifab copy routines.
  !
  type(copyassoc), pointer, save, private :: the_copyassoc_head => Null()

  integer, save, private :: the_copyassoc_cnt = 0  ! Count of copyassocs on list.
  integer, save, private :: the_copyassoc_max = 25 ! Maximum # copyassocs allowed on list.
  !
  ! Global list of fluxassoc's used by ml_crse_contrib()
  !
  type(fluxassoc), pointer, save, private :: the_fluxassoc_head => Null()

  type layout
     integer                   :: la_type =  LA_UNDF
     type(layout_rep), pointer :: lap     => Null()
  end type layout

  !! Defines the box distribution and box connectivity of a boxarray.
  type layout_rep
     integer                         :: dim    = 0            ! Spatial dimension: 1, 2, or 3.
     integer                         :: id     = 0
     integer                         :: nboxes = 0            ! Total count of boxes.
     integer                         :: nlocal = 0            ! Number of boxes we own.
     type(box)                       :: pd                    ! Problem Domain.
     logical, pointer                :: pmask(:)    => Null() ! Periodic mask.
     integer, pointer                :: prc(:)      => Null() ! Which processor owns each box.
     integer, pointer                :: idx(:)      => Null() ! Global indices of boxes we own.
     type(boxarray)                  :: bxa
     type(boxassoc), pointer         :: bxasc       => Null()
     type(fgassoc), pointer          :: fgasc       => Null()
     type(syncassoc), pointer        :: snasc       => Null()
     type(coarsened_layout), pointer :: crse_la     => Null()
     type(pn_layout), pointer        :: pn_children => Null()
     ! Box Hashing
     integer                         :: crsn                = -1
     integer                         :: plo(MAX_SPACEDIM)   = 0
     integer                         :: phi(MAX_SPACEDIM)   = 0
     integer                         :: vshft(MAX_SPACEDIM) = 0
     type(box_hash_bin), pointer     :: bins(:,:,:)         => Null()
  end type layout_rep

  !! A layout that is derived by coarsening an existing layout,
  !! The processor distribution and the number of boxes will be
  !! the same as for the parent layout.  The intent is to be used
  !! in multigrid solvers that keep coarsened grids on the same
  !! processor as their parent in the hierarchy.
  type coarsened_layout
     integer                         :: dim = 0
     integer, pointer                :: crse(:) => Null()
     type(layout)                    :: la
     type(coarsened_layout), pointer :: next => Null()
  end type coarsened_layout

  type pn_layout
     integer                  :: dim = 0
     integer, pointer         :: refr(:) => Null()
     type(layout)             :: la
     type(pn_layout), pointer :: next => Null()
  end type pn_layout

  integer, private :: g_layout_next_id = 0

  interface built_q
     module procedure layout_built_q
     module procedure boxassoc_built_q
     module procedure fgassoc_built_q
     module procedure syncassoc_built_q
     module procedure copyassoc_built_q
     module procedure fluxassoc_built_q
  end interface

  interface build
     module procedure layout_build_ba
  end interface

  interface destroy
     module procedure layout_destroy
  end interface

  interface print
     module procedure layout_print
  end interface

  interface local
     module procedure layout_local
  end interface

  interface remote
     module procedure layout_remote
  end interface

  interface nboxes
     module procedure layout_nboxes
  end interface

  interface nlocal
     module procedure layout_nlocal
  end interface

  interface local_index
     module procedure layout_local_index
  end interface

  interface global_index
     module procedure layout_global_index
  end interface

  interface get_box
     module procedure layout_get_box
  end interface

  interface get_boxarray
     module procedure layout_boxarray
  end interface

  interface not_equal
     module procedure layout_not_equal
  end interface
  interface operator(.NE.)
     module procedure layout_not_equal
  end interface

  interface equal
     module procedure layout_equal
  end interface
  interface operator(.EQ.)
     module procedure layout_equal
  end interface

  interface get_proc
     module procedure layout_get_proc
     module procedure layout_get_proc_v
  end interface

  interface get_pd
     module procedure layout_get_pd
  end interface

  interface get_dim
     module procedure layout_dim
  end interface

  interface get_pmask
     module procedure layout_get_pmask
  end interface

  private :: greater_i, layout_next_id, layout_rep_build, layout_rep_destroy

  type(mem_stats), private, save :: la_ms
  type(mem_stats), private, save :: bxa_ms
  type(mem_stats), private, save :: avx_ms
  type(mem_stats), private, save :: fgx_ms
  type(mem_stats), private, save :: snx_ms
  type(mem_stats), private, save :: cpx_ms
  type(mem_stats), private, save :: flx_ms

contains

  subroutine layout_set_sfc_threshold(v)
    integer, intent(in) :: v
    sfc_threshold = v
  end subroutine layout_set_sfc_threshold
  pure function layout_get_sfc_threshold() result(r)
    integer :: r
    r = sfc_threshold
  end function layout_get_sfc_threshold

  subroutine layout_set_copyassoc_max(v)
    use bl_error_module
    integer, intent(in) :: v
    call bl_assert(the_copyassoc_max .gt. 1, "the_copyassoc_max must be > 1")
    the_copyassoc_max = v
  end subroutine layout_set_copyassoc_max
  pure function layout_get_copyassoc_max() result(r)
    integer :: r
    r = the_copyassoc_max
  end function layout_get_copyassoc_max

  subroutine layout_set_verbosity(v)
    integer, intent(in) :: v
    verbose = v
  end subroutine layout_set_verbosity
  pure function layout_get_verbosity() result(r)
    integer :: r
    r = verbose
  end function layout_get_verbosity

  subroutine layout_set_mapping(mapping)
    integer, intent(in) :: mapping
    def_mapping = mapping
  end subroutine layout_set_mapping
  pure function layout_get_mapping() result(r)
    integer :: r
    r = def_mapping
  end function layout_get_mapping

  subroutine layout_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    la_ms = ms
  end subroutine layout_set_mem_stats
  subroutine boxassoc_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    bxa_ms = ms
  end subroutine boxassoc_set_mem_stats
  subroutine fgassoc_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    fgx_ms = ms
  end subroutine fgassoc_set_mem_stats
  subroutine syncassoc_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    snx_ms = ms
  end subroutine syncassoc_set_mem_stats
  subroutine copyassoc_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    cpx_ms = ms
  end subroutine copyassoc_set_mem_stats
  subroutine fluxassoc_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    flx_ms = ms
  end subroutine fluxassoc_set_mem_stats

  function layout_mem_stats() result(r)
    type(mem_stats) :: r
    r = la_ms
  end function layout_mem_stats
  function boxassoc_mem_stats() result(r)
    type(mem_stats) :: r
    r = bxa_ms
  end function boxassoc_mem_stats
  function fgassoc_mem_stats() result(r)
    type(mem_stats) :: r
    r = fgx_ms
  end function fgassoc_mem_stats
  function syncassoc_mem_stats() result(r)
    type(mem_stats) :: r
    r = snx_ms
  end function syncassoc_mem_stats
  function copyassoc_mem_stats() result(r)
    type(mem_stats) :: r
    r = cpx_ms
  end function copyassoc_mem_stats
  function fluxassoc_mem_stats() result(r)
    type(mem_stats) :: r
    r = flx_ms
  end function fluxassoc_mem_stats

  function layout_next_id() result(r)
    integer :: r
    g_layout_next_id = g_layout_next_id + 1
    r = g_layout_next_id
  end function layout_next_id

  pure function layout_not_equal(a,b) result(r)
    type(layout), intent(in) :: a, b
    logical :: r
    r = .not. associated(a%lap, b%lap)
  end function layout_not_equal

  pure function layout_equal(a,b) result(r)
    type(layout), intent(in) :: a, b
    logical :: r
    r = associated(a%lap, b%lap)
  end function layout_equal

  pure function layout_built_q(la) result(r)
    logical :: r
    type(layout), intent(in) :: la
    r = associated(la%lap)
  end function layout_built_q

  pure function layout_dim(la) result(r)
    integer :: r
    type(layout), intent(in) :: la
    r = la%lap%dim
  end function layout_dim

  pure function layout_nboxes(la) result(r)
    integer :: r
    type(layout), intent(in) :: la
    r = la%lap%nboxes
  end function layout_nboxes

  pure function layout_nlocal(la) result(r)
    integer :: r
    type(layout), intent(in) :: la
    r = la%lap%nlocal
  end function layout_nlocal
  !
  ! Given a global index "i" into the boxarray on which
  ! the layout is built, returns the corresponding local
  ! index "r" at which that index is stored in la&lap%idx(:).
  ! It's an error if the layout does not "own" index "r"
  ! or if the index "i" is not in the range of boxes for the
  ! layout. nboxes is the number of boxes in the associated layout.
  !
  function layout_local_index(la,i) result(r)

    use bl_error_module

    integer,      intent(in) :: i
    type(layout), intent(in) :: la
    integer                  :: r

    call bl_assert(i >= 1 .and. i <= la%lap%nboxes, "layout_local_index: invalid global index")

    if (parallel_nprocs() == 1) then
       call bl_assert(la%lap%idx(i) == i, "layout_local_index: how did this happen?")
       r = i
    else
       r = bsearch(la%lap%idx,i)
       call bl_assert(r >= 1, "layout_local_index: no corresponding local index")
    endif

  contains
      !
      ! Returns -1 if "val" is not found in array "arr".
      !
      ! "arr" is assumed to be sorted from smallest to largest.
      !
      pure function bsearch (arr,val) result (r)
        integer, intent(in) :: arr(:), val
        integer             :: r, lo, hi, mid
        r  = -1
        lo = lbound(arr,1)
        hi = ubound(arr,1)
        do while (lo <= hi)
           mid = (lo + hi) / 2
           if (arr(mid) == val) then
              r = mid
              exit
           else if (arr(mid) > val) then
              hi = mid - 1
           else
              lo = mid + 1
           end if
        end do
      end function bsearch

  end function layout_local_index

  pure function layout_global_index(la,i) result(r)
    integer,      intent(in) :: i
    type(layout), intent(in) :: la
    integer                  :: r
    r = la%lap%idx(i)
  end function layout_global_index

  pure function layout_get_pd(la) result(r)
    type(box) :: r
    type(layout), intent(in) :: la
    r = la%lap%pd
  end function layout_get_pd

  function layout_boxarray(la) result(r)
    type(layout), intent(in) :: la
    type(boxarray) :: r
    r = la%lap%bxa
  end function layout_boxarray

  pure function layout_get_pmask(la) result(r)
    type(layout), intent(in) :: la
    logical :: r(la%lap%dim)
    r = la%lap%pmask
  end function layout_get_pmask

  subroutine layout_rep_build(lap, ba, pd, pmask, mapping, explicit_mapping)
    use bl_error_module
    type(layout_rep), intent(out) :: lap
    type(boxarray), intent(in) :: ba
    type(box), intent(in) :: pd
    logical, intent(in) :: pmask(:)
    integer, intent(in), optional :: mapping
    integer, intent(in), optional :: explicit_mapping(:)
    integer :: lmapping, i, j

    lmapping = def_mapping; if ( present(mapping) ) lmapping = mapping
    if ( present(explicit_mapping) ) then
       if ( present(mapping) ) then
          if ( mapping /= LA_EXPLICIT ) then
             call bl_error("layout_rep_build():explicit_mapping doesn't match mapping")
          end if
       end if
       lmapping = LA_EXPLICIT
    end if

    call boxarray_build_copy(lap%bxa, ba)

    lap%dim    = get_dim(lap%bxa)
    lap%nboxes = nboxes(lap%bxa)
    lap%id     = layout_next_id()
    lap%pd     = pd
    allocate(lap%pmask(lap%dim))
    lap%pmask = pmask

    allocate(lap%prc(lap%nboxes))

    select case (lmapping)
    case (LA_EXPLICIT)
       if ( .not. present(explicit_mapping) ) then
          call bl_error("layout_rep_build(): mapping explicit but no explicit_mapping")
       end if
       if ( size(lap%prc) /= size(explicit_mapping) ) then
          call bl_error("layout_rep_build(): incommensurate explicit mapping")
       end if
       lap%prc = explicit_mapping
    case (LA_LOCAL)
       lap%prc = parallel_myproc()
    case (LA_ROUNDROBIN)
       call layout_roundrobin(lap%prc, ba)
    case (LA_KNAPSACK)
       call layout_knapsack(lap%prc, ba)
    case default
       call bl_error("layout_rep_build(): unknown mapping:", lmapping)
    end select

    lap%nlocal = 0
    do i = 1, lap%nboxes
       if (parallel_myproc() == lap%prc(i)) lap%nlocal = lap%nlocal + 1
    end do

    allocate(lap%idx(lap%nlocal))

    j = 1
    do i = 1, lap%nboxes
       if (parallel_myproc() == lap%prc(i)) then
          lap%idx(j) = i
          j = j + 1
       end if
    end do
  end subroutine layout_rep_build

  subroutine layout_flush_copyassoc_cache ()
    type(copyassoc), pointer :: cp, ncp
    integer                  :: i

    i  =  1
    cp => the_copyassoc_head
    do while ( associated(cp) )
       ncp => cp%next
       call copyassoc_destroy(cp)
       deallocate(cp)
       cp => ncp
       i  =  i + 1
    end do

    the_copyassoc_cnt  =  0
    the_copyassoc_head => Null()
  end subroutine layout_flush_copyassoc_cache

  recursive subroutine layout_rep_destroy(lap, la_type)
    type(layout_rep), pointer :: lap
    integer, intent(in) :: la_type
    type(coarsened_layout), pointer :: clp, oclp
    type(pn_layout), pointer :: pnp, opnp
    type(boxassoc),  pointer :: bxa, obxa
    type(fgassoc),   pointer :: fgxa, ofgxa
    type(syncassoc), pointer :: snxa, osnxa
    type(fluxassoc), pointer :: fla, nfla, pfla
    if ( la_type /= LA_CRSN ) then
       deallocate(lap%prc)
       deallocate(lap%idx)
    end if
    call destroy(lap%bxa)
    if ( (la_type == LA_BASE) .or. (la_type == LA_PCHD) ) then
       deallocate(lap%pmask)
    end if
    clp => lap%crse_la
    do while ( associated(clp) )
       oclp => clp%next
       deallocate(clp%crse)
       call layout_rep_destroy(clp%la%lap, LA_CRSN)
       deallocate(clp)
       clp => oclp
    end do
    pnp => lap%pn_children
    do while ( associated(pnp) )
       opnp => pnp%next
       deallocate(pnp%refr)
       call layout_rep_destroy(pnp%la%lap, LA_PCHD)
       deallocate(pnp)
       pnp  => opnp
    end do
    !
    ! Get rid of boxassocs.
    !
    bxa => lap%bxasc
    do while ( associated(bxa) )
       obxa => bxa%next
       call boxassoc_destroy(bxa)
       deallocate(bxa)
       bxa => obxa
    end do
    !
    ! Get rid of fgassocs.
    !
    fgxa => lap%fgasc
    do while ( associated(fgxa) )
       ofgxa => fgxa%next
       call fgassoc_destroy(fgxa)
       deallocate(fgxa)
       fgxa => ofgxa
    end do
    !
    ! Get rid of syncassocs.
    !
    snxa => lap%snasc
    do while ( associated(snxa) )
       osnxa => snxa%next
       call syncassoc_destroy(snxa)
       deallocate(snxa)
       snxa => osnxa
    end do
    !
    ! Remove any boxarray hash.
    !
    call clear_box_hash_bin(lap)
    !
    ! Remove all fluxassoc's associated with this layout_rep.
    !
    fla  => the_fluxassoc_head
    pfla => Null()
    do while ( associated(fla) )
       nfla => fla%next
       if ( associated(lap, fla%lap_src) .or. associated(lap, fla%lap_dst) ) then
          if ( associated(fla, the_fluxassoc_head) ) then
             the_fluxassoc_head => fla%next
          else
             pfla%next => nfla
          end if
          call fluxassoc_destroy(fla)
          deallocate(fla)
       else
          if ( .not. associated(pfla) ) then
             pfla => the_fluxassoc_head
          else
             pfla => pfla%next
          end if
       end if
       fla => nfla
    end do

    deallocate(lap)
  end subroutine layout_rep_destroy

  subroutine layout_build_ba(la, ba, pd, pmask, mapping, explicit_mapping)
    type(layout)  , intent(  out) :: la
    type(boxarray), intent(in   ) :: ba
    type(box)     , intent(in   ) :: pd
    logical, intent(in), optional :: pmask(:)
    integer, intent(in), optional :: mapping
    integer, intent(in), optional :: explicit_mapping(:)

    logical :: lpmask(get_dim(ba))
    lpmask = .false.; if ( present(pmask) ) lpmask = pmask

    allocate(la%lap)
    la%la_type = LA_BASE
    call layout_rep_build(la%lap, ba, pd, lpmask, mapping, explicit_mapping)
    call mem_stats_alloc(la_ms)
  end subroutine layout_build_ba

  subroutine layout_destroy(la)
    use bl_error_module
    type(layout), intent(inout) :: la
    if ( la%la_type /= LA_BASE ) call bl_error("layout_destroy(): confused")
    call layout_rep_destroy(la%lap, LA_BASE)
    call mem_stats_dealloc(la_ms)
  end subroutine layout_destroy

  subroutine layout_build_pn(lapn, la, ba, rr, mapping, explicit_mapping)
    use bl_error_module
    type(layout), intent(out)   :: lapn
    type(layout), intent(inout) :: la
    type(boxarray), intent(in) :: ba
    integer, intent(in) :: rr(:)
    integer, intent(in), optional :: mapping
    integer, intent(in), optional :: explicit_mapping(:)
    type(pn_layout), pointer :: pla
    type(box) :: rpd

    if ( size(rr) /= la%lap%dim ) then
       call bl_error("layout_build_pn(): incommensurate refinement ratio")
    end if

    ! This is wrong: I need to make sure that the boxarray and the
    ! refinement match.  This will be OK until we do regridding

    pla => la%lap%pn_children
    do while ( associated(pla) )
       if ( all(pla%refr == rr) ) then
          lapn = pla%la
          return
       end if
       pla => pla%next
    end do

    ! Should also check for the proper nestedness

    allocate(pla)
    allocate(pla%refr(la%lap%dim))
    pla%dim = la%lap%dim
    pla%refr = rr

    pla%la%la_type = LA_PCHD

    allocate(pla%la%lap)
    rpd = refine(la%lap%pd, pla%refr)
    call layout_rep_build(pla%la%lap, ba, rpd, la%lap%pmask, &
        mapping, explicit_mapping)

    ! install the new coarsened layout into the layout
    pla%next => la%lap%pn_children
    la%lap%pn_children => pla
    lapn = pla%la
  end subroutine layout_build_pn

  subroutine layout_build_coarse(lac, la, cr)
    use bl_error_module
    type(layout), intent(out)   :: lac
    type(layout), intent(inout) :: la
    integer, intent(in) :: cr(:)
    type(coarsened_layout), pointer :: cla

    if ( size(cr) /= la%lap%dim ) then
       call bl_error("layout_build_coarse(): incommensurate cr")
    end if
    ! check if la already has this coarsened_layout
    cla => la%lap%crse_la
    do while ( associated(cla) )
       if ( all(cla%crse == cr) ) then
          lac = cla%la
          return
       end if
       cla => cla%next
    end do

    ! not, found, must build another
    allocate(cla)
    allocate(cla%crse(la%lap%dim))
    cla%dim = la%lap%dim
    cla%crse = cr

    ! allocate(cla%la)
    cla%la%la_type = LA_CRSN
    allocate(cla%la%lap)
    cla%la%lap%dim = la%lap%dim
    cla%la%lap%id  = layout_next_id()
    cla%la%lap%nboxes = la%lap%nboxes
    cla%la%lap%nlocal = la%lap%nlocal
    cla%la%lap%pd = coarsen(la%lap%pd, cla%crse)
    cla%la%lap%pmask => la%lap%pmask

    cla%la%lap%prc => la%lap%prc
    cla%la%lap%idx => la%lap%idx

    call boxarray_build_copy(cla%la%lap%bxa, la%lap%bxa)

    call boxarray_coarsen(cla%la%lap%bxa, cla%crse)

    ! install the new coarsened layout into the layout
    cla%next => la%lap%crse_la
    la%lap%crse_la => cla
    lac = cla%la
  end subroutine layout_build_coarse

  pure function layout_remote(la, i) result(r)
    type(layout), intent(in) :: la
    integer, intent(in) :: i
    logical :: r
    r = la%lap%prc(i) /= parallel_myproc()
  end function layout_remote

  pure function layout_local(la, i) result(r)
    type(layout), intent(in) :: la
    integer, intent(in) :: i
    logical :: r
    r = la%lap%prc(i) == parallel_myproc()
  end function layout_local

  pure function layout_get_box(la, i) result(r)
    type(layout), intent(in) :: la
    integer, intent(in) :: i
    type(box) :: r
    r = get_box(la%lap%bxa, i)
  end function layout_get_box

  pure function layout_get_proc(la, i) result(r)
    type(layout), intent(in) :: la
    integer, intent(in) :: i
    integer :: r
    r = la%lap%prc(i)
  end function layout_get_proc

  pure function layout_get_proc_v(la) result(r)
    type(layout), intent(in) :: la
    integer :: r(size(la%lap%prc))
    r = la%lap%prc
  end function layout_get_proc_v

  function greater_i(a,b) result(r)
    logical :: r
    integer, intent(in) :: a, b
    r = a > b
  end function greater_i

  subroutine layout_roundrobin(prc, bxs)
    use fab_module
    use sort_i_module
    use bl_error_module
    integer,        intent(out), dimension(:) :: prc
    type(boxarray), intent(in )               :: bxs

    integer          :: i
    integer, pointer :: luc(:)
    integer, allocatable :: ibxs(:), idx(:)

    if ( nboxes(bxs) /= size(prc) ) call bl_error('layout_roundrobin: how did this happen?')

    allocate(ibxs(nboxes(bxs)), idx(nboxes(bxs)))

    do i = 1, size(ibxs,1)
       ibxs(i) = volume(get_box(bxs,i))
    end do

    call sort(ibxs, idx, greater_i)

    luc => least_used_cpus()

    do i = 0, size(prc,1)-1
       prc(i) = luc(mod(idx(i)-1,parallel_nprocs()))
    end do

    deallocate(luc)

  end subroutine layout_roundrobin

  subroutine layout_knapsack(prc, bxs)
    use fab_module
    use bl_error_module
    use knapsack_module
    integer,        intent(out), dimension(:) :: prc
    type(boxarray), intent(in )               :: bxs

    integer              :: i
    integer, pointer     :: luc(:)
    integer, allocatable :: ibxs(:), tprc(:)
    type(box), pointer   :: pbxs(:)

    if ( nboxes(bxs) /= size(prc) ) call bl_error('layout_knapsack: how did this happen?')

    allocate(ibxs(nboxes(bxs)), tprc(nboxes(bxs)))

    do i = 1, size(ibxs,1)
       ibxs(i) = volume(get_box(bxs,i))
    end do

    if ( (nboxes(bxs)/parallel_nprocs()) > sfc_threshold ) then
       !
       ! Use Morton space-filling-curve distribution if we have "enough" grids.
       !
       pbxs => dataptr(bxs)
       call sfc_i(tprc, ibxs, pbxs, parallel_nprocs())
    else
       !
       ! knapsack_i() sorts boxes so that CPU 0 contains largest volume & CPU nprocs-1 the least.
       !
       call knapsack_i(tprc, ibxs, parallel_nprocs())
    end if

    luc => least_used_cpus()

    do i = 1, size(prc,1)
       prc(i) = luc(tprc(i))
    end do

    deallocate(luc)

  end subroutine layout_knapsack

  function layout_efficiency(la, np) result(r)
    type(layout), intent(in) :: la
    real(kind=dp_t) :: r
    integer, intent(in), optional :: np
    real(kind=dp_t) :: weights(la%lap%nboxes)
    real(kind=dp_t) :: p_max_weight
    integer :: i, lnp
    lnp = parallel_nprocs(); if ( present(np) ) lnp = np
    do i = 1, nboxes(la%lap%bxa)
       weights(i) = dvolume(get_box(la%lap%bxa,i))
    end do
    p_max_weight = -Huge(p_max_weight)
    do i = 0, lnp-1
       p_max_weight = max(p_max_weight, sum(weights, mask = la%lap%prc==i))
    end do
    r = sum(weights)/lnp/p_max_weight
  end function layout_efficiency

  subroutine layout_print(la, str, unit, skip)
    use bl_IO_module
    type(layout), intent(in) :: la
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: skip
    integer :: un, i
    un = unit_stdout(unit)
    call unit_skip(unit, skip)
    write(unit=un, fmt='("LAYOUT[(* ")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(": ", A)') str
    else
       write(unit=un, fmt='()')
    end if
    if ( .not. associated(la%lap) ) then
       call unit_skip(unit, skip)
       write(unit=un, fmt='(" empty *)]")')
    else
       call unit_skip(unit, skip)
       write(unit=un, fmt='(" ID = ",i0)', advance = 'no') la%lap%id
       call unit_skip(unit, skip)
       write(unit=un, fmt='(" DIM     = ",i2)') la%lap%dim
       call unit_skip(unit, skip)
       write(unit=un, fmt='(" NBOXES  = ",i3)') la%lap%nboxes
       call unit_skip(unit, skip)
       write(unit=un, fmt='(" PD      = ",i2)', advance = 'no')
       call print(la%lap%pd, unit = unit)
       do i = 1, nboxes(la)
          call unit_skip(unit=un, skip = unit_get_skip(skip) + 1)
          write(unit=un, fmt = '(I0,": ")', advance = 'no') i
          call print(get_box(la,i), unit = un, advance = 'no')
          write(unit=un, fmt = '(" ",I0)') get_proc(la,i)
       end do
       write(unit=un, fmt = '(" *)]")')
    end if
  end subroutine layout_print

  pure function boxassoc_check(bxa, ng, nodal, cross, idim) result(r)
    type(boxassoc), intent(in) :: bxa
    integer,        intent(in) :: ng
    logical,        intent(in) :: nodal(:)
    logical,        intent(in) :: cross
    integer, intent(in), optional :: idim
    logical                    :: r
    r = (bxa%grwth == ng) .and. all(bxa%nodal .eqv. nodal) .and. (bxa%cross .eqv. cross)
    if (present(idim)) then
       r = r .and. (bxa%idim == idim)
    end if
  end function boxassoc_check

  pure function fgassoc_check(fgxa, ng) result(r)
    type(fgassoc), intent(in) :: fgxa
    integer,       intent(in) :: ng
    logical                   :: r
    r = (fgxa%grwth == ng)
  end function fgassoc_check

  pure function syncassoc_check(snxa, ng, nodal, lall) result(r)
    type(syncassoc), intent(in) :: snxa
    integer,         intent(in) :: ng
    logical,         intent(in) :: nodal(:)
    logical,         intent(in) :: lall
    logical                     :: r
    r = (snxa%grwth == ng) .and. all(snxa%nodal .eqv. nodal) .and. (snxa%lall .eqv. lall)
  end function syncassoc_check

  function layout_boxassoc(la, ng, nodal, cross, idim) result(r)
    type(boxassoc)               :: r
    type(layout) , intent(inout) :: la
    integer, intent(in)          :: ng
    logical, intent(in)          :: nodal(:)
    logical, intent(in)          :: cross
    integer, intent(in),optional :: idim

    type(boxassoc), pointer :: bp

    bp => la%lap%bxasc
    do while ( associated(bp) )
       if ( boxassoc_check(bp, ng, nodal, cross, idim) ) then
          r = bp
          return
       end if
       bp => bp%next
    end do
    !
    ! Have to build one.
    !
    allocate (bp)
    call boxassoc_build(bp, la%lap, ng, nodal, cross, idim=idim)
    bp%next      => la%lap%bxasc
    la%lap%bxasc => bp
    r = bp
  end function layout_boxassoc

  function layout_fgassoc(la, ng) result(r)
    type(fgassoc)                :: r
    type(layout) , intent(inout) :: la
    integer, intent(in)          :: ng

    type(fgassoc), pointer :: fgp

    fgp => la%lap%fgasc
    do while ( associated(fgp) )
       if ( fgassoc_check(fgp, ng) ) then
          r = fgp
          return
       end if
       fgp => fgp%next
    end do
    !
    ! Have to go build one.
    !
    allocate (fgp)
    call fgassoc_build(fgp, la, ng)
    fgp%next     => la%lap%fgasc
    la%lap%fgasc => fgp
    r = fgp
  end function layout_fgassoc

  function layout_syncassoc(la, ng, nodal, lall) result(r)
    type(syncassoc)              :: r
    type(layout) , intent(inout) :: la
    integer, intent(in)          :: ng
    logical, intent(in)          :: nodal(:)
    logical, intent(in)          :: lall
    type(syncassoc), pointer     :: sp

    sp => la%lap%snasc
    do while ( associated(sp) )
       if ( syncassoc_check(sp, ng, nodal, lall) ) then
          r = sp
          return
       end if
       sp => sp%next
    end do
    !
    ! Have to go looking for it.
    !
    allocate (sp)
    call syncassoc_build(sp, la%lap, ng, nodal, lall)
    sp%next      => la%lap%snasc
    la%lap%snasc => sp
    r = sp
  end function layout_syncassoc

  pure function boxassoc_built_q(bxasc) result(r)
    logical :: r
    type(boxassoc), intent(in) :: bxasc
    r = bxasc%dim /= 0
  end function boxassoc_built_q

  pure function fgassoc_built_q(fgasc) result(r)
    logical :: r
    type(fgassoc), intent(in) :: fgasc
    r = fgasc%dim /= 0
  end function fgassoc_built_q

  pure function syncassoc_built_q(snasc) result(r)
    logical :: r
    type(syncassoc), intent(in) :: snasc
    r = snasc%dim /= 0
  end function syncassoc_built_q

  subroutine boxarray_bndry_periodic(bxai, dmn, b, nodal, pmask, ng, shfts, cross, idim)
    type(boxarray), intent(out) :: bxai
    type(box),      intent(in)  :: dmn, b
    logical,        intent(in)  :: nodal(:), pmask(:)
    integer,        intent(in)  :: ng
    integer,        intent(out) :: shfts(:,:)
    logical,        intent(in)  :: cross
    integer,        intent(in)  :: idim

    integer               :: i, cnt
    type(box)             :: bxs(3**get_dim(b)), gbx
    type(box),allocatable :: bv(:)
    integer               :: shft(3**get_dim(b),get_dim(b))
    integer               :: upbx(1:get_dim(b)), lwbx(1:get_dim(b))
    type(boxarray)        :: tba
    type(list_box)        :: bl

    if ( cross ) then
       gbx = box_nodalize(b,nodal)
       do i = 1, gbx%dim
          if (idim .ne. 0 .and. idim .ne. i) cycle
          !
          ! lo face
          !
          upbx    = upb(gbx)
          lwbx    = lwb(gbx)
          upbx(i) = lwbx(i) - 1
          lwbx(i) = lwbx(i) - ng
          call push_back(bl, make_box(lwbx,upbx))
          !
          ! hi face
          !
          upbx    = upb(gbx)
          lwbx    = lwb(gbx)
          lwbx(i) = upbx(i) + 1
          upbx(i) = upbx(i) + ng
          call push_back(bl, make_box(lwbx,upbx))
       end do
       call build(tba, bl, sort = .false.)
       call destroy(bl)
    else
       call boxarray_box_boundary_n(tba, box_nodalize(b,nodal), ng)       
    end if

    shfts = 0

    call box_periodic_shift(dmn, b, nodal, pmask, ng, shft, cnt, bxs)

    if ( cnt > 0 ) then
       allocate(bv(nboxes(tba)+cnt))
       do i = 1, nboxes(tba)
          bv(i) = get_box(tba,i)
       end do
       bv(nboxes(tba)+1:nboxes(tba)+cnt) = bxs(1:cnt)
       shfts(nboxes(tba)+1:nboxes(tba)+cnt,:) = shft(1:cnt,:)
       call destroy(tba)
       call boxarray_build_v(tba, bv, sort = .false.)
    end if

    bxai = tba

  end subroutine boxarray_bndry_periodic

  subroutine sumassoc_build_innards(bxasc, la, lap, latmp, ng, anynodal, lcnt_r, cnt_s, cnt_r, pvol, parr)
    type(boxassoc),   intent(inout) :: bxasc
    type(layout),     intent(inout) :: la
    type(layout_rep), intent(in)    :: lap
    type(layout),     intent(inout) :: latmp
    integer,          intent(in)    :: ng
    logical,          intent(in)    :: anynodal
    integer,          intent(inout) :: lcnt_r, cnt_s, cnt_r
    integer,          intent(inout) :: pvol(0:,:), parr(0:,:)

    integer                         :: i, j, ii, jj, i_r, i_s, li_r
    type(boxarray)                  :: bxa, bxai
    type(box)                       :: sbx, dbx
    type(local_copy_desc), pointer  :: n_cpy(:) => Null()
    type(comm_dsc), pointer         :: n_snd(:) => Null(), n_rcv(:) => Null()
    type(box_intersector), pointer  :: bi(:)
    integer                         :: shft(2*3**la%lap%dim,la%lap%dim), sh(MAX_SPACEDIM+1)

    bxa = get_boxarray(la)

    li_r = 1; i_r = 1; i_s = 1
    !
    ! Consider all copies I <- J.
    !
    do i = 1, nboxes(bxa)
       call boxarray_bndry_periodic(bxai, lap%pd, get_box(bxa,i), bxasc%nodal, lap%pmask, ng, shft, .false., bxasc%idim)
       do ii = 1, nboxes(bxai)
          if ( anynodal ) then
             bi => layout_get_box_intersector(latmp, get_box(bxai,ii))
          else
             bi => layout_get_box_intersector(la,    get_box(bxai,ii))
          end if
          do jj = 1, size(bi)
             j   = bi(jj)%i
             if ( remote(la,i) .and. remote(la,j) ) cycle
             dbx = bi(jj)%bx
             sbx = shift(dbx,-shft(ii,:))
             if ( local(la,i) .and. local(la, j) ) then
                if ( li_r > size(bxasc%l_con%cpy) ) then
                   allocate(n_cpy(size(bxasc%l_con%cpy) * 2))
                   n_cpy(1:li_r-1) = bxasc%l_con%cpy(1:li_r-1)
                   deallocate(bxasc%l_con%cpy)
                   bxasc%l_con%cpy => n_cpy
                end if
                lcnt_r                    = lcnt_r + 1
                bxasc%l_con%cpy(li_r)%nd  = j
                bxasc%l_con%cpy(li_r)%ns  = i
                bxasc%l_con%cpy(li_r)%sbx = sbx
                bxasc%l_con%cpy(li_r)%dbx = dbx
                li_r                      = li_r + 1
             else if ( local(la, j) ) then
                if ( i_r > size(bxasc%r_con%rcv) ) then
                   allocate(n_rcv(size(bxasc%r_con%rcv) * 2))
                   n_rcv(1:i_r-1) = bxasc%r_con%rcv(1:i_r-1)
                   deallocate(bxasc%r_con%rcv)
                   bxasc%r_con%rcv => n_rcv
                end if
                cnt_r                    = cnt_r + 1
                parr(lap%prc(i), 1)      = parr(lap%prc(i), 1) + 1
                pvol(lap%prc(i), 1)      = pvol(lap%prc(i), 1) + volume(dbx)
                bxasc%r_con%rcv(i_r)%nd  = j
                bxasc%r_con%rcv(i_r)%ns  = i
                bxasc%r_con%rcv(i_r)%sbx = sbx
                bxasc%r_con%rcv(i_r)%dbx = dbx
                bxasc%r_con%rcv(i_r)%pr  = get_proc(la, i)
                sh                       = 1
                sh(1:bxasc%dim)          = extent(dbx)
                bxasc%r_con%rcv(i_r)%sh  = sh
                i_r                      = i_r + 1
             else if ( local(la, i) ) then
                if ( i_s > size(bxasc%r_con%snd) ) then
                   allocate(n_snd(size(bxasc%r_con%snd) * 2))
                   n_snd(1:i_s-1) = bxasc%r_con%snd(1:i_s-1)
                   deallocate(bxasc%r_con%snd)
                   bxasc%r_con%snd => n_snd
                end if
                cnt_s                    = cnt_s + 1
                parr(lap%prc(j), 2)      = parr(lap%prc(j), 2) + 1
                pvol(lap%prc(j), 2)      = pvol(lap%prc(j), 2) + volume(dbx)
                bxasc%r_con%snd(i_s)%nd  = j
                bxasc%r_con%snd(i_s)%ns  = i
                bxasc%r_con%snd(i_s)%sbx = sbx
                bxasc%r_con%snd(i_s)%dbx = dbx
                bxasc%r_con%snd(i_s)%pr  = get_proc(la, j)
                i_s                      = i_s + 1
             end if
          end do
          deallocate(bi)
       end do
       call boxarray_destroy(bxai)
    end do
  end subroutine sumassoc_build_innards

  subroutine boxassoc_build_innards(bxasc, la, lap, latmp, ng, anynodal, cross, lcnt_r, cnt_s, cnt_r, pvol, parr)
    type(boxassoc),   intent(inout) :: bxasc
    type(layout),     intent(inout) :: la
    type(layout_rep), intent(in)    :: lap
    type(layout),     intent(inout) :: latmp
    integer,          intent(in)    :: ng
    logical,          intent(in)    :: anynodal
    logical,          intent(in)    :: cross
    integer,          intent(inout) :: lcnt_r, cnt_s, cnt_r
    integer,          intent(inout) :: pvol(0:,:), parr(0:,:)

    integer                         :: i, j, ii, jj, i_r, i_s, li_r
    type(boxarray)                  :: bxa, bxai
    type(box)                       :: sbx, dbx
    type(local_copy_desc), pointer  :: n_cpy(:) => Null()
    type(comm_dsc), pointer         :: n_snd(:) => Null(), n_rcv(:) => Null()
    type(box_intersector), pointer  :: bi(:)
    integer                         :: shft(2*3**la%lap%dim,la%lap%dim), sh(MAX_SPACEDIM+1)

    bxa = get_boxarray(la)

    li_r = 1; i_r = 1; i_s = 1
    !
    ! Consider all copies I <- J.
    !
    do i = 1, nboxes(bxa)
       call boxarray_bndry_periodic(bxai, lap%pd, get_box(bxa,i), bxasc%nodal, lap%pmask, ng, shft, cross, bxasc%idim)
       do ii = 1, nboxes(bxai)
          if ( anynodal ) then
             bi => layout_get_box_intersector(latmp, get_box(bxai,ii))
          else
             bi => layout_get_box_intersector(la,    get_box(bxai,ii))
          end if
          do jj = 1, size(bi)
             j   = bi(jj)%i
             if ( remote(la,i) .and. remote(la,j) ) cycle
             sbx = bi(jj)%bx
             dbx = shift(sbx,-shft(ii,:))
             if ( local(la,i) .and. local(la, j) ) then
                if ( li_r > size(bxasc%l_con%cpy) ) then
                   allocate(n_cpy(size(bxasc%l_con%cpy) * 2))
                   n_cpy(1:li_r-1) = bxasc%l_con%cpy(1:li_r-1)
                   deallocate(bxasc%l_con%cpy)
                   bxasc%l_con%cpy => n_cpy
                end if
                lcnt_r                    = lcnt_r + 1
                bxasc%l_con%cpy(li_r)%nd  = i
                bxasc%l_con%cpy(li_r)%ns  = j
                bxasc%l_con%cpy(li_r)%sbx = sbx
                bxasc%l_con%cpy(li_r)%dbx = dbx
                li_r                      = li_r + 1
             else if ( local(la, j) ) then
                if ( i_s > size(bxasc%r_con%snd) ) then
                   allocate(n_snd(size(bxasc%r_con%snd) * 2))
                   n_snd(1:i_s-1) = bxasc%r_con%snd(1:i_s-1)
                   deallocate(bxasc%r_con%snd)
                   bxasc%r_con%snd => n_snd
                end if
                cnt_s                    = cnt_s + 1
                parr(lap%prc(i), 2)      = parr(lap%prc(i), 2) + 1
                pvol(lap%prc(i), 2)      = pvol(lap%prc(i), 2) + volume(sbx)
                bxasc%r_con%snd(i_s)%nd  = i
                bxasc%r_con%snd(i_s)%ns  = j
                bxasc%r_con%snd(i_s)%sbx = sbx
                bxasc%r_con%snd(i_s)%dbx = dbx
                bxasc%r_con%snd(i_s)%pr  = get_proc(la, i)
                i_s                      = i_s + 1
             else if ( local(la, i) ) then
                if ( i_r > size(bxasc%r_con%rcv) ) then
                   allocate(n_rcv(size(bxasc%r_con%rcv) * 2))
                   n_rcv(1:i_r-1) = bxasc%r_con%rcv(1:i_r-1)
                   deallocate(bxasc%r_con%rcv)
                   bxasc%r_con%rcv => n_rcv
                end if
                cnt_r                    = cnt_r + 1
                parr(lap%prc(j), 1)      = parr(lap%prc(j), 1) + 1
                pvol(lap%prc(j), 1)      = pvol(lap%prc(j), 1) + volume(sbx)
                bxasc%r_con%rcv(i_r)%nd  = i
                bxasc%r_con%rcv(i_r)%ns  = j
                bxasc%r_con%rcv(i_r)%sbx = sbx
                bxasc%r_con%rcv(i_r)%dbx = dbx
                bxasc%r_con%rcv(i_r)%pr  = get_proc(la, j)
                sh                       = 1
                sh(1:bxasc%dim)          = extent(sbx)
                bxasc%r_con%rcv(i_r)%sh  = sh
                i_r                      = i_r + 1
             end if
          end do
          deallocate(bi)
       end do
       call boxarray_destroy(bxai)
    end do
  end subroutine boxassoc_build_innards

  subroutine boxassoc_build(bxasc, lap, ng, nodal, cross, do_sum_boundary, idim)
    use bl_prof_module
    use bl_error_module

    integer,          intent(in)           :: ng
    logical,          intent(in)           :: nodal(:)
    type(layout_rep), intent(in), target   :: lap
    type(boxassoc),   intent(inout)        :: bxasc
    logical,          intent(in)           :: cross
    logical,          intent(in), optional :: do_sum_boundary
    integer,          intent(in), optional :: idim

    integer                         :: pv, rpv, spv, pi_r, pi_s, pcnt_r, pcnt_s
    type(boxarray)                  :: bxa, batmp
    type(layout)                    :: la, latmp
    integer                         :: lcnt_r, cnt_r, cnt_s, i_r, i_s, np, i
    integer, parameter              :: ChunkSize = 500
    integer, allocatable            :: pvol(:,:), ppvol(:,:), parr(:,:)
    type(local_copy_desc), pointer  :: n_cpy(:) => Null()
    type(comm_dsc), pointer         :: n_snd(:) => Null(), n_rcv(:) => Null()
    logical                         :: anynodal, ldsm
    type(bl_prof_timer), save       :: bpt

    if ( built_q(bxasc) ) call bl_error("boxassoc_build(): already built")

    call build(bpt, "boxassoc_build")

    la%lap       => lap
    bxa          =  get_boxarray(la)
    bxasc%dim    =  get_dim(bxa)
    bxasc%grwth  =  ng
    bxasc%nboxes =  nboxes(bxa)
    bxasc%cross  =  cross
    np           =  parallel_nprocs()
    anynodal     =  any( nodal .eqv. .true. )

    ldsm = .false.; if ( present(do_sum_boundary) ) ldsm = do_sum_boundary

    bxasc%idim = 0; if ( present(idim) ) bxasc%idim = idim
    if (bxasc%idim < 0 .or. bxasc%idim > bxasc%dim) then
       call bl_error("BOXASSOC_BUILD:", bxasc%idim)
    end if
    
    if (bxasc%idim .ne. 0) then
       bxasc%cross = .true.   ! idim > 0 is a special case of cross
    end if

    allocate(bxasc%nodal(bxasc%dim))
    allocate(parr(0:np-1,2))
    allocate(pvol(0:np-1,2))
    allocate(ppvol(0:np-1,2))
    allocate(bxasc%l_con%cpy(ChunkSize))
    allocate(bxasc%r_con%snd(ChunkSize))
    allocate(bxasc%r_con%rcv(ChunkSize))

    if ( anynodal ) then
       !
       ! Build a temporary layout to be used in intersection tests below.
       !
       call copy(batmp, bxa)
       call boxarray_nodalize(batmp, nodal)
       call build(latmp, batmp, boxarray_bbox(batmp), mapping = LA_LOCAL)  ! LA_LOCAL ==> bypass processor distribution calculation.
       call destroy(batmp)
    end if

    bxasc%nodal = nodal

    parr = 0; pvol = 0; lcnt_r = 0; cnt_r = 0; cnt_s = 0

    if ( ldsm ) then
       if (bxasc%idim .ne. 0) then
          call bl_error("BOXASSOC_BUILD: sumassoc_build_innards not supported for nonzero idim")
       end if
       call sumassoc_build_innards(bxasc, la, lap, latmp, ng, anynodal, &
            lcnt_r, cnt_s, cnt_r, pvol, parr)
    else
       call boxassoc_build_innards(bxasc, la, lap, latmp, ng, anynodal, bxasc%cross, &
            lcnt_r, cnt_s, cnt_r, pvol, parr)
    end if

    if ( anynodal ) then
       call destroy(latmp)
    else
       call clear_box_hash_bin(la%lap)
    end if

    bxasc%l_con%ncpy = lcnt_r
    bxasc%r_con%nsnd = cnt_s
    bxasc%r_con%nrcv = cnt_r
    !
    ! Trim off unused space.
    !
    allocate(n_cpy(lcnt_r))
    n_cpy(1:lcnt_r) = bxasc%l_con%cpy(1:lcnt_r)
    deallocate(bxasc%l_con%cpy)
    bxasc%l_con%cpy => n_cpy

    allocate(n_snd(cnt_s))
    n_snd(1:cnt_s) = bxasc%r_con%snd(1:cnt_s)
    deallocate(bxasc%r_con%snd)
    bxasc%r_con%snd => n_snd

    allocate(n_rcv(cnt_r))
    n_rcv(1:cnt_r) = bxasc%r_con%rcv(1:cnt_r)
    deallocate(bxasc%r_con%rcv)
    bxasc%r_con%rcv => n_rcv
    !
    ! This region packs the src/recv boxes into processor order
    !
    do i = 0, np-1
       ppvol(i,1) = sum(pvol(0:i-1,1))
       ppvol(i,2) = sum(pvol(0:i-1,2))
    end do
    !
    ! Pack Receives maintaining original ordering
    !
    do i_r = 1, cnt_r
       i = bxasc%r_con%rcv(i_r)%pr
       bxasc%r_con%rcv(i_r)%pv = ppvol(i,1)
       pv = volume(bxasc%r_con%rcv(i_r)%dbx)
       ppvol(i,1) = ppvol(i,1) + pv
    end do
    !
    ! Pack Sends maintaining original ordering
    !
    do i_s = 1, cnt_s
       i = bxasc%r_con%snd(i_s)%pr
       bxasc%r_con%snd(i_s)%pv = ppvol(i,2)
       pv = volume(bxasc%r_con%snd(i_s)%dbx)
       ppvol(i,2) = ppvol(i,2) + pv
    end do
    !
    ! Now compute the volume of data that each processor expects.
    !
    pcnt_r = count(parr(:,1) /= 0 )
    pcnt_s = count(parr(:,2) /= 0 )
    bxasc%r_con%nrp  = pcnt_r
    bxasc%r_con%nsp  = pcnt_s
    bxasc%r_con%rvol = sum(pvol(:,1))
    bxasc%r_con%svol = sum(pvol(:,2))
    allocate(bxasc%r_con%str(pcnt_s))
    allocate(bxasc%r_con%rtr(pcnt_r))
    pi_r = 1; pi_s = 1; rpv  = 0; spv  = 0
    do i = 0, size(pvol,dim=1)-1
       if ( pvol(i,1) /= 0 ) then
          bxasc%r_con%rtr(pi_r)%sz = pvol(i,1)
          bxasc%r_con%rtr(pi_r)%pr = i
          bxasc%r_con%rtr(pi_r)%pv = rpv
          rpv  = rpv + pvol(i,1)
          pi_r = pi_r + 1
       end if
       if ( pvol(i,2) /= 0 ) then
          bxasc%r_con%str(pi_s)%sz = pvol(i,2)
          bxasc%r_con%str(pi_s)%pr = i
          bxasc%r_con%str(pi_s)%pv = spv
          spv  = spv + pvol(i,2)
          pi_s = pi_s + 1
       end if
    end do

    call mem_stats_alloc(bxa_ms, bxasc%r_con%nsnd + bxasc%r_con%nrcv)

    call destroy(bpt)

  end subroutine boxassoc_build

  subroutine fgassoc_build(fgasc, la, ng)
    use bl_prof_module
    use bl_error_module

    integer,       intent(in   ) :: ng
    type(layout),  intent(inout) :: la     ! Only modified by layout_get_box_intersector()
    type(fgassoc), intent(inout) :: fgasc

    integer                        :: i, j
    type(box)                      :: bx
    type(list_box)                 :: bl, pieces, leftover
    type(box_intersector), pointer :: bi(:)
    type(bl_prof_timer), save      :: bpt

    if ( built_q(fgasc) ) call bl_error("fgassoc_build(): already built")

    call build(bpt, "fgassoc_build")

    fgasc%dim   = la%lap%dim
    fgasc%grwth = ng
    !
    ! Build list of all ghost cells not covered by valid region.
    !
    do i = 1, nboxes(la%lap%bxa)
       bx = get_box(la%lap%bxa,i)
       call boxarray_box_diff(fgasc%ba, grow(bx,ng), bx)
       do j = 1, nboxes(fgasc%ba)
          call push_back(bl, get_box(fgasc%ba,j))
       end do
       call boxarray_destroy(fgasc%ba)
    end do

    call boxarray_build_l(fgasc%ba, bl, sort = .false.)

    call list_destroy_box(bl)

    do i = 1, nboxes(fgasc%ba)
       bx = get_box(fgasc%ba,i)
       bi => layout_get_box_intersector(la, bx)
       do j = 1, size(bi)
          call push_back(pieces, bi(j)%bx)
       end do
       deallocate(bi)
       leftover = boxlist_boxlist_diff(bx, pieces)
       call splice(bl, leftover)
       call list_destroy_box(pieces)
    end do

    call clear_box_hash_bin(la%lap)
    !
    ! Remove any overlaps on remaining cells.
    !
    call boxarray_destroy(fgasc%ba)
    call boxarray_build_l(fgasc%ba, bl, sort = .false.)
    call list_destroy_box(bl)
    call boxarray_to_domain(fgasc%ba)
    !
    ! Chop things up a bit so no one processor gets too much data.
    !
    select case ( la%lap%dim ) 
    case (3)
       call boxarray_maxsize(fgasc%ba, 32)
    case (2)
       call boxarray_maxsize(fgasc%ba, 128)
    end select

    call mem_stats_alloc(fgx_ms, volume(fgasc%ba))
    
    call destroy(bpt)

  end subroutine fgassoc_build

  subroutine internal_sync_unique_cover(la, ng, nodal, lall, filled)
    use bl_error_module

    type(layout), intent(in)  :: la
    integer, intent(in)       :: ng
    logical, intent(in)       :: nodal(:)
    logical, intent(in)       :: lall
    type(local_conn), pointer :: filled(:)

    type(box)                      :: jbx, abx
    integer                        :: i, j, k, ii, jj, cnt
    integer                        :: shft(2*3**(la%lap%dim),la%lap%dim)
    integer, parameter             :: ChunkSize = 500
    type(local_copy_desc)          :: lcd
    type(local_copy_desc), pointer :: n_cpy(:) => Null()
    type(list_box)                 :: lb1, lb2
    type(boxarray)                 :: bxa, batmp
    type(box_intersector), pointer :: bi(:)
    type(layout)                   :: latmp

    bxa = get_boxarray(la)
    !
    ! Build a temporary layout to be used in intersection tests below.
    !
    call copy(batmp, bxa)
    call boxarray_nodalize(batmp, nodal)
    if ( lall ) call boxarray_grow(batmp,ng)
    call build(latmp, batmp, boxarray_bbox(batmp), mapping = LA_LOCAL)  ! LA_LOCAL ==> bypass processor distribution calculation.
    call destroy(batmp)

    allocate(filled(nboxes(bxa)))
    do i = 1, nboxes(bxa)
       filled(i)%ncpy = 0
       allocate(filled(i)%cpy(ChunkSize))
    end do

    do j = 1, nboxes(bxa)
       jbx = get_box(latmp, j)
       call box_internal_sync_shift(la%lap%pd, jbx, la%lap%pmask, nodal, shft, cnt)
       do jj = 1, cnt
          bi => layout_get_box_intersector(latmp, shift(jbx,shft(jj,:)))
          do ii = 1, size(bi)
             !
             ! Only consider j -> i for i >= j
             !
             if ( bi(ii)%i < j ) cycle
             i   = bi(ii)%i
             abx = bi(ii)%bx
             !
             ! Do not overwrite ourselves.
             !
             if ( i == j .and. all(shft(jj,:) == 0) ) cycle
             !
             ! Find parts of abx that haven't been written to already.
             !
             do k = 1, filled(i)%ncpy
                call push_back(lb1, filled(i)%cpy(k)%dbx)
             end do
             lb2 = boxlist_boxlist_diff(abx, lb1)
             do while ( .not. empty(lb2) )
                filled(i)%ncpy = filled(i)%ncpy + 1
                if ( filled(i)%ncpy > size(filled(i)%cpy) ) then
                   allocate(n_cpy(size(filled(i)%cpy) * 2))
                   n_cpy(1:filled(i)%ncpy-1) = filled(i)%cpy(1:filled(i)%ncpy-1)
                   deallocate(filled(i)%cpy)
                   filled(i)%cpy => n_cpy
                end if
                lcd%ns  = j
                lcd%nd  = i
                lcd%sbx = shift(front(lb2), -shft(jj,:))
                lcd%dbx = front(lb2)
                filled(i)%cpy(filled(i)%ncpy) = lcd
                call pop_front(lb2)
             end do
             call destroy(lb1)
             call destroy(lb2)
          end do
          deallocate(bi)
       end do
    end do
    call destroy(latmp)

  end subroutine internal_sync_unique_cover

  subroutine syncassoc_build(snasc, lap, ng, nodal, lall)
    use bl_prof_module
    use bl_error_module

    integer,          intent(in)         :: ng
    logical,          intent(in)         :: nodal(:)
    type(layout_rep), intent(in), target :: lap
    type(syncassoc),  intent(inout)      :: snasc
    logical,          intent(in)         :: lall

    integer                        :: i, j, ii, jj, pv, rpv, spv, pi_r, pi_s, pcnt_r, pcnt_s, np
    type(box)                      :: dbx, sbx
    type(boxarray)                 :: bxa
    type(layout)                   :: la
    integer                        :: lcnt_r, li_r, cnt_r, cnt_s, i_r, i_s, sh(MAX_SPACEDIM+1)
    integer, parameter             :: ChunkSize = 500
    integer, allocatable           :: pvol(:,:), ppvol(:,:), parr(:,:)
    type(local_copy_desc), pointer :: n_cpy(:) => Null()
    type(comm_dsc), pointer        :: n_snd(:) => Null(), n_rcv(:) => Null()
    type(local_conn), pointer      :: filled(:)
    type(bl_prof_timer), save      :: bpt

    if ( built_q(snasc) ) call bl_error("syncassoc_build(): already built")

    call build(bpt, "syncassoc_build")

    la%lap       => lap
    bxa          =  get_boxarray(la)
    snasc%dim    =  get_dim(bxa)
    snasc%grwth  =  ng
    snasc%nboxes =  nboxes(bxa)
    np           =  parallel_nprocs()

    allocate(snasc%nodal(snasc%dim))
    allocate(parr(0:np-1,2))
    allocate(pvol(0:np-1,2))
    allocate(ppvol(0:np-1,2))
    allocate(snasc%l_con%cpy(ChunkSize))
    allocate(snasc%r_con%snd(ChunkSize))
    allocate(snasc%r_con%rcv(ChunkSize))

    snasc%lall  = lall
    snasc%nodal = nodal

    call internal_sync_unique_cover(la, snasc%grwth, snasc%nodal, snasc%lall, filled)

    parr = 0; pvol = 0; lcnt_r = 0; cnt_r = 0; cnt_s = 0; li_r = 1; i_r = 1; i_s = 1

    do jj = 1, nboxes(bxa)
       if ( filled(jj)%ncpy > 0 ) then
          do ii = 1, filled(jj)%ncpy
             i   = filled(jj)%cpy(ii)%nd
             j   = filled(jj)%cpy(ii)%ns
             sbx = filled(jj)%cpy(ii)%sbx
             dbx = filled(jj)%cpy(ii)%dbx
             if ( local(la, i) .and. local(la, j) ) then
                if ( li_r > size(snasc%l_con%cpy) ) then
                   allocate(n_cpy(size(snasc%l_con%cpy) * 2))
                   n_cpy(1:li_r-1) = snasc%l_con%cpy(1:li_r-1)
                   deallocate(snasc%l_con%cpy)
                   snasc%l_con%cpy => n_cpy
                end if
                lcnt_r                    = lcnt_r + 1
                snasc%l_con%cpy(li_r)%nd  = i
                snasc%l_con%cpy(li_r)%ns  = j
                snasc%l_con%cpy(li_r)%sbx = sbx
                snasc%l_con%cpy(li_r)%dbx = dbx
                li_r                      = li_r + 1
             else if ( local(la, j) ) then ! must send
                if ( i_s > size(snasc%r_con%snd) ) then
                   allocate(n_snd(size(snasc%r_con%snd) * 2))
                   n_snd(1:i_s-1) = snasc%r_con%snd(1:i_s-1)
                   deallocate(snasc%r_con%snd)
                   snasc%r_con%snd => n_snd
                end if
                cnt_s                    = cnt_s + 1
                parr(lap%prc(i), 2)      = parr(lap%prc(i), 2) + 1
                pvol(lap%prc(i), 2)      = pvol(lap%prc(i), 2) + volume(dbx)
                snasc%r_con%snd(i_s)%nd  = i
                snasc%r_con%snd(i_s)%ns  = j
                snasc%r_con%snd(i_s)%sbx = sbx
                snasc%r_con%snd(i_s)%dbx = dbx
                snasc%r_con%snd(i_s)%pr  = get_proc(la, i)
                i_s                      = i_s + 1
             else if ( local(la, i) ) then  ! must recv
                if ( i_r > size(snasc%r_con%rcv) ) then
                   allocate(n_rcv(size(snasc%r_con%rcv) * 2))
                   n_rcv(1:i_r-1) = snasc%r_con%rcv(1:i_r-1)
                   deallocate(snasc%r_con%rcv)
                   snasc%r_con%rcv => n_rcv
                end if
                cnt_r                    = cnt_r + 1
                parr(lap%prc(j), 1)      = parr(lap%prc(j), 1) + 1
                pvol(lap%prc(j), 1)      = pvol(lap%prc(j), 1) + volume(dbx)
                snasc%r_con%rcv(i_r)%nd  = i
                snasc%r_con%rcv(i_r)%ns  = j
                snasc%r_con%rcv(i_r)%sbx = sbx
                snasc%r_con%rcv(i_r)%dbx = dbx
                snasc%r_con%rcv(i_r)%pr  = get_proc(la, j)
                sh                       = 1
                sh(1:snasc%dim)          = extent(dbx)
                snasc%r_con%rcv(i_r)%sh  = sh
                i_r                      = i_r + 1
             end if
          end do
       end if
    end do

    do i = 1, nboxes(bxa)
       deallocate(filled(i)%cpy)
    end do
    deallocate(filled)

    snasc%l_con%ncpy = lcnt_r
    snasc%r_con%nsnd = cnt_s
    snasc%r_con%nrcv = cnt_r
    !
    ! Trim off unused space.
    !
    allocate(n_cpy(lcnt_r))
    n_cpy(1:lcnt_r) = snasc%l_con%cpy(1:lcnt_r)
    deallocate(snasc%l_con%cpy)
    snasc%l_con%cpy => n_cpy

    allocate(n_snd(cnt_s))
    n_snd(1:cnt_s) = snasc%r_con%snd(1:cnt_s)
    deallocate(snasc%r_con%snd)
    snasc%r_con%snd => n_snd

    allocate(n_rcv(cnt_r))
    n_rcv(1:cnt_r) = snasc%r_con%rcv(1:cnt_r)
    deallocate(snasc%r_con%rcv)
    snasc%r_con%rcv => n_rcv
    !
    ! This region packs the src/recv boxes into processor order
    !
    do i = 0, np-1
       ppvol(i,1) = sum(pvol(0:i-1,1))
       ppvol(i,2) = sum(pvol(0:i-1,2))
    end do
    !
    ! Pack Receives maintaining original ordering
    !
    do i_r = 1, cnt_r
       i = snasc%r_con%rcv(i_r)%pr
       snasc%r_con%rcv(i_r)%pv = ppvol(i,1)
       pv = volume(snasc%r_con%rcv(i_r)%dbx)
       ppvol(i,1) = ppvol(i,1) + pv
    end do
    !
    ! Pack Sends maintaining original ordering
    !
    do i_s = 1, cnt_s
       i = snasc%r_con%snd(i_s)%pr
       snasc%r_con%snd(i_s)%pv = ppvol(i,2)
       pv = volume(snasc%r_con%snd(i_s)%dbx)
       ppvol(i,2) = ppvol(i,2) + pv
    end do
    !
    ! Now compute the volume of data the each processor expects
    !
    pcnt_r = count(parr(:,1) /= 0 )
    pcnt_s = count(parr(:,2) /= 0 )
    snasc%r_con%nrp  = pcnt_r
    snasc%r_con%nsp  = pcnt_s
    snasc%r_con%rvol = sum(pvol(:,1))
    snasc%r_con%svol = sum(pvol(:,2))
    allocate(snasc%r_con%str(pcnt_s))
    allocate(snasc%r_con%rtr(pcnt_r))
    pi_r = 1; pi_s = 1; rpv  = 0; spv  = 0
    do i = 0, size(pvol,dim=1)-1
       if ( pvol(i,1) /= 0 ) then
          snasc%r_con%rtr(pi_r)%sz = pvol(i,1)
          snasc%r_con%rtr(pi_r)%pr = i
          snasc%r_con%rtr(pi_r)%pv = rpv
          rpv  = rpv + pvol(i,1)
          pi_r = pi_r + 1
       end if
       if ( pvol(i,2) /= 0 ) then
          snasc%r_con%str(pi_s)%sz = pvol(i,2)
          snasc%r_con%str(pi_s)%pr = i
          snasc%r_con%str(pi_s)%pv = spv
          spv  = spv + pvol(i,2)
          pi_s = pi_s + 1
       end if
    end do

    call mem_stats_alloc(snx_ms, snasc%r_con%nsnd + snasc%r_con%nrcv)

    call destroy(bpt)

  end subroutine syncassoc_build

  function boxassoc_bytes(bxasc) result(nbytes)
    use bl_error_module
    integer :: nbytes
    type(boxassoc), intent(in) :: bxasc
    if ( .not. built_q(bxasc) ) call bl_error("boxassoc_bytes(): not built")

    nbytes = 0
    nbytes = nbytes + 4 * 9;
    nbytes = nbytes + bxasc%l_con%ncpy * 4 * 16
    nbytes = nbytes + bxasc%r_con%nsnd * 4 * 25
    nbytes = nbytes + bxasc%r_con%nrcv * 4 * 25
    nbytes = nbytes + bxasc%r_con%nsp  * 4 * 3
    nbytes = nbytes + bxasc%r_con%nrp  * 4 * 3

  end function boxassoc_bytes

  subroutine boxassoc_destroy(bxasc)
    use bl_error_module
    type(boxassoc), intent(inout) :: bxasc
    if ( .not. built_q(bxasc) ) call bl_error("boxassoc_destroy(): not built")
    call mem_stats_dealloc(bxa_ms, bxasc%r_con%nsnd + bxasc%r_con%nrcv)
    deallocate(bxasc%nodal)
    deallocate(bxasc%l_con%cpy)
    deallocate(bxasc%r_con%snd)
    deallocate(bxasc%r_con%rcv)
    deallocate(bxasc%r_con%str)
    deallocate(bxasc%r_con%rtr)
  end subroutine boxassoc_destroy

  subroutine fgassoc_destroy(fgasc)
    use bl_error_module
    type(fgassoc), intent(inout) :: fgasc
    if ( .not. built_q(fgasc) ) call bl_error("fgassoc_destroy(): not built")
    call mem_stats_dealloc(fgx_ms, volume(fgasc%ba))
    call destroy(fgasc%ba)
  end subroutine fgassoc_destroy

  subroutine syncassoc_destroy(snasc)
    use bl_error_module
    type(syncassoc), intent(inout) :: snasc
    if ( .not. built_q(snasc) ) call bl_error("syncassoc_destroy(): not built")
    call mem_stats_dealloc(snx_ms, snasc%r_con%nsnd + snasc%r_con%nrcv)
    deallocate(snasc%nodal)
    deallocate(snasc%l_con%cpy)
    deallocate(snasc%r_con%snd)
    deallocate(snasc%r_con%rcv)
    deallocate(snasc%r_con%str)
    deallocate(snasc%r_con%rtr)
  end subroutine syncassoc_destroy

  subroutine copyassoc_build(cpasc, la_dst, la_src, nd_dst, nd_src)
    use bl_prof_module
    use bl_error_module

    type(copyassoc),  intent(inout) :: cpasc
    type(layout),     intent(in)    :: la_src, la_dst
    logical,          intent(in)    :: nd_dst(:), nd_src(:)

    integer                        :: i, j, pv, rpv, spv, pi_r, pi_s, pcnt_r, pcnt_s
    integer                        :: sh(MAX_SPACEDIM+1), jj
    type(box)                      :: bx
    type(boxarray)                 :: batmp
    type(layout)                   :: lasrctmp
    integer                        :: lcnt_r, li_r, cnt_r, cnt_s, i_r, i_s
    integer, allocatable           :: pvol(:,:), ppvol(:,:), parr(:,:)
    type(local_copy_desc), pointer :: n_cpy(:) => Null()
    type(comm_dsc), pointer        :: n_snd(:) => Null(), n_rcv(:) => Null()
    integer, parameter             :: ChunkSize = 500
    type(box_intersector), pointer :: bi(:)
    type(bl_prof_timer), save      :: bpt

    if ( built_q(cpasc) ) call bl_error("copyassoc_build(): already built")

    call build(bpt, "copyassoc_build")

    call boxarray_build_copy(cpasc%ba_src, get_boxarray(la_src))
    call boxarray_build_copy(cpasc%ba_dst, get_boxarray(la_dst))

    cpasc%dim  = get_dim(cpasc%ba_src)
    cpasc%hash = layout_boxarray_hash(cpasc%ba_src) + layout_boxarray_hash(cpasc%ba_dst)

    allocate(cpasc%prc_src(la_src%lap%nboxes))
    allocate(cpasc%prc_dst(la_dst%lap%nboxes))
    allocate(cpasc%nd_dst(la_dst%lap%dim))
    allocate(cpasc%nd_src(la_src%lap%dim))
    allocate(parr(0:parallel_nprocs()-1,2))
    allocate(pvol(0:parallel_nprocs()-1,2))
    allocate(ppvol(0:parallel_nprocs()-1,2))
    allocate(cpasc%l_con%cpy(ChunkSize))
    allocate(cpasc%r_con%snd(ChunkSize))
    allocate(cpasc%r_con%rcv(ChunkSize))

    cpasc%prc_src = la_src%lap%prc
    cpasc%prc_dst = la_dst%lap%prc

    cpasc%nd_dst = nd_dst
    cpasc%nd_src = nd_src
    !
    ! Build a temporary layout to be used in intersection tests below.
    !
    call copy(batmp, cpasc%ba_src)
    call boxarray_nodalize(batmp, nd_src)
    call build(lasrctmp, batmp, boxarray_bbox(batmp), mapping = LA_LOCAL)  ! LA_LOCAL ==> bypass processor distribution calculation.
    call destroy(batmp)

    parr = 0; pvol = 0; lcnt_r = 0; cnt_r = 0; cnt_s = 0; li_r = 1; i_r = 1; i_s = 1
    !
    ! Consider all copies I <- J.
    !
    do i = 1, nboxes(cpasc%ba_dst)
       bi => layout_get_box_intersector(lasrctmp, box_nodalize(get_box(cpasc%ba_dst,i), nd_dst))
       do jj = 1, size(bi)
          j = bi(jj)%i
          if ( remote(la_dst,i) .and. remote(la_src,j) ) cycle
          bx = bi(jj)%bx
          if ( local(la_dst, i) .and. local(la_src, j) ) then
             if ( li_r > size(cpasc%l_con%cpy) ) then
                allocate(n_cpy(size(cpasc%l_con%cpy) * 2))
                n_cpy(1:li_r-1) = cpasc%l_con%cpy(1:li_r-1)
                deallocate(cpasc%l_con%cpy)
                cpasc%l_con%cpy => n_cpy
             end if
             lcnt_r                    = lcnt_r + 1
             cpasc%l_con%cpy(li_r)%nd  = i
             cpasc%l_con%cpy(li_r)%ns  = j
             cpasc%l_con%cpy(li_r)%sbx = bx
             cpasc%l_con%cpy(li_r)%dbx = bx
             li_r                      = li_r + 1
          else if ( local(la_src, j) ) then
             if ( i_s > size(cpasc%r_con%snd) ) then
                allocate(n_snd(size(cpasc%r_con%snd) * 2))
                n_snd(1:i_s-1) = cpasc%r_con%snd(1:i_s-1)
                deallocate(cpasc%r_con%snd)
                cpasc%r_con%snd => n_snd
             end if
             cnt_s                      = cnt_s + 1
             parr(la_dst%lap%prc(i), 2) = parr(la_dst%lap%prc(i), 2) + 1
             pvol(la_dst%lap%prc(i), 2) = pvol(la_dst%lap%prc(i), 2) + volume(bx)
             cpasc%r_con%snd(i_s)%nd    = i
             cpasc%r_con%snd(i_s)%ns    = j
             cpasc%r_con%snd(i_s)%sbx   = bx
             cpasc%r_con%snd(i_s)%dbx   = bx
             cpasc%r_con%snd(i_s)%pr    = get_proc(la_dst,i)
             i_s                        = i_s + 1
          else if ( local(la_dst, i) ) then
             if ( i_r > size(cpasc%r_con%rcv) ) then
                allocate(n_rcv(size(cpasc%r_con%rcv) * 2))
                n_rcv(1:i_r-1) = cpasc%r_con%rcv(1:i_r-1)
                deallocate(cpasc%r_con%rcv)
                cpasc%r_con%rcv => n_rcv
             end if
             cnt_r                      = cnt_r + 1
             parr(la_src%lap%prc(j), 1) = parr(la_src%lap%prc(j), 1) + 1
             pvol(la_src%lap%prc(j), 1) = pvol(la_src%lap%prc(j), 1) + volume(bx)
             cpasc%r_con%rcv(i_r)%nd    = i
             cpasc%r_con%rcv(i_r)%ns    = j
             cpasc%r_con%rcv(i_r)%sbx   = bx
             cpasc%r_con%rcv(i_r)%dbx   = bx
             cpasc%r_con%rcv(i_r)%pr    = get_proc(la_src,j)
             sh                         = 1
             sh(1:cpasc%dim)            = extent(bx)
             cpasc%r_con%rcv(i_r)%sh    = sh
             i_r                        = i_r + 1
          end if
       end do
       deallocate(bi)
    end do

    call destroy(lasrctmp)

    cpasc%l_con%ncpy = lcnt_r
    cpasc%r_con%nsnd = cnt_s
    cpasc%r_con%nrcv = cnt_r
    !
    ! Trim off unused space.
    !
    allocate(n_cpy(lcnt_r))
    n_cpy(1:lcnt_r) = cpasc%l_con%cpy(1:lcnt_r)
    deallocate(cpasc%l_con%cpy)
    cpasc%l_con%cpy => n_cpy

    allocate(n_snd(cnt_s))
    n_snd(1:cnt_s) = cpasc%r_con%snd(1:cnt_s)
    deallocate(cpasc%r_con%snd)
    cpasc%r_con%snd => n_snd

    allocate(n_rcv(cnt_r))
    n_rcv(1:cnt_r) = cpasc%r_con%rcv(1:cnt_r)
    deallocate(cpasc%r_con%rcv)
    cpasc%r_con%rcv => n_rcv
    !
    ! This region packs the src/recv boxes into processor order
    !
    do i = 0, parallel_nprocs()-1
       ppvol(i,1) = sum(pvol(0:i-1,1))
       ppvol(i,2) = sum(pvol(0:i-1,2))
    end do
    !
    ! Pack Receives maintaining original ordering
    !
    do i_r = 1, cnt_r
       i = cpasc%r_con%rcv(i_r)%pr
       cpasc%r_con%rcv(i_r)%pv = ppvol(i,1)
       pv = volume(cpasc%r_con%rcv(i_r)%dbx)
       ppvol(i,1) = ppvol(i,1) + pv
    end do
    !
    ! Pack Sends maintaining original ordering
    !
    do i_s = 1, cnt_s
       i = cpasc%r_con%snd(i_s)%pr
       cpasc%r_con%snd(i_s)%pv = ppvol(i,2)
       pv = volume(cpasc%r_con%snd(i_s)%dbx)
       ppvol(i,2) = ppvol(i,2) + pv
    end do
    !
    ! Now compute the volume of data the each processor expects
    !
    pcnt_r = count(parr(:,1) /= 0 )
    pcnt_s = count(parr(:,2) /= 0 )
    cpasc%r_con%nrp  = pcnt_r
    cpasc%r_con%nsp  = pcnt_s
    cpasc%r_con%rvol = sum(pvol(:,1))
    cpasc%r_con%svol = sum(pvol(:,2))
    allocate(cpasc%r_con%str(pcnt_s))
    allocate(cpasc%r_con%rtr(pcnt_r))
    pi_r = 1; pi_s = 1; rpv  = 0; spv  = 0
    do i = 0, size(pvol,dim=1)-1
       if ( pvol(i,1) /= 0 ) then
          cpasc%r_con%rtr(pi_r)%sz = pvol(i,1)
          cpasc%r_con%rtr(pi_r)%pr = i
          cpasc%r_con%rtr(pi_r)%pv = rpv
          rpv  = rpv + pvol(i,1)
          pi_r = pi_r + 1
       end if
       if ( pvol(i,2) /= 0 ) then
          cpasc%r_con%str(pi_s)%sz = pvol(i,2)
          cpasc%r_con%str(pi_s)%pr = i
          cpasc%r_con%str(pi_s)%pv = spv
          spv  = spv + pvol(i,2)
          pi_s = pi_s + 1
       end if
    end do

    call mem_stats_alloc(cpx_ms, cpasc%r_con%nsnd + cpasc%r_con%nrcv)

    call destroy(bpt)

  end subroutine copyassoc_build

  subroutine fluxassoc_build(flasc, la_dst, la_src, nd_dst, nd_src, side, crse_domain, ir)
    use bl_prof_module
    use bl_error_module

    type(fluxassoc),  intent(inout) :: flasc
    type(layout),     intent(in)    :: la_dst
    type(layout),     intent(inout) :: la_src
    logical,          intent(in)    :: nd_dst(:), nd_src(:)
    integer,          intent(in)    :: side
    type(box),        intent(in)    :: crse_domain
    integer,          intent(in)    :: ir(:)

    integer                        :: i, j, pv, rpv, spv, pi_r, pi_s, np, dir, dm
    integer                        :: sh(MAX_SPACEDIM+1)
    integer                        :: lo_dom(la_dst%lap%dim), hi_dom(la_dst%lap%dim), loflux(la_dst%lap%dim)
    type(box)                      :: fbox, isect
    type(layout)                   :: lasrctmp
    type(boxarray)                 :: bxa_src, bxa_dst, batmp
    integer                        :: lcnt_r, li_r, cnt_r, cnt_s, i_r, i_s, ii
    integer, allocatable           :: pvol(:,:), ppvol(:,:), parr(:,:), mpvol(:,:)
    type(local_copy_desc), pointer :: n_cpy(:) => Null()
    type(comm_dsc), pointer        :: n_snd(:) => Null(), n_rcv(:) => Null()
    type(box), pointer             :: pfbxs(:) => Null()
    type(box_intersector), pointer :: bi(:)
    integer, parameter             :: ChunkSize = 500
    logical                        :: anynodal
    type(bl_prof_timer), save      :: bpt

    if ( built_q(flasc) ) call bl_error("fluxassoc_build(): already built")

    call build(bpt, "fluxassoc_build")

    dm                =  la_dst%lap%dim
    np                =  parallel_nprocs()
    dir               =  iabs(side)
    lo_dom            =  lwb(crse_domain)
    hi_dom            =  upb(crse_domain)+1
    bxa_src           =  get_boxarray(la_src)
    bxa_dst           =  get_boxarray(la_dst)
    flasc%dim         =  dm
    flasc%side        =  side
    flasc%crse_domain =  crse_domain
    flasc%lap_dst     => la_dst%lap
    flasc%lap_src     => la_src%lap
    flasc%ir(1:dm)    =  ir(1:dm)
    flasc%flux%dim    =  dm
    flasc%mask%dim    =  dm
    anynodal          =  any( nd_src .eqv. .true. )

    allocate(flasc%nd_dst(dm), flasc%nd_src(dm))

    flasc%nd_dst = nd_dst
    flasc%nd_src = nd_src

    allocate(parr(0:np-1,2))
    allocate(pvol(0:np-1,2))
    allocate(mpvol(0:np-1,2))
    allocate(ppvol(0:np-1,2))
    allocate(flasc%flux%l_con%cpy(ChunkSize))
    allocate(flasc%flux%r_con%snd(ChunkSize))
    allocate(flasc%flux%r_con%rcv(ChunkSize))
    allocate(flasc%mask%r_con%snd(ChunkSize))
    allocate(flasc%mask%r_con%rcv(ChunkSize))
    allocate(flasc%fbxs(ChunkSize))

    if ( anynodal ) then
       !
       ! Build a temporary layout to be used in intersection tests below.
       !
       call copy(batmp, bxa_src)
       call boxarray_nodalize(batmp, nd_src)
       call build(lasrctmp, batmp, boxarray_bbox(batmp), mapping = LA_LOCAL)  ! LA_LOCAL ==> bypass processor distribution calculation.
       call destroy(batmp)
    end if

    parr = 0; pvol = 0; mpvol = 0; lcnt_r = 0; cnt_r = 0; cnt_s = 0; li_r = 1; i_r = 1; i_s = 1

    do j = 1, nboxes(bxa_dst)

       if ( anynodal ) then
          bi => layout_get_box_intersector(lasrctmp, box_nodalize(get_box(bxa_dst,j), nd_dst))
       else
          bi => layout_get_box_intersector(la_src,   box_nodalize(get_box(bxa_dst,j), nd_dst))
       end if

       do ii = 1, size(bi)
          i = bi(ii)%i

          if ( remote(la_dst,j) .and. remote(la_src,i) ) cycle

          if ( anynodal ) then
             fbox = get_box(lasrctmp,i)
          else
             fbox = get_box(la_src,  i)
          end if
          isect  = bi(ii)%bx
          loflux = lwb(fbox)

          if ( la_dst%lap%pmask(dir) .or. (loflux(dir) /= lo_dom(dir) .and. loflux(dir) /= hi_dom(dir)) ) then

             if ( local(la_dst,j) .and. local(la_src,i) ) then
                if ( li_r > size(flasc%flux%l_con%cpy) ) then
                   allocate(n_cpy(size(flasc%flux%l_con%cpy) * 2))
                   n_cpy(1:li_r-1) = flasc%flux%l_con%cpy(1:li_r-1)
                   deallocate(flasc%flux%l_con%cpy)
                   flasc%flux%l_con%cpy => n_cpy
                end if
                lcnt_r                         = lcnt_r + 1
                flasc%flux%l_con%cpy(li_r)%nd  = j
                flasc%flux%l_con%cpy(li_r)%ns  = i
                flasc%flux%l_con%cpy(li_r)%sbx = isect
                flasc%flux%l_con%cpy(li_r)%dbx = isect
                li_r                           = li_r + 1
             else if ( local(la_src,i) ) then
                if ( i_s > size(flasc%flux%r_con%snd) ) then
                   allocate(n_snd(size(flasc%flux%r_con%snd) * 2))
                   n_snd(1:i_s-1) = flasc%flux%r_con%snd(1:i_s-1)
                   deallocate(flasc%flux%r_con%snd)
                   flasc%flux%r_con%snd => n_snd
                   allocate(n_snd(size(flasc%mask%r_con%snd) * 2))
                   n_snd(1:i_s-1) = flasc%mask%r_con%snd(1:i_s-1)
                   deallocate(flasc%mask%r_con%snd)
                   flasc%mask%r_con%snd => n_snd
                end if
                cnt_s                         = cnt_s + 1
                parr(la_dst%lap%prc(j), 2)    = parr(la_dst%lap%prc(j), 2) + 1
                pvol(la_dst%lap%prc(j), 2)    = pvol(la_dst%lap%prc(j), 2) + volume(isect)
                flasc%flux%r_con%snd(i_s)%nd  = j
                flasc%flux%r_con%snd(i_s)%ns  = i
                flasc%flux%r_con%snd(i_s)%sbx = isect
                flasc%flux%r_con%snd(i_s)%dbx = isect
                flasc%flux%r_con%snd(i_s)%pr  = get_proc(la_dst,j)
                isect%lo(1:dm)                = isect%lo(1:dm) * ir(1:dm)
                isect%hi(1:dm)                = isect%hi(1:dm) * ir(1:dm)
                mpvol(la_dst%lap%prc(j), 2)   = mpvol(la_dst%lap%prc(j), 2) + volume(isect)
                flasc%mask%r_con%snd(i_s)%nd  = j
                flasc%mask%r_con%snd(i_s)%ns  = i
                flasc%mask%r_con%snd(i_s)%sbx = isect
                flasc%mask%r_con%snd(i_s)%dbx = isect
                flasc%mask%r_con%snd(i_s)%pr  = get_proc(la_dst,j)
                i_s                           = i_s + 1
             else
                if ( i_r > size(flasc%flux%r_con%rcv) ) then
                   allocate(n_rcv(size(flasc%flux%r_con%rcv) * 2))
                   n_rcv(1:i_r-1) = flasc%flux%r_con%rcv(1:i_r-1)
                   deallocate(flasc%flux%r_con%rcv)
                   flasc%flux%r_con%rcv => n_rcv
                   allocate(n_rcv(size(flasc%mask%r_con%rcv) * 2))
                   n_rcv(1:i_r-1) = flasc%mask%r_con%rcv(1:i_r-1)
                   deallocate(flasc%mask%r_con%rcv)
                   flasc%mask%r_con%rcv => n_rcv
                   allocate(pfbxs(size(flasc%fbxs) * 2))
                   pfbxs(1:i_r-1) = flasc%fbxs(1:i_r-1)
                   deallocate(flasc%fbxs)
                   flasc%fbxs => pfbxs
                end if
                cnt_r                         = cnt_r + 1
                parr(la_src%lap%prc(i), 1)    = parr(la_src%lap%prc(i), 1) + 1
                pvol(la_src%lap%prc(i), 1)    = pvol(la_src%lap%prc(i), 1) + volume(isect)
                flasc%flux%r_con%rcv(i_r)%nd  = j
                flasc%flux%r_con%rcv(i_r)%ns  = i
                flasc%flux%r_con%rcv(i_r)%sbx = isect
                flasc%flux%r_con%rcv(i_r)%dbx = isect
                flasc%flux%r_con%rcv(i_r)%pr  = get_proc(la_src,i)
                sh                            = 1
                sh(1:dm)                      = extent(isect)
                flasc%flux%r_con%rcv(i_r)%sh  = sh
                flasc%fbxs(i_r)               = fbox
                isect%lo(1:dm)                = isect%lo(1:dm) * ir(1:dm)
                isect%hi(1:dm)                = isect%hi(1:dm) * ir(1:dm)
                mpvol(la_src%lap%prc(i), 1)   = mpvol(la_src%lap%prc(i), 1) + volume(isect)
                flasc%mask%r_con%rcv(i_r)%nd  = j
                flasc%mask%r_con%rcv(i_r)%ns  = i
                flasc%mask%r_con%rcv(i_r)%sbx = isect
                flasc%mask%r_con%rcv(i_r)%dbx = isect
                flasc%mask%r_con%rcv(i_r)%pr  = get_proc(la_src,i)
                sh                            = 1
                sh(1:dm)                      = extent(isect)
                flasc%mask%r_con%rcv(i_r)%sh  = sh
                i_r                           = i_r + 1
             end if
          end if
       end do
       deallocate(bi)
    end do

    if ( anynodal ) then
       call destroy(lasrctmp)
    else
       call clear_box_hash_bin(la_src%lap)
    end if

    flasc%flux%l_con%ncpy = lcnt_r
    flasc%flux%r_con%nsnd = cnt_s
    flasc%flux%r_con%nrcv = cnt_r
    flasc%mask%r_con%nsnd = cnt_s
    flasc%mask%r_con%nrcv = cnt_r
    !
    ! Trim off unused space.
    !
    allocate(n_cpy(lcnt_r))
    n_cpy(1:lcnt_r) = flasc%flux%l_con%cpy(1:lcnt_r)
    deallocate(flasc%flux%l_con%cpy)
    flasc%flux%l_con%cpy => n_cpy

    allocate(n_snd(cnt_s))
    n_snd(1:cnt_s) = flasc%flux%r_con%snd(1:cnt_s)
    deallocate(flasc%flux%r_con%snd)
    flasc%flux%r_con%snd => n_snd
    allocate(n_snd(cnt_s))
    n_snd(1:cnt_s) = flasc%mask%r_con%snd(1:cnt_s)
    deallocate(flasc%mask%r_con%snd)
    flasc%mask%r_con%snd => n_snd

    allocate(n_rcv(cnt_r))
    n_rcv(1:cnt_r) = flasc%flux%r_con%rcv(1:cnt_r)
    deallocate(flasc%flux%r_con%rcv)
    flasc%flux%r_con%rcv => n_rcv
    allocate(n_rcv(cnt_r))
    n_rcv(1:cnt_r) = flasc%mask%r_con%rcv(1:cnt_r)
    deallocate(flasc%mask%r_con%rcv)
    flasc%mask%r_con%rcv => n_rcv
    allocate(pfbxs(cnt_r))
    pfbxs(1:cnt_r) = flasc%fbxs(1:cnt_r)
    deallocate(flasc%fbxs)
    flasc%fbxs => pfbxs

    do i = 0, np-1
       ppvol(i,1) = sum(pvol(0:i-1,1))
       ppvol(i,2) = sum(pvol(0:i-1,2))
    end do

    do i_r = 1, cnt_r
       i = flasc%flux%r_con%rcv(i_r)%pr
       flasc%flux%r_con%rcv(i_r)%pv = ppvol(i,1)
       pv = volume(flasc%flux%r_con%rcv(i_r)%dbx)
       ppvol(i,1) = ppvol(i,1) + pv
    end do

    do i_s = 1, cnt_s
       i = flasc%flux%r_con%snd(i_s)%pr
       flasc%flux%r_con%snd(i_s)%pv = ppvol(i,2)
       pv = volume(flasc%flux%r_con%snd(i_s)%dbx)
       ppvol(i,2) = ppvol(i,2) + pv
    end do

    flasc%flux%r_con%nrp  = count(parr(:,1) /= 0 )
    flasc%flux%r_con%nsp  = count(parr(:,2) /= 0 )
    flasc%flux%r_con%rvol = sum(pvol(:,1))
    flasc%flux%r_con%svol = sum(pvol(:,2))
    allocate(flasc%flux%r_con%rtr(flasc%flux%r_con%nrp))
    allocate(flasc%flux%r_con%str(flasc%flux%r_con%nsp))
    pi_r = 1; pi_s = 1; rpv = 0; spv = 0
    do i = 0, size(pvol,dim=1)-1
       if ( pvol(i,1) /= 0 ) then
          flasc%flux%r_con%rtr(pi_r)%sz = pvol(i,1)
          flasc%flux%r_con%rtr(pi_r)%pr = i
          flasc%flux%r_con%rtr(pi_r)%pv = rpv
          rpv  = rpv + pvol(i,1)
          pi_r = pi_r + 1
       end if
       if ( pvol(i,2) /= 0 ) then
          flasc%flux%r_con%str(pi_s)%sz = pvol(i,2)
          flasc%flux%r_con%str(pi_s)%pr = i
          flasc%flux%r_con%str(pi_s)%pv = spv
          spv  = spv + pvol(i,2)
          pi_s = pi_s + 1
       end if
    end do

    do i = 0, np-1
       ppvol(i,1) = sum(mpvol(0:i-1,1))
       ppvol(i,2) = sum(mpvol(0:i-1,2))
    end do

    do i_r = 1, cnt_r
       i = flasc%mask%r_con%rcv(i_r)%pr
       flasc%mask%r_con%rcv(i_r)%pv = ppvol(i,1)
       pv = volume(flasc%mask%r_con%rcv(i_r)%dbx)
       ppvol(i,1) = ppvol(i,1) + pv
    end do

    do i_s = 1, cnt_s
       i = flasc%mask%r_con%snd(i_s)%pr
       flasc%mask%r_con%snd(i_s)%pv = ppvol(i,2)
       pv = volume(flasc%mask%r_con%snd(i_s)%dbx)
       ppvol(i,2) = ppvol(i,2) + pv
    end do

    flasc%mask%r_con%nrp  = count(parr(:,1) /= 0 )
    flasc%mask%r_con%nsp  = count(parr(:,2) /= 0 )
    flasc%mask%r_con%rvol = sum(mpvol(:,1))
    flasc%mask%r_con%svol = sum(mpvol(:,2))
    allocate(flasc%mask%r_con%rtr(flasc%mask%r_con%nrp))
    allocate(flasc%mask%r_con%str(flasc%mask%r_con%nsp))
    pi_r = 1; pi_s = 1; rpv = 0; spv = 0
    do i = 0, size(mpvol,dim=1)-1
       if ( mpvol(i,1) /= 0 ) then
          flasc%mask%r_con%rtr(pi_r)%sz = mpvol(i,1)
          flasc%mask%r_con%rtr(pi_r)%pr = i
          flasc%mask%r_con%rtr(pi_r)%pv = rpv
          rpv  = rpv + mpvol(i,1)
          pi_r = pi_r + 1
       end if
       if ( mpvol(i,2) /= 0 ) then
          flasc%mask%r_con%str(pi_s)%sz = mpvol(i,2)
          flasc%mask%r_con%str(pi_s)%pr = i
          flasc%mask%r_con%str(pi_s)%pv = spv
          spv  = spv + mpvol(i,2)
          pi_s = pi_s + 1
       end if
    end do
    !
    ! Recall that fluxassoc contains two copyassocs.
    !
    call mem_stats_alloc(cpx_ms, flasc%flux%r_con%nsnd + flasc%flux%r_con%nrcv)
    call mem_stats_alloc(cpx_ms, flasc%flux%r_con%nsnd + flasc%flux%r_con%nrcv)
    call mem_stats_alloc(flx_ms, size(flasc%fbxs))

    call destroy(bpt)

  end subroutine fluxassoc_build

  subroutine copyassoc_destroy(cpasc)
    use bl_error_module
    type(copyassoc), intent(inout) :: cpasc
    if ( .not. built_q(cpasc) ) call bl_error("copyassoc_destroy(): not built")
    call mem_stats_dealloc(cpx_ms, cpasc%r_con%nsnd + cpasc%r_con%nrcv)
    call destroy(cpasc%ba_src)
    call destroy(cpasc%ba_dst)
    deallocate(cpasc%nd_dst)
    deallocate(cpasc%nd_src)
    deallocate(cpasc%l_con%cpy)
    deallocate(cpasc%r_con%snd)
    deallocate(cpasc%r_con%rcv)
    deallocate(cpasc%r_con%str)
    deallocate(cpasc%r_con%rtr)
    deallocate(cpasc%prc_src)
    deallocate(cpasc%prc_dst)
    cpasc%dim = 0
  end subroutine copyassoc_destroy

  subroutine fluxassoc_destroy(flasc)
    use bl_error_module
    type(fluxassoc), intent(inout) :: flasc
    if ( .not. built_q(flasc) ) call bl_error("fluxassoc_destroy(): not built")
    call mem_stats_dealloc(flx_ms, size(flasc%fbxs))
    call mem_stats_dealloc(cpx_ms, flasc%flux%r_con%nsnd + flasc%flux%r_con%nrcv)
    call mem_stats_dealloc(cpx_ms, flasc%flux%r_con%nsnd + flasc%flux%r_con%nrcv)

    deallocate(flasc%flux%l_con%cpy)
    deallocate(flasc%flux%r_con%snd)
    deallocate(flasc%flux%r_con%rcv)
    deallocate(flasc%mask%r_con%snd)
    deallocate(flasc%mask%r_con%rcv)
    deallocate(flasc%flux%r_con%rtr)
    deallocate(flasc%flux%r_con%str)
    deallocate(flasc%mask%r_con%rtr)
    deallocate(flasc%mask%r_con%str)

    deallocate(flasc%fbxs)
    deallocate(flasc%nd_dst)
    deallocate(flasc%nd_src)
    flasc%dim = 0
  end subroutine fluxassoc_destroy

  pure function layout_boxarray_hash(ba) result(r)
    type(boxarray), intent(in) :: ba
    integer                    :: r, i
    r = 0
    do i = 1, nboxes(ba)
       r = r + lwb(get_box(ba,i),1) + upb(get_box(ba,i),1)
    end do
  end function layout_boxarray_hash

  function copyassoc_check(cpasc, la_dst, la_src, nd_dst, nd_src) result(r)
    type(copyassoc), intent(in) :: cpasc
    type(layout),    intent(in) :: la_src, la_dst
    logical,         intent(in) :: nd_dst(:), nd_src(:)
    logical                     :: r
    integer                     :: hash
    r = all(cpasc%nd_dst .eqv. nd_dst) .and. all(cpasc%nd_src .eqv. nd_src)
    if (.not. r) return
    hash = layout_boxarray_hash(get_boxarray(la_src)) + layout_boxarray_hash(get_boxarray(la_dst))
    r = (hash .eq. cpasc%hash)
    if (.not. r) return
    r = boxarray_same_q(cpasc%ba_src, get_boxarray(la_src))
    if (.not. r) return
    r = boxarray_same_q(cpasc%ba_dst, get_boxarray(la_dst))
    if (.not. r) return
    r =  all(cpasc%prc_src == la_src%lap%prc)
    if (.not. r) return
    r =  all(cpasc%prc_dst == la_dst%lap%prc)
  end function copyassoc_check

  pure function fluxassoc_check(flasc, la_dst, la_src, nd_dst, nd_src, side, crse_domain, ir) result(r)
    logical                     :: r
    type(fluxassoc), intent(in) :: flasc
    type(layout),    intent(in) :: la_src, la_dst
    logical,         intent(in) :: nd_dst(:), nd_src(:)
    integer,         intent(in) :: side
    type(box),       intent(in) :: crse_domain
    integer,         intent(in) :: ir(:)
    r = associated(flasc%lap_dst, la_dst%lap) .and. associated(flasc%lap_src, la_src%lap)
    if (.not. r) return
    r = all(flasc%nd_dst .eqv. nd_dst) .and. all(flasc%nd_src .eqv. nd_src)
    if (.not. r) return
    r = (flasc%side .eq. side) .and. equal(flasc%crse_domain, crse_domain)
    if (.not. r) return
    r = all(flasc%ir(1:flasc%dim) .eq. ir(1:flasc%dim))
  end function fluxassoc_check

  function layout_copyassoc(la_dst, la_src, nd_dst, nd_src) result(r)

    use bl_error_module

    type(copyassoc)                :: r
    type(layout),    intent(inout) :: la_dst
    type(layout),    intent(in)    :: la_src
    logical,         intent(in)    :: nd_dst(:), nd_src(:)
    type(copyassoc), pointer       :: cp, ncp
    integer                        :: i, cnt, which
    !
    ! Do we have one stored?
    !
    i  =  1
    cp => the_copyassoc_head
    do while ( associated(cp) )
       if ( copyassoc_check(cp, la_dst, la_src, nd_dst, nd_src) ) then
          cp%reused = cp%reused + 1
          r = cp
          return
       end if
       cp => cp%next
       i  =  i + 1
    end do
    !
    ! Gotta build one.
    !
    allocate(cp)
    call copyassoc_build(cp, la_dst, la_src, nd_dst, nd_src)
    the_copyassoc_cnt = the_copyassoc_cnt + 1
    cp%next => the_copyassoc_head
    the_copyassoc_head => cp
    r = cp

    if ( the_copyassoc_cnt .gt. the_copyassoc_max ) then
       !
       ! We want to remove a least used copyassoc, but we never remove the
       ! first one.  That's the one we've (potentially) just built and will be
       ! passing back for use.  We want to get the latest entry in the list
       ! having the lowest reuse count.
       !
       i     =  1
       cnt   = Huge(1)
       which = -1
       cp => the_copyassoc_head
       do while ( associated(cp) )
          if ( i .gt. 1 .and. cp%reused .le. cnt ) then
             cnt   = cp%reused
             which = i
          end if
          cp => cp%next
          i  =  i + 1
       end do

       call bl_assert(which .le. the_copyassoc_cnt, 'which not in range')
       !
       ! Get rid of "which".
       !
       if ( which .gt. 1 ) then
          i =  1
          cp => the_copyassoc_head
          do while ( associated(cp) )
             if ( i .eq. (which-1) ) then
                ncp     => cp%next
                cp%next => ncp%next
                call copyassoc_destroy(ncp)
                deallocate(ncp)
                the_copyassoc_cnt = the_copyassoc_cnt - 1
                exit
             end if
             cp => cp%next
             i  =  i + 1
          end do
       end if
    end if

  end function layout_copyassoc

  function layout_fluxassoc(la_dst, la_src, nd_dst, nd_src, side, crse_domain, ir) result(r)
    type(fluxassoc)                :: r
    type(layout),    intent(in)    :: la_dst
    type(layout),    intent(inout) :: la_src
    logical,         intent(in)    :: nd_dst(:), nd_src(:)
    type(fluxassoc), pointer       :: fl
    integer,         intent(in)    :: side
    type(box),       intent(in)    :: crse_domain
    integer,         intent(in)    :: ir(:)
    !
    ! Do we have one stored?
    !
    fl => the_fluxassoc_head
    do while ( associated(fl) )
       if ( fluxassoc_check(fl, la_dst, la_src, nd_dst, nd_src, side, crse_domain, ir) ) then
          r = fl
          return
       end if
       fl => fl%next
    end do
    !
    ! Gotta build one.
    !
    allocate(fl)
    call fluxassoc_build(fl, la_dst, la_src, nd_dst, nd_src, side, crse_domain, ir)
    fl%next => the_fluxassoc_head
    the_fluxassoc_head => fl
    r = fl
  end function layout_fluxassoc

  pure function copyassoc_built_q(cpasc) result(r)
    logical :: r
    type(copyassoc), intent(in) :: cpasc
    r = cpasc%dim /= 0
  end function copyassoc_built_q

  pure function fluxassoc_built_q(flasc) result(r)
    logical :: r
    type(fluxassoc), intent(in) :: flasc
    r = flasc%dim /= 0
  end function fluxassoc_built_q

  subroutine clear_box_hash_bin(lap)
    type(layout_rep), intent(inout) :: lap

    integer :: i,j,k

    if ( associated(lap%bins) ) then

       do k = lbound(lap%bins,3), ubound(lap%bins,3)
          do j = lbound(lap%bins,2), ubound(lap%bins,2)
             do i = lbound(lap%bins,1), ubound(lap%bins,1)
                if ( associated(lap%bins(i,j,k)%iv) ) then
                   deallocate(lap%bins(i,j,k)%iv)
                end if
             end do
          end do
       end do
       deallocate(lap%bins)

       lap%bins => Null()
    end if

  end subroutine clear_box_hash_bin

  subroutine init_box_hash_bin(la, crsn)
    use bl_prof_module
    use bl_error_module
    type(layout), intent(inout) :: la
    integer, intent(in), optional :: crsn
    type(boxarray) :: ba
    integer, dimension(MAX_SPACEDIM) :: ext, vsz
    integer :: dm, n, i, j, k, cnt, full
    type(box) :: bx, cbx
    integer :: lcrsn
    integer :: sz
    type(box_hash_bin), pointer :: bins(:,:,:)
    integer, pointer :: ipv(:)
    type(bl_prof_timer), save :: bpt
    call build(bpt, "i_bx_hash")

    dm = la%lap%dim
    ba = get_boxarray(la)
    vsz = 0; vsz(1:dm) = -Huge(1)
    do n = 1, nboxes(ba)
       vsz(1:dm) = max(vsz(1:dm),extent(get_box(ba,n)))
    end do
    if ( present(crsn) ) then
       lcrsn = crsn
    else
       lcrsn = maxval(vsz)
    end if
    la%lap%crsn = lcrsn
    bx  = boxarray_bbox(ba)
    cbx = coarsen(bx, lcrsn)
    la%lap%plo = 0; la%lap%plo(1:dm) = lwb(cbx)
    la%lap%phi = 0; la%lap%phi(1:dm) = upb(cbx)
    la%lap%vshft = int_coarsen(vsz, lcrsn+1)

    allocate(la%lap%bins(la%lap%plo(1):la%lap%phi(1),la%lap%plo(2):la%lap%phi(2),la%lap%plo(3):la%lap%phi(3)))

    bins => la%lap%bins

    do n = 1, nboxes(ba)
       ext = 0; ext(1:dm) = int_coarsen(lwb(get_box(ba,n)), lcrsn)
       if ( .not. contains(cbx, ext(1:dm)) ) then
          call bl_error("init_box_hash_bin(): not contained!")
       end if
       if ( .not. associated(bins(ext(1),ext(2),ext(3))%iv) ) then
          allocate(ipv(1))
          ipv(1) = n
       else
          sz = size(bins(ext(1),ext(2),ext(3))%iv)
          allocate(ipv(sz+1))
          ipv(1:sz) = bins(ext(1),ext(2),ext(3))%iv(1:sz)
          ipv(sz+1) = n
          deallocate(bins(ext(1),ext(2),ext(3))%iv)
       end if
       bins(ext(1),ext(2),ext(3))%iv => ipv
    end do

    if ( .false. .and. parallel_IOProcessor() ) then
       cnt = size(bins,dim=1)*size(bins,dim=2)*size(bins,dim=3); full = 0
       do k = la%lap%plo(3),la%lap%phi(3)
          do j = la%lap%plo(2),la%lap%phi(2)
             do i = la%lap%plo(1),la%lap%phi(1)
                if ( associated(bins(i,j,k)%iv) ) full = full + 1
             end do
          end do
       end do
!       if ( cnt > 0 ) then
!          write(6,'(A I6 A I3)') '*** init_box_hash_bin: bins: ', cnt, ' %full: ', INT(REAL(full)/cnt*100)
!       endif
    end if

    call destroy(bpt)
  end subroutine init_box_hash_bin

  function layout_get_box_intersector(la, bx) result(bi)
    use bl_error_module
    type(box_intersector), pointer :: bi(:)
    type(layout), intent(inout) :: la
    type(box), intent(in) :: bx
    type(box_hash_bin), pointer :: bins(:,:,:)
    integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
    integer :: dm
    type(box) :: bx1
    integer :: i, j, k, n
    type(boxarray) :: ba
    integer, parameter :: ChunkSize = 50
    integer :: cnt
    type(box_intersector), pointer :: tbi(:)

    !$OMP CRITICAL(initboxhashbin)
    if (.not. associated(la%lap%bins)) call init_box_hash_bin(la)
    !$OMP END CRITICAL(initboxhashbin)

    allocate(bi(ChunkSize))

    dm = la%lap%dim
    ba = get_boxarray(la)
    bins => la%lap%bins
    bx1 = coarsen(bx, la%lap%crsn)
    lo = 0; lo(1:dm) = lwb(bx1)
    hi = 0; hi(1:dm) = upb(bx1)
    cnt = 0
    select case ( dm ) 
    case (3)
       do k = max(lo(3)-la%lap%vshft(3)-1,la%lap%plo(3)), min(hi(3)+la%lap%vshft(3), la%lap%phi(3))
          do j = max(lo(2)-la%lap%vshft(2)-1,la%lap%plo(2)), min(hi(2)+la%lap%vshft(2), la%lap%phi(2))
             do i = max(lo(1)-la%lap%vshft(1)-1,la%lap%plo(1)), min(hi(1)+la%lap%vshft(1), la%lap%phi(1))

                if ( associated(bins(i,j,k)%iv) ) then
                   do n = 1, size(bins(i,j,k)%iv)
                      bx1 = intersection(bx, get_box(ba,bins(i,j,k)%iv(n)))
                      if ( empty(bx1) ) cycle
                      cnt = cnt + 1

                      if (cnt > size(bi)) then
                         allocate(tbi(size(bi) * 2))
                         tbi(1:cnt-1) = bi(1:cnt-1)
                         deallocate(bi)
                         bi => tbi
                      end if

                      bi(cnt)%i  = bins(i,j,k)%iv(n)
                      bi(cnt)%bx = bx1
                   end do
                end if

             end do
          end do
       end do
    case (2)
       k = 0
       do j = max(lo(2)-la%lap%vshft(2)-1,la%lap%plo(2)), min(hi(2)+la%lap%vshft(2), la%lap%phi(2))
          do i = max(lo(1)-la%lap%vshft(1)-1,la%lap%plo(1)), min(hi(1)+la%lap%vshft(1), la%lap%phi(1))

             if ( associated(bins(i,j,k)%iv) ) then
                do n = 1, size(bins(i,j,k)%iv)
                   bx1 = intersection(bx, get_box(ba,bins(i,j,k)%iv(n)))
                   if ( empty(bx1) ) cycle
                   cnt = cnt + 1

                   if (cnt > size(bi)) then
                      allocate(tbi(size(bi) * 2))
                      tbi(1:cnt-1) = bi(1:cnt-1)
                      deallocate(bi)
                      bi => tbi
                   end if

                   bi(cnt)%i  = bins(i,j,k)%iv(n)
                   bi(cnt)%bx = bx1
                end do
             end if

          end do
       end do
    case (1)
       k = 0
       j = 0
       do i = max(lo(1)-la%lap%vshft(1)-1,la%lap%plo(1)), min(hi(1)+la%lap%vshft(1), la%lap%phi(1))

          if ( associated(bins(i,j,k)%iv) ) then
             do n = 1, size(bins(i,j,k)%iv)
                bx1 = intersection(bx, get_box(ba,bins(i,j,k)%iv(n)))
                if ( empty(bx1) ) cycle
                cnt = cnt + 1

                if (cnt > size(bi)) then
                   allocate(tbi(size(bi) * 2))
                   tbi(1:cnt-1) = bi(1:cnt-1)
                   deallocate(bi)
                   bi => tbi
                end if

                bi(cnt)%i  = bins(i,j,k)%iv(n)
                bi(cnt)%bx = bx1
             end do
          end if

       end do
    end select

    allocate(tbi(cnt))
    tbi(1:cnt) = bi(1:cnt)
    deallocate(bi)
    bi => tbi

  end function layout_get_box_intersector
  !
  ! What's in "bx" excluding what's in get_boxarray(la)
  !
  subroutine layout_boxarray_diff(ba, bx, la)
    type(boxarray), intent(out  )   :: ba
    type(layout),   intent(inout)   :: la
    type(box),      intent(in   )   :: bx
    type(list_box)                  :: bl1, bl
    integer                         :: i
    type(box_intersector), pointer  :: bi(:)
    call build(bl1)
    bi => layout_get_box_intersector(la, bx)
    do i = 1, size(bi)
       call push_back(bl1, bi(i)%bx)
    end do
    deallocate(bi)
    bl = boxlist_boxlist_diff(bx, bl1)
    call boxarray_build_l(ba, bl)
    call destroy(bl)
    call destroy(bl1)
  end subroutine layout_boxarray_diff

  ! build special copyassoc for bndry_reg copy to other
  subroutine copyassoc_build_br_to_other(cpasc, la_dst, la_src, nd_dst, nd_src)
    type(copyassoc), intent(inout) :: cpasc
    type(layout), intent(in) :: la_dst, la_src
    logical, intent(in) :: nd_dst(:), nd_src(:)

    type(box) :: bx, bx_src
    integer :: i0, i, j, pv, rpv, spv, pi_r, pi_s, pcnt_r, pcnt_s
    integer :: lcnt_r, li_r, cnt_r, cnt_s, i_r, i_s, sh(MAX_SPACEDIM+1)
    type(local_copy_desc), pointer :: n_cpy(:) => Null()
    type(comm_dsc), pointer :: n_snd(:) => Null(), n_rcv(:) => Null()
    integer, allocatable :: pvol(:,:), ppvol(:,:), parr(:,:)
    integer, parameter :: ChunkSize = 500

    cpasc%dim = get_dim(la_dst)
    
    allocate(cpasc%nd_dst(cpasc%dim));  cpasc%nd_dst = nd_dst
    allocate(cpasc%nd_src(cpasc%dim));  cpasc%nd_src = nd_src

    allocate(cpasc%l_con%cpy(ChunkSize))
    allocate(cpasc%r_con%snd(ChunkSize))
    allocate(cpasc%r_con%rcv(ChunkSize))

    allocate(cpasc%prc_src(la_src%lap%nboxes))
    allocate(cpasc%prc_dst(la_dst%lap%nboxes))
    cpasc%prc_src = la_src%lap%prc
    cpasc%prc_dst = la_dst%lap%prc

    allocate(parr(0:parallel_nprocs()-1,2))
    allocate(pvol(0:parallel_nprocs()-1,2))
    allocate(ppvol(0:parallel_nprocs()-1,2))

    parr = 0; pvol = 0; lcnt_r = 0; cnt_r = 0; cnt_s = 0; li_r = 1; i_r = 1; i_s = 1

    i0 = 1
    do j = 1, nboxes(la_src) 
       bx_src = get_box(la_src, j)

       do i = i0, nboxes(la_dst)
          bx = get_box(la_dst, i)

          if ( .not. box_contains(bx_src, bx)) then
             i0 = i
             exit
          else
             if ( remote(la_dst,i) .and. remote(la_src,j) ) then
                cycle
             else if ( local(la_dst, i) .and. local(la_src, j) ) then
                if ( li_r > size(cpasc%l_con%cpy) ) then
                   allocate(n_cpy(size(cpasc%l_con%cpy) * 2))
                   n_cpy(1:li_r-1) = cpasc%l_con%cpy(1:li_r-1)
                   deallocate(cpasc%l_con%cpy)
                   cpasc%l_con%cpy => n_cpy
                end if
                lcnt_r                    = lcnt_r + 1
                cpasc%l_con%cpy(li_r)%nd  = i
                cpasc%l_con%cpy(li_r)%ns  = j
                cpasc%l_con%cpy(li_r)%sbx = bx
                cpasc%l_con%cpy(li_r)%dbx = bx
                li_r                      = li_r + 1
             else if ( local(la_src, j) ) then ! dst is remote
                if ( i_s > size(cpasc%r_con%snd) ) then
                   allocate(n_snd(size(cpasc%r_con%snd) * 2))
                   n_snd(1:i_s-1) = cpasc%r_con%snd(1:i_s-1)
                   deallocate(cpasc%r_con%snd)
                   cpasc%r_con%snd => n_snd
                end if
                cnt_s                      = cnt_s + 1
                parr(la_dst%lap%prc(i), 2) = parr(la_dst%lap%prc(i), 2) + 1
                pvol(la_dst%lap%prc(i), 2) = pvol(la_dst%lap%prc(i), 2) + volume(bx)
                cpasc%r_con%snd(i_s)%nd    = i
                cpasc%r_con%snd(i_s)%ns    = j
                cpasc%r_con%snd(i_s)%sbx   = bx
                cpasc%r_con%snd(i_s)%dbx   = bx
                cpasc%r_con%snd(i_s)%pr    = get_proc(la_dst,i)
                i_s                        = i_s + 1
             else ! src is remote, dst is local
                if ( i_r > size(cpasc%r_con%rcv) ) then
                   allocate(n_rcv(size(cpasc%r_con%rcv) * 2))
                   n_rcv(1:i_r-1) = cpasc%r_con%rcv(1:i_r-1)
                   deallocate(cpasc%r_con%rcv)
                   cpasc%r_con%rcv => n_rcv
                end if
                cnt_r                      = cnt_r + 1
                parr(la_src%lap%prc(j), 1) = parr(la_src%lap%prc(j), 1) + 1
                pvol(la_src%lap%prc(j), 1) = pvol(la_src%lap%prc(j), 1) + volume(bx)
                cpasc%r_con%rcv(i_r)%nd    = i
                cpasc%r_con%rcv(i_r)%ns    = j
                cpasc%r_con%rcv(i_r)%sbx   = bx
                cpasc%r_con%rcv(i_r)%dbx   = bx
                cpasc%r_con%rcv(i_r)%pr    = get_proc(la_src,j)
                sh                         = 1
                sh(1:cpasc%dim)            = extent(bx)
                cpasc%r_con%rcv(i_r)%sh    = sh
                i_r                        = i_r + 1
             end if
          end if
       end do
    end do

    cpasc%l_con%ncpy = lcnt_r
    cpasc%r_con%nsnd = cnt_s
    cpasc%r_con%nrcv = cnt_r
    !
    ! Trim off unused space.
    !
    allocate(n_cpy(lcnt_r))
    n_cpy(1:lcnt_r) = cpasc%l_con%cpy(1:lcnt_r)
    deallocate(cpasc%l_con%cpy)
    cpasc%l_con%cpy => n_cpy

    allocate(n_snd(cnt_s))
    n_snd(1:cnt_s) = cpasc%r_con%snd(1:cnt_s)
    deallocate(cpasc%r_con%snd)
    cpasc%r_con%snd => n_snd

    allocate(n_rcv(cnt_r))
    n_rcv(1:cnt_r) = cpasc%r_con%rcv(1:cnt_r)
    deallocate(cpasc%r_con%rcv)
    cpasc%r_con%rcv => n_rcv
    !
    ! This region packs the src/recv boxes into processor order
    !
    do i = 0, parallel_nprocs()-1
       ppvol(i,1) = sum(pvol(0:i-1,1))
       ppvol(i,2) = sum(pvol(0:i-1,2))
    end do
    !
    ! Pack Receives maintaining original ordering
    !
    do i_r = 1, cnt_r
       i = cpasc%r_con%rcv(i_r)%pr
       cpasc%r_con%rcv(i_r)%pv = ppvol(i,1)
       pv = volume(cpasc%r_con%rcv(i_r)%dbx)
       ppvol(i,1) = ppvol(i,1) + pv
    end do
    !
    ! Pack Sends maintaining original ordering
    !
    do i_s = 1, cnt_s
       i = cpasc%r_con%snd(i_s)%pr
       cpasc%r_con%snd(i_s)%pv = ppvol(i,2)
       pv = volume(cpasc%r_con%snd(i_s)%dbx)
       ppvol(i,2) = ppvol(i,2) + pv
    end do
    !
    ! Now compute the volume of data the each processor expects
    !
    pcnt_r = count(parr(:,1) /= 0 )
    pcnt_s = count(parr(:,2) /= 0 )
    cpasc%r_con%nrp  = pcnt_r
    cpasc%r_con%nsp  = pcnt_s
    cpasc%r_con%rvol = sum(pvol(:,1))
    cpasc%r_con%svol = sum(pvol(:,2))
    allocate(cpasc%r_con%str(pcnt_s))
    allocate(cpasc%r_con%rtr(pcnt_r))
    pi_r = 1; pi_s = 1; rpv  = 0; spv  = 0
    do i = 0, size(pvol,dim=1)-1
       if ( pvol(i,1) /= 0 ) then
          cpasc%r_con%rtr(pi_r)%sz = pvol(i,1)
          cpasc%r_con%rtr(pi_r)%pr = i
          cpasc%r_con%rtr(pi_r)%pv = rpv
          rpv  = rpv + pvol(i,1)
          pi_r = pi_r + 1
       end if
       if ( pvol(i,2) /= 0 ) then
          cpasc%r_con%str(pi_s)%sz = pvol(i,2)
          cpasc%r_con%str(pi_s)%pr = i
          cpasc%r_con%str(pi_s)%pv = spv
          spv  = spv + pvol(i,2)
          pi_s = pi_s + 1
       end if
    end do

    call mem_stats_alloc(cpx_ms, cpasc%r_con%nsnd + cpasc%r_con%nrcv)

    deallocate(pvol, ppvol, parr)
  end subroutine copyassoc_build_br_to_other

end module layout_module
