module layout_module

  use parallel
  use boxarray_module
  use knapsack_module
  use bl_mem_stat_module
  use bl_prof_module

  implicit none

  integer, private, parameter :: LA_UNDF = 0
  integer, private, parameter :: LA_BASE = 1
  integer, private, parameter :: LA_CRSN = 2
  integer, private, parameter :: LA_PCHD = 3
  integer, private, parameter :: LA_PPRT = 4
  integer, private, parameter :: LA_DERV = 5

  integer, parameter :: LA_KNAPSACK   = 101
  integer, parameter :: LA_ROUNDROBIN = 102
  integer, parameter :: LA_LOCAL      = 103
  integer, parameter :: LA_EXPLICIT   = 104

  integer, private :: def_mapping = LA_KNAPSACK

  type comm_dsc
     integer   :: nd = 0        ! dst box number
     integer   :: ns = 0        ! src box number
     type(box) :: dbx           ! dst sub-box
     type(box) :: sbx           ! src sub-box
     integer   :: pv = 0        ! number of points in buf prior to this
     integer   :: av = 0        ! number of points in buf including this one
     integer   :: sh(MAX_SPACEDIM+1) = 0 ! shape for data from rcvbuf
     integer   :: s1(1) = 0     ! shape for data from rcvbuf
     integer   :: pr            ! Processors number of src or dest
  end type comm_dsc

  type trns_dsc
     integer :: sz = 0          ! Size of chunk 
     integer :: pv = -1         ! Number points in buf prior to this
     integer :: pr = MPI_ANY_SOURCE ! src or destination processor
  end type trns_dsc

  type remote_conn
     integer :: svol = 0        ! numpts in snd volume
     integer :: rvol = 0        ! numpts in rcv volume
     integer :: nsnd = 0        ! Number of snd chunks, redundant
     integer :: nrcv = 0        ! Number of rcv chunks, redundant
     type(comm_dsc), pointer :: snd(:) => Null()
     type(comm_dsc), pointer :: rcv(:) => Null()
     integer :: nrp  = 0        ! Number of processes receiving from, redundant
     integer :: nsp  = 0        ! Number of processes sending to, redundant
     type(trns_dsc), pointer :: str(:) => Null()
     type(trns_dsc), pointer :: rtr(:) => Null()
  end type remote_conn

  type local_copy_desc
     integer   :: ns = 0         ! Source box in layout
     integer   :: nd = 0         ! Destination box in layout
     type(box) :: sbx            ! Sub-box for this copy
     type(box) :: dbx            ! Sub-box for this copy
  end type local_copy_desc

  type local_conn
     integer :: ncpy            ! Number of cpy chunks, redundant
     type(local_copy_desc), pointer :: cpy(:) => Null()
  end type local_conn

  type boxassoc
     integer :: dim    = 0                  ! spatial dimension 1, 2, or 3
     integer :: nboxes = 0                  ! number of boxes
     integer :: grwth  = 0                  ! growth factor
     logical, pointer :: nodal(:) => Null() ! nodal flag
     logical :: cross = .false.             ! cross/full stencil?
     type(local_conn)  :: l_con
     type(remote_conn) :: r_con
     type(boxassoc), pointer :: next => Null()
  end type boxassoc

  type syncassoc
     integer :: dim    = 0                  ! spatial dimension 1, 2, or 3
     integer :: nboxes = 0                  ! number of boxes
     integer :: grwth  = 0                  ! growth factor
     logical :: lall   = .false.            ! use valid region or everything
     logical, pointer :: nodal(:) => Null() ! nodal flag
     type(local_conn)  :: l_con
     type(remote_conn) :: r_con
     type(syncassoc), pointer :: next => Null()
  end type syncassoc

  type copyassoc
     integer           :: dim = 0             ! spatial dimension 1, 2, or 3
     logical, pointer  :: nd_dst(:) => Null() ! dst nodal flag
     logical, pointer  :: nd_src(:) => Null() ! src nodal flag
     type(local_conn)  :: l_con
     type(remote_conn) :: r_con
     type(copyassoc),  pointer :: next    => Null()
     type(layout_rep), pointer :: lap_dst => Null()
     type(layout_rep), pointer :: lap_src => Null()
  end type copyassoc

  type box_intersector
     integer :: i
     type(box) :: bx
  end type box_intersector

  type box_hash_bin
     integer, pointer :: iv(:) => Null()
  end type box_hash_bin

  !
  ! Global list of copyassoc's used by multifab copy routines.
  !
  type(copyassoc), pointer, save :: the_copyassoc_head => Null()

  type layout
     integer :: la_type = LA_UNDF
     type(layout_rep), pointer :: lap => Null()
  end type layout

  !! Defines the box distribution and box connectivity
  !! of a boxarray
  type layout_rep
     integer :: dim = 0         ! spatial dimension 1, 2, or 3
     integer :: id  = 0
     integer :: nboxes = 0
     type(box) :: pd            ! Problem Domain 
     logical, pointer  :: pmask(:) => Null() ! periodic mask
     integer, pointer, dimension(:) :: prc => Null()
     type(boxarray) :: bxa
     type(boxassoc), pointer :: bxasc => Null()
     type(syncassoc), pointer :: snasc => Null()
     type(coarsened_layout), pointer :: crse_la => Null()
     type(pn_layout), pointer :: pn_children => Null()
     type(derived_layout), pointer :: dlay => Null()
     ! Box Hashing
     integer :: crsn = -1
     integer :: plo(MAX_SPACEDIM) = 0
     integer :: phi(MAX_SPACEDIM) = 0
     integer :: vshft(MAX_SPACEDIM) = 0
     type(box_hash_bin), pointer :: bins(:,:,:) => Null()
  end type layout_rep

  !! A layout that is derived by coarsening an existing layout,
  !! The processor distribution and the number of boxes will be
  !! the same as for the parent layout.  The intent is to be used
  !! in multigrid solvers that keep coarsened grids on the same
  !! processor as their parent in the hierarchy.
  type coarsened_layout
     integer :: dim = 0
     integer, pointer :: crse(:) => Null()
     type(layout) :: la
     type(coarsened_layout), pointer :: next => Null()
  end type coarsened_layout

  type pn_layout
     integer :: dim = 0
     integer, pointer :: refr(:) => Null()
     type(layout) :: la
     type(pn_layout), pointer :: next => Null()
  end type pn_layout

  type derived_layout
     integer :: dim = 0
     type(layout) :: la
     type(derived_layout), pointer :: next => Null()
  end type derived_layout

  integer, private :: g_layout_next_id = 0;

  interface built_q
     module procedure layout_built_q
     module procedure boxassoc_built_q
     module procedure syncassoc_built_q
     module procedure copyassoc_built_q
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

  interface nboxes
     module procedure layout_nboxes
  end interface

  interface remote
     module procedure layout_remote
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

  private layout_next_id
  private layout_rep_build
  private layout_rep_destroy

  type(mem_stats), private, save :: bxa_ms
  type(mem_stats), private, save :: la_ms

contains

  subroutine layout_set_mapping(mapping)
    integer, intent(in) :: mapping
    def_mapping = mapping
  end subroutine layout_set_mapping
  function layout_get_mapping() result(r)
    integer :: r
    r = def_mapping
  end function layout_get_mapping

  subroutine boxassoc_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    bxa_ms = ms
  end subroutine boxassoc_set_mem_stats
  subroutine layout_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    la_ms = ms
  end subroutine layout_set_mem_stats

  function boxassoc_mem_stats() result(r)
    type(mem_stats) :: r
    r = bxa_ms
  end function boxassoc_mem_stats
  function layout_mem_stats() result(r)
    type(mem_stats) :: r
    r = la_ms
  end function layout_mem_stats

  function layout_next_id() result(r)
    integer :: r
    g_layout_next_id = g_layout_next_id + 1
    r = g_layout_next_id
  end function layout_next_id

  function layout_not_equal(a,b) result(r)
    type(layout), intent(in) :: a, b
    logical :: r
    r = .not. associated(a%lap, b%lap)
  end function layout_not_equal

  function layout_equal(a,b) result(r)
    type(layout), intent(in) :: a, b
    logical :: r
    r = associated(a%lap, b%lap)
  end function layout_equal

  function layout_built_q(la) result(r)
    logical :: r
    type(layout), intent(in) :: la
    r = associated(la%lap)
  end function layout_built_q

  function layout_dim(la) result(r)
    integer :: r
    type(layout), intent(in) :: la
    r = la%lap%dim
  end function layout_dim

  function layout_nboxes(la) result(r)
    integer :: r
    type(layout), intent(in) :: la
    r = la%lap%nboxes
  end function layout_nboxes

  function layout_get_pd(la) result(r)
    type(box) :: r
    type(layout), intent(in) :: la
    r = la%lap%pd
  end function layout_get_pd

  function layout_boxarray(la) result(r)
    type(layout), intent(in) :: la
    type(boxarray) :: r
    r = la%lap%bxa
  end function layout_boxarray

  function layout_get_pmask(la) result(r)
    type(layout), intent(in) :: la
    logical :: r(la%lap%dim)
    r = la%lap%pmask
  end function layout_get_pmask

  subroutine layout_rep_build(lap, ba, pd, pmask, mapping, explicit_mapping)
    type(layout_rep), intent(out) :: lap
    type(boxarray), intent(in) :: ba
    type(box), intent(in) :: pd
    logical, intent(in) :: pmask(:)
    integer, intent(in), optional :: mapping
    integer, intent(in), optional :: explicit_mapping(:)
    integer :: i, j
    integer :: lmapping

    lmapping = def_mapping; if ( present(mapping) ) lmapping = mapping
    if ( present(explicit_mapping) ) then
       if ( present(mapping) ) then
          if ( mapping /= LA_EXPLICIT ) then
             call bl_error("LAYOUT_REP_BUILD:explicit_mapping doesn't match mapping")
          end if
       end if
       lmapping = LA_EXPLICIT
    end if
    call boxarray_build_copy(lap%bxa, ba)
    lap%dim    = lap%bxa%dim
    lap%nboxes = lap%bxa%nboxes
    lap%id     = layout_next_id()
    lap%pd     = pd
    allocate(lap%pmask(lap%dim))
    lap%pmask = pmask

    allocate(lap%prc(lap%nboxes))
    select case (lmapping)
    case (LA_EXPLICIT)
       if ( .not. present(explicit_mapping) ) then
          call bl_error("LAYOUT_REP_BUILD: mapping explicit but no explicit_mapping")
       end if
       if ( size(lap%prc) /= size(explicit_mapping) ) then
          call bl_error("LAYOUT_REP_BUILD: incommesurate explicit mapping")
       end if
       lap%prc = explicit_mapping
    case (LA_LOCAL)
       lap%prc = parallel_myproc()
    case (LA_ROUNDROBIN)
       call layout_roundrobin(lap%prc, ba%bxs)
    case (LA_KNAPSACK)
       call layout_knapsack(lap%prc, ba%bxs)
    case default
       call bl_error("LAYOUT_REP_BUILD: unknown mapping:", lmapping)
    end select

  end subroutine layout_rep_build

  recursive subroutine layout_rep_destroy(lap, la_type)
    type(layout_rep), pointer :: lap
    integer, intent(in) :: la_type
    type(coarsened_layout), pointer :: clp, oclp
    type(pn_layout), pointer :: pnp, opnp
    type(derived_layout), pointer :: dla, odla
    type(boxassoc),  pointer :: bxa, obxa
    type(syncassoc), pointer :: snxa, osnxa
    type(copyassoc), pointer :: cpa, ncpa, pcpa
    integer :: i, j, k
    if ( la_type /= LA_CRSN ) then
       deallocate(lap%prc)
    end if
    if ( la_type /= LA_DERV ) then
       call destroy(lap%bxa)
    end if
    if ( la_type == LA_BASE .or. la_type == LA_PCHD ) then
       deallocate(lap%pmask)
    end if
    clp => lap%crse_la
    do while ( associated(clp) )
       oclp => clp%next
       deallocate(clp%crse)
       call layout_rep_destroy(clp%la%lap, LA_CRSN)
!      deallocate(clp%la%lap)
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
    dla => lap%dlay
    do while ( associated(dla) )
       odla => dla%next
       call layout_rep_destroy(dla%la%lap, LA_DERV)
       deallocate(dla)
       dla  => odla
    end do
    !
    ! Get rid of boxassocs
    !
    bxa => lap%bxasc
    do while ( associated(bxa) )
       obxa => bxa%next
       call boxassoc_destroy(bxa)
       deallocate(bxa)
       bxa => obxa
    end do
    !
    ! Get rid of syncassocs
    !
    snxa => lap%snasc
    do while ( associated(snxa) )
       osnxa => snxa%next
       call syncassoc_destroy(snxa)
       deallocate(snxa)
       snxa => osnxa
    end do
    !
    ! remove any boxarray hash
    !
    if ( associated(lap%bins) ) then
       do k = lbound(lap%bins,3), ubound(lap%bins,3)
          do j = lbound(lap%bins,2), ubound(lap%bins,2)
             do i = lbound(lap%bins,1), ubound(lap%bins,1)
                deallocate(lap%bins(i,j,k)%iv)
             end do
          end do
       end do
       deallocate(lap%bins)
    end if
    !
    ! Remove all copyassoc's associated with this layout_rep.
    !
    cpa  => the_copyassoc_head
    pcpa => Null()
    do while ( associated(cpa) )
       ncpa => cpa%next
       if ( associated(lap, cpa%lap_src) .or. associated(lap, cpa%lap_dst) ) then
          if ( associated(cpa, the_copyassoc_head) ) then
             the_copyassoc_head => cpa%next
          else
             pcpa%next => ncpa
          end if
          call copyassoc_destroy(cpa)
          deallocate(cpa)
       else
          if ( .not. associated(pcpa) ) then
             pcpa => the_copyassoc_head
          else
             pcpa => pcpa%next
          end if
       end if
       cpa => ncpa
    end do
    deallocate(lap)
  end subroutine layout_rep_destroy

  subroutine layout_build_ba(la, ba, pd, pmask, mapping, explicit_mapping)
    type(layout), intent(out) :: la
    type(boxarray), intent(in) :: ba
    type(box), intent(in), optional :: pd
    logical, intent(in), optional :: pmask(:)
    integer, intent(in), optional :: mapping
    integer, intent(in), optional :: explicit_mapping(:)
    type(box) :: lpd
    logical :: lpmask(ba%dim)
    lpmask = .false.; if ( present(pmask) ) lpmask = pmask
    if ( present(pd) ) then
       lpd = pd
    else
       lpd = boxarray_bbox(ba)
    end if
    allocate(la%lap)
    la%la_type = LA_BASE
    call layout_rep_build(la%lap, ba, lpd, lpmask, mapping, explicit_mapping)
  end subroutine layout_build_ba

  subroutine layout_destroy(la)
    type(layout), intent(inout) :: la
    if ( la%la_type /= LA_BASE ) call bl_error("LAYOUT_DESTROY: confused")
    call layout_rep_destroy(la%lap, LA_BASE)
!   deallocate(la%lap)
  end subroutine layout_destroy

  subroutine layout_build_pn(lapn, la, ba, rr, mapping, explicit_mapping)
    type(layout), intent(out)   :: lapn
    type(layout), intent(inout) :: la
    type(boxarray), intent(in) :: ba
    integer, intent(in) :: rr(:)
    integer, intent(in), optional :: mapping
    integer, intent(in), optional :: explicit_mapping(:)
    type(pn_layout), pointer :: pla
    type(box) :: rpd

    if ( size(rr) /= la%lap%dim ) then
       call bl_error("LAYOUT_BUILD_PN: incommensurate refinement ratio")
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

  subroutine layout_build_derived(lad, la, prc, root)
    type(layout), intent(out) :: lad
    type(layout), intent(inout) :: la
    integer, intent(in), optional :: prc(:)
    integer, intent(in), optional :: root
    integer :: l_root
    type(derived_layout), pointer :: dla
    integer :: l_prc(la%lap%nboxes)
    integer :: i, j

    l_root = -1
    if ( present(prc) ) then
       if ( present(root) ) &
            call bl_error("LAYOUT_BUILD_DERIVED: not both root and prc")
       if ( size(prc) /= la%lap%nboxes ) &
            call bl_error("LAYOUT_BUILD_DERIVED: incommensurate prc")
    else if ( present(root) ) then
       l_root = root
    else if ( .not. present(prc) ) then
       l_root = parallel_IOProcessorNode()
    end if
    if ( l_root == -1 ) then
       !! handle prc case
       l_prc = prc
    else
       !! handle root case
       l_prc = l_root
    end if

    dla => la%lap%dlay
    do while ( associated(dla) )
       if ( all(dla%la%lap%prc == l_prc) ) then
          lad = dla%la
          return
       end if
       dla => dla%next
    end do

    ! not found
    allocate(dla)
    dla%dim = la%lap%dim

    dla%la%la_type = LA_DERV
    allocate(dla%la%lap)
    dla%la%lap%dim = la%lap%dim
    dla%la%lap%nboxes = la%lap%nboxes
    dla%la%lap%id  = layout_next_id()
    dla%la%lap%pd = la%lap%pd
    dla%la%lap%pmask => la%lap%pmask

    allocate(dla%la%lap%prc(size(l_prc)))
    dla%la%lap%prc = l_prc
    dla%la%lap%bxa = la%lap%bxa

    ! install the new derived into the layout
    dla%next => la%lap%dlay
    la%lap%dlay => dla
    lad = dla%la

  end subroutine layout_build_derived

  subroutine layout_build_coarse(lac, la, cr)
    type(layout), intent(out)   :: lac
    type(layout), intent(inout) :: la
    integer, intent(in) :: cr(:)
    type(coarsened_layout), pointer :: cla

    if ( size(cr) /= la%lap%dim ) then
       call bl_error("LAYOUT_BUILD_COARSE: incommensurate cr")
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
    cla%la%lap%pd = coarsen(la%lap%pd, cla%crse)
    cla%la%lap%pmask => la%lap%pmask

    cla%la%lap%prc => la%lap%prc

    call boxarray_build_v(cla%la%lap%bxa, la%lap%bxa%bxs)

    call boxarray_coarsen(cla%la%lap%bxa, cla%crse)

    ! install the new coarsened layout into the layout
    cla%next => la%lap%crse_la
    la%lap%crse_la => cla
    lac = cla%la
  end subroutine layout_build_coarse

  function layout_remote(la, i) result(r)
    type(layout), intent(in) :: la
    integer, intent(in) :: i
    logical :: r
    r = la%lap%prc(i) /= parallel_myproc()
  end function layout_remote

  function layout_local(la, i) result(r)
    type(layout), intent(in) :: la
    integer, intent(in) :: i
    logical :: r
    r = la%lap%prc(i) == parallel_myproc()
  end function layout_local

  function layout_get_box(la, i) result(r)
    type(layout), intent(in) :: la
    integer, intent(in) :: i
    type(box) :: r
    r = get_box(la%lap%bxa, i)
  end function layout_get_box

  function layout_get_proc(la, i) result(r)
    type(layout), intent(in) :: la
    integer, intent(in) :: i
    integer :: r
    r = la%lap%prc(i)
  end function layout_get_proc

  function layout_get_proc_v(la) result(r)
    type(layout), intent(in) :: la
    integer :: r(size(la%lap%prc))
    r = la%lap%prc
  end function layout_get_proc_v

  subroutine layout_roundrobin(prc, bxs)
    integer, intent(out), dimension(:) :: prc
    type(box), intent(in), dimension(:) :: bxs
    integer :: i
    prc = mod((/(i,i=0,size(bxs)-1)/),parallel_nprocs())
  end subroutine layout_roundrobin

  subroutine layout_knapsack(prc, bxs)
    use knapsack_module
    integer, intent(out), dimension(:) :: prc
    type(box), intent(in), dimension(:) :: bxs
    integer :: ibxs(size(bxs))
    ibxs = volume(bxs)
    call knapsack_i(prc, ibxs, parallel_nprocs())
  end subroutine layout_knapsack

  function layout_efficiency(la, np) result(r)
    type(layout), intent(in) :: la
    real(kind=dp_t) :: r
    integer, intent(in), optional :: np
    real(kind=dp_t) :: weights(la%lap%nboxes)
    real(kind=dp_t) :: p_max_weight
    integer :: i, lnp
    lnp = parallel_nprocs(); if ( present(np) ) lnp = np
    weights = box_dvolume(la%lap%bxa%bxs)
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
       write(unit=un, fmt='(" NBOXES  = ",i2)') la%lap%nboxes
       call unit_skip(unit, skip)
       write(unit=un, fmt='(" PD      = ",i2)', advance = 'no')
       call print(la%lap%pd, unit = unit)
       do i = 1, nboxes(la)
          call unit_skip(unit, skip = unit_get_skip(skip) + 1)
          write(unit=un, fmt = '(I0,": ")', advance = 'no') i
          call print(get_box(la,i), unit = unit, advance = 'no')
          write(unit=un, fmt = '(" ",I0)') get_proc(la,i)
       end do
       write(unit=un, fmt = '(" *)]")')
    end if
  end subroutine layout_print

  function boxassoc_check(bxa, ng, nodal, cross) result(r)
    logical :: r
    type(boxassoc), intent(in) :: bxa
    integer, intent(in) :: ng
    logical, intent(in), optional :: nodal(:)
    logical, intent(in), optional :: cross
    if ( present(nodal) ) then
       r = bxa%grwth == ng .and. all(bxa%nodal .eqv. nodal)
    else
       r = bxa%grwth == ng .and. all(bxa%nodal .eqv. .false.)
    end if
    if ( present(cross) ) then
       r = r .and. ( bxa%cross .eqv. cross )
    else
       r = r .and. ( bxa%cross .eqv. .false. )
    end if
  end function boxassoc_check

  function syncassoc_check(snxa, ng, nodal, lall) result(r)
    logical                     :: r
    type(syncassoc), intent(in) :: snxa
    integer, intent(in)         :: ng
    logical, intent(in)         :: nodal(:)
    logical, intent(in)         :: lall
    r = snxa%grwth == ng .and. all(snxa%nodal .eqv. nodal) .and. (snxa%lall .eqv. lall)
  end function syncassoc_check

  function layout_boxassoc(la, ng, nodal, cross) result(r)
    type(boxassoc) :: r
    type(layout) , intent(inout) :: la
    integer, intent(in) :: ng
    logical, intent(in), optional :: nodal(:)
    logical, intent(in), optional :: cross
    type(boxassoc), pointer :: bp

    bp => la%lap%bxasc
    do while ( associated(bp) )
       if ( boxassoc_check(bp, ng, nodal, cross) ) then
          r = bp
          return
       end if
       bp => bp%next
    end do
    !
    ! Didn't find; so have to go looking for it.
    !
    allocate (bp)
    call boxassoc_build(bp, la%lap, ng, nodal, cross)
    bp%next => la%lap%bxasc
    la%lap%bxasc => bp
    r = bp
  end function layout_boxassoc

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
    ! Didn't find; so have to go looking for it.
    !
    allocate (sp)
    call syncassoc_build(sp, la%lap, ng, nodal, lall)
    sp%next => la%lap%snasc
    la%lap%snasc => sp
    r = sp
  end function layout_syncassoc

  function boxassoc_built_q(bxasc) result(r)
    logical :: r
    type(boxassoc), intent(in) :: bxasc
    r = bxasc%dim /= 0
  end function boxassoc_built_q

  function syncassoc_built_q(snasc) result(r)
    logical :: r
    type(syncassoc), intent(in) :: snasc
    r = snasc%dim /= 0
  end function syncassoc_built_q

  subroutine boxarray_bndry_periodic(bxai, dmn, b, nodal, pmask, ng, shfts, cross)
    type(boxarray), intent(out)          :: bxai
    type(box),      intent(in)           :: dmn, b
    logical,        intent(in)           :: nodal(:), pmask(:)
    integer,        intent(in)           :: ng
    integer,        intent(out)          :: shfts(:,:)
    logical,        intent(in), optional :: cross

    integer               :: i, j, cnt, emptylo(3), emptyhi(3)
    type(box)             :: bxs(3**b%dim), gbx, emptybx
    type(box),allocatable :: bv(:)
    integer               :: shft(3**b%dim,b%dim)
    type(boxarray)        :: tba, cba, dba
    logical               :: lcross
    type(bl_prof_timer), save :: bpt

    call build(bpt, "ba_bndry_periodic")

    lcross = .false.; if ( present(cross) ) lcross = cross

    call boxarray_box_boundary_n(tba, box_nodalize(b,nodal), ng)

    if ( lcross ) then
       call boxarray_box_corners(cba, box_nodalize(b,nodal), ng)
       call boxarray_diff(tba, cba)
    end if

    shfts = 0

    call box_periodic_shift(dmn, b, nodal, pmask, ng, shft, cnt, bxs)

    if ( cnt > 0 ) then
       if ( lcross ) then
          emptylo = 0
          emptyhi = -1
          emptybx = make_box(emptylo(1:b%dim), emptyhi(1:b%dim))
          gbx     = grow(box_nodalize(b,nodal), ng)
          do i = 1, cnt
             call boxarray_build_bx(dba, intersection(gbx, shift(bxs(i), -shft(i,:))))
             call boxarray_diff(dba, cba)
             if ( nboxes(dba) > 1) call bl_error("BOXARRAY_BNDRY_PERIODIC: nboxes(dba) > 1")
             if ( empty(dba) ) then
                bxs(i) = emptybx
             else
                bxs(i) = shift(dba%bxs(1), shft(i,:))
             endif
             call destroy(dba)
          end do
       end if
       allocate(bv(tba%nboxes+cnt))
       bv(1:tba%nboxes) = tba%bxs(1:tba%nboxes)
       bv(tba%nboxes+1:tba%nboxes+cnt) = bxs(1:cnt)
       shfts(tba%nboxes+1:tba%nboxes+cnt,:) = shft(1:cnt,:)
       call destroy(tba)
       call boxarray_build_v(tba, bv, sort = .false.)
    end if

    if ( lcross ) call destroy(cba)

    bxai = tba

    call destroy(bpt)

  end subroutine boxarray_bndry_periodic

  subroutine boxassoc_build(bxasc, lap, ng, nodal, cross)

    integer,          intent(in)           :: ng
    logical,          intent(in), optional :: nodal(:)
    type(layout_rep), intent(in), target   :: lap
    type(boxassoc),   intent(inout)        :: bxasc
    logical,          intent(in), optional :: cross

    integer                        :: i, j, ii, pv, rpv, spv, pi_r, pi_s, pcnt_r, pcnt_s
    integer                        :: shft(2*3**lap%dim,lap%dim), sh(MAX_SPACEDIM+1)
    type(box)                      :: abx, bx
    type(boxarray)                 :: bxa, bxai
    type(layout)                   :: la
    integer                        :: lcnt_r, li_r, cnt_r, cnt_s, i_r, i_s
    integer                        :: lcnt_r_max, cnt_r_max, cnt_s_max
    integer, parameter             :: chunksize = 100
    integer, allocatable           :: pvol(:,:), ppvol(:,:), parr(:,:)
    type(local_copy_desc), pointer :: n_cpy(:) => Null()
    type(comm_dsc), pointer        :: n_snd(:) => Null(), n_rcv(:) => Null()
    logical                        :: first
    type(bl_prof_timer), save      :: bpt
    type(box), allocatable :: abxx(:)
    logical, allocatable   :: is_empty(:)

    if ( built_q(bxasc) ) call bl_error("BOXASSOC_BUILD: already built")

    call build(bpt, "boxassoc_build")

    la%lap       => lap
    bxa          =  get_boxarray(la)
    bxasc%dim    =  bxa%dim
    bxasc%grwth  =  ng
    bxasc%nboxes =  bxa%nboxes

    allocate(bxasc%nodal(bxasc%dim))
    allocate(parr(0:parallel_nprocs()-1,2))
    allocate(pvol(0:parallel_nprocs()-1,2))
    allocate(ppvol(0:parallel_nprocs()-1,2))
    allocate(bxasc%l_con%cpy(chunksize))
    allocate(bxasc%r_con%snd(chunksize))
    allocate(bxasc%r_con%rcv(chunksize))

    allocate(is_empty(26),abxx(26))

    bxasc%nodal = .false.; if ( present(nodal) ) bxasc%nodal = nodal

    parr = 0; pvol = 0; lcnt_r = 0; cnt_r = 0; cnt_s = 0; li_r = 1; i_r = 1; i_s = 1
    !
    ! Consider all copies I <- J.
    !
    do i = 1, bxa%nboxes
       first = .true.
       do j = 1, bxa%nboxes
          if ( remote(la,i) .and. remote(la,j) ) cycle
          if ( first ) then
             call boxarray_bndry_periodic(bxai, lap%pd, bxa%bxs(i), bxasc%nodal, lap%pmask, ng, shft, cross)
             first = .false.
             if ( size(is_empty) < bxai%nboxes) then
                deallocate(is_empty,abxx)
                allocate(is_empty(bxai%nboxes),abxx(bxai%nboxes))
             end if
          end if
          bx = box_nodalize(get_box(bxa, j), nodal)
          call box_intersection_and_empty(abxx, is_empty, bx, bxai%bxs)
          do ii = 1, bxai%nboxes; if ( is_empty(ii) ) cycle
             abx = abxx(ii)
             if ( local(la,i) .and. local(la, j) ) then
                if ( li_r > size(bxasc%l_con%cpy) ) then
                   allocate(n_cpy(size(bxasc%l_con%cpy) + chunksize))
                   n_cpy(1:li_r-1) = bxasc%l_con%cpy(1:li_r-1)
                   deallocate(bxasc%l_con%cpy)
                   bxasc%l_con%cpy => n_cpy
                end if
                lcnt_r                    = lcnt_r + 1
                bxasc%l_con%cpy(li_r)%nd  = i
                bxasc%l_con%cpy(li_r)%ns  = j
                bxasc%l_con%cpy(li_r)%sbx = abx
                bxasc%l_con%cpy(li_r)%dbx = shift(abx,-shft(ii,:))
                li_r                      = li_r + 1
             else if ( local(la, j) ) then
                if ( i_s > size(bxasc%r_con%snd) ) then
                   allocate(n_snd(size(bxasc%r_con%snd) + chunksize))
                   n_snd(1:i_s-1) = bxasc%r_con%snd(1:i_s-1)
                   deallocate(bxasc%r_con%snd)
                   bxasc%r_con%snd => n_snd
                end if
                cnt_s                    = cnt_s + 1
                parr(lap%prc(i), 2)      = parr(lap%prc(i), 2) + 1
                pvol(lap%prc(i), 2)      = pvol(lap%prc(i), 2) + volume(abx)
                bxasc%r_con%snd(i_s)%nd  = i
                bxasc%r_con%snd(i_s)%ns  = j
                bxasc%r_con%snd(i_s)%sbx = abx
                bxasc%r_con%snd(i_s)%dbx = shift(abx,-shft(ii,:))
                bxasc%r_con%snd(i_s)%pr  = get_proc(la, i)
                bxasc%r_con%snd(i_s)%s1  = volume(abx)
                i_s                      = i_s + 1
             else if ( local(la, i) ) then
                if ( i_r > size(bxasc%r_con%rcv) ) then
                   allocate(n_rcv(size(bxasc%r_con%rcv) + chunksize))
                   n_rcv(1:i_r-1) = bxasc%r_con%rcv(1:i_r-1)
                   deallocate(bxasc%r_con%rcv)
                   bxasc%r_con%rcv => n_rcv
                end if
                cnt_r                    = cnt_r + 1
                parr(lap%prc(j), 1)      = parr(lap%prc(j), 1) + 1
                pvol(lap%prc(j), 1)      = pvol(lap%prc(j), 1) + volume(abx)
                bxasc%r_con%rcv(i_r)%nd  = i
                bxasc%r_con%rcv(i_r)%ns  = j
                bxasc%r_con%rcv(i_r)%sbx = abx
                bxasc%r_con%rcv(i_r)%dbx = shift(abx,-shft(ii,:))
                bxasc%r_con%rcv(i_r)%pr  = get_proc(la, j)
                sh                       = 1
                sh(1:bxasc%dim)          = extent(abx)
                bxasc%r_con%rcv(i_r)%sh  = sh
                i_r                      = i_r + 1
             end if
          end do
       end do
       if ( .not. first ) call destroy(bxai)
    end do

    if ( .false. ) then
       call parallel_reduce(lcnt_r_max, lcnt_r, MPI_MAX, proc = parallel_IOProcessorNode())
       call parallel_reduce(cnt_s_max,   cnt_s, MPI_MAX, proc = parallel_IOProcessorNode())
       call parallel_reduce(cnt_r_max,   cnt_r, MPI_MAX, proc = parallel_IOProcessorNode())
       if ( parallel_IOProcessor() ) then
          print*, '*** chunksize = ', chunksize
          print*, '*** max(lcnt_r) = ', lcnt_r_max
          print*, '*** max(cnt_s) = ', cnt_s_max
          print*, '*** max(cnt_r) = ', cnt_r_max
       end if
    end if

    bxasc%l_con%ncpy = lcnt_r
    bxasc%r_con%nsnd = cnt_s
    bxasc%r_con%nrcv = cnt_r

    allocate(n_cpy(lcnt_r))
    n_cpy(1:lcnt_r) = bxasc%l_con%cpy(1:lcnt_r)
    deallocate(bxasc%l_con%cpy)
    bxasc%l_con%cpy => n_cpy

    allocate(n_snd(cnt_s))
    n_snd(1:cnt_s)  = bxasc%r_con%snd(1:cnt_s)
    deallocate(bxasc%r_con%snd)
    bxasc%r_con%snd => n_snd

    allocate(n_rcv(cnt_r))
    n_rcv(1:cnt_r)  = bxasc%r_con%rcv(1:cnt_r)
    deallocate(bxasc%r_con%rcv)
    bxasc%r_con%rcv => n_rcv
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
       i = bxasc%r_con%rcv(i_r)%pr
       bxasc%r_con%rcv(i_r)%pv = ppvol(i,1)
       pv = volume(bxasc%r_con%rcv(i_r)%dbx)
       bxasc%r_con%rcv(i_r)%av = bxasc%r_con%rcv(i_r)%pv + pv
       ppvol(i,1) = ppvol(i,1) + pv
    end do
    !
    ! Pack Sends maintaining original ordering
    !
    do i_s = 1, cnt_s
       i = bxasc%r_con%snd(i_s)%pr
       bxasc%r_con%snd(i_s)%pv = ppvol(i,2)
       pv = volume(bxasc%r_con%snd(i_s)%dbx)
       bxasc%r_con%snd(i_s)%av = bxasc%r_con%snd(i_s)%pv + pv
       ppvol(i,2) = ppvol(i,2) + pv
    end do
    !
    ! Now compute the volume of data the each processor expects
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
    call mem_stats_alloc(bxa_ms)

    call destroy(bpt)

  end subroutine boxassoc_build

  subroutine internal_sync_unique_cover(la, ng, nodal, lall, filled)

    type(layout), intent(in)  :: la
    integer, intent(in)       :: ng
    logical, intent(in)       :: nodal(:)
    logical, intent(in)       :: lall
    type(local_conn), pointer :: filled(:)

    type(box)                      :: ibx, jbx, abx
    integer                        :: i, j, k, jj, cnt
    integer                        :: shft(2*3**(la%lap%dim),la%lap%dim)
    integer, parameter             :: chunksize = 10
    type(local_copy_desc)          :: lcd
    type(local_copy_desc), pointer :: n_cpy(:) => Null()
    type(list_box)                 :: lb1, lb2
    type(boxarray)                 :: bxa, ba1, ba2
    type(box), allocatable         :: bxs(:)

    bxa = get_boxarray(la)

    allocate(filled(bxa%nboxes))

    do i = 1, bxa%nboxes
       filled(i)%ncpy = 0
       allocate(filled(i)%cpy(chunksize))
    end do

    do j = 1, bxa%nboxes
       jbx = box_nodalize(bxa%bxs(j), nodal)
       if ( lall ) jbx = grow(jbx,ng)
       call box_internal_sync_shift(la%lap%pd, jbx, la%lap%pmask, nodal, shft, cnt)
       do i = j, bxa%nboxes
          ibx = box_nodalize(bxa%bxs(i), nodal)
          if ( lall ) ibx = grow(ibx,ng)
          do jj = 1, cnt
             !
             ! Do not overwrite ourselves.
             !
             if ( i == j .and. all(shft(jj,:) == 0) ) cycle
             abx = intersection(ibx, shift(jbx,shft(jj,:)))
             if ( empty(abx) ) cycle
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
                   allocate(n_cpy(size(filled(i)%cpy) + chunksize))
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
       end do
    end do
    !
    ! Test that we're a unique cover; i.e. no overlap.  Is there a better way to do this?
    !
    if ( .false. ) then
       do i = 1, bxa%nboxes
          if ( filled(i)%ncpy > 0 ) then
             allocate(bxs(filled(i)%ncpy))
             do j = 1, filled(i)%ncpy
                bxs(j) = filled(i)%cpy(j)%dbx
             end do
             call boxarray_add_clean_boxes(ba1, bxs, simplify = .false.)
             call boxarray_build_v(ba2, bxs, sort = .false.)
             if ( .not. boxarray_same_q(ba1, ba2) ) then
                print*, "*** NOT a unique covering !!!"
                call print(ba1, "ba1")
                call print(ba2, "ba2")
                stop
             end if
             deallocate(bxs)
             call destroy(ba1)
             call destroy(ba2)
          end if
       end do
    end if

  end subroutine internal_sync_unique_cover

  subroutine syncassoc_build(snasc, lap, ng, nodal, lall)

    integer,          intent(in)         :: ng
    logical,          intent(in)         :: nodal(:)
    type(layout_rep), intent(in), target :: lap
    type(syncassoc),  intent(inout)      :: snasc
    logical,          intent(in)         :: lall

    integer                        :: i, j, ii, jj, pv, rpv, spv, pi_r, pi_s, pcnt_r, pcnt_s
    type(box)                      :: dbx, sbx
    type(boxarray)                 :: bxa
    type(layout)                   :: la
    integer                        :: lcnt_r_max, cnt_r_max, cnt_s_max, cnt
    integer                        :: lcnt_r, li_r, cnt_r, cnt_s, i_r, i_s, sh(MAX_SPACEDIM+1)
    integer, parameter             :: chunksize = 100
    integer, allocatable           :: pvol(:,:), ppvol(:,:), parr(:,:)
    type(local_copy_desc), pointer :: n_cpy(:) => Null()
    type(comm_dsc), pointer        :: n_snd(:) => Null(), n_rcv(:) => Null()
    type(local_conn), pointer      :: filled(:)
    type(bl_prof_timer), save      :: bpt

    if ( built_q(snasc) ) call bl_error("SYNCASSOC_BUILD: already built")

    call build(bpt, "syncassoc_build")

    la%lap       => lap
    bxa          =  get_boxarray(la)
    snasc%dim    =  bxa%dim
    snasc%grwth  =  ng
    snasc%nboxes =  bxa%nboxes

    allocate(snasc%nodal(snasc%dim))
    allocate(parr(0:parallel_nprocs()-1,2))
    allocate(pvol(0:parallel_nprocs()-1,2))
    allocate(ppvol(0:parallel_nprocs()-1,2))
    allocate(snasc%l_con%cpy(chunksize))
    allocate(snasc%r_con%snd(chunksize))
    allocate(snasc%r_con%rcv(chunksize))

    snasc%lall  = lall
    snasc%nodal = nodal

    call internal_sync_unique_cover(la, snasc%grwth, snasc%nodal, snasc%lall, filled)

    parr = 0; pvol = 0; lcnt_r = 0; cnt_r = 0; cnt_s = 0; li_r = 1; i_r = 1; i_s = 1

    do jj = 1, bxa%nboxes
       if ( filled(jj)%ncpy > 0 ) then
          do ii = 1, filled(jj)%ncpy
             i   = filled(jj)%cpy(ii)%nd
             j   = filled(jj)%cpy(ii)%ns
             sbx = filled(jj)%cpy(ii)%sbx
             dbx = filled(jj)%cpy(ii)%dbx
             if ( local(la, i) .and. local(la, j) ) then
                if ( li_r > size(snasc%l_con%cpy) ) then
                   allocate(n_cpy(size(snasc%l_con%cpy) + chunksize))
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
                   allocate(n_snd(size(snasc%r_con%snd) + chunksize))
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
                snasc%r_con%snd(i_s)%s1  = volume(dbx)
                i_s                      = i_s + 1
             else if ( local(la, i) ) then  ! must recv
                if ( i_r > size(snasc%r_con%rcv) ) then
                   allocate(n_rcv(size(snasc%r_con%rcv) + chunksize))
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

    do i = 1, bxa%nboxes
       deallocate(filled(i)%cpy)
    end do
    deallocate(filled)

    snasc%l_con%ncpy = lcnt_r
    snasc%r_con%nsnd = cnt_s
    snasc%r_con%nrcv = cnt_r

    allocate(n_cpy(lcnt_r))
    n_cpy(1:lcnt_r) = snasc%l_con%cpy(1:lcnt_r)
    deallocate(snasc%l_con%cpy)
    snasc%l_con%cpy => n_cpy

    allocate(n_snd(cnt_s))
    n_snd(1:cnt_s)  = snasc%r_con%snd(1:cnt_s)
    deallocate(snasc%r_con%snd)
    snasc%r_con%snd => n_snd

    allocate(n_rcv(cnt_r))
    n_rcv(1:cnt_r)  = snasc%r_con%rcv(1:cnt_r)
    deallocate(snasc%r_con%rcv)
    snasc%r_con%rcv => n_rcv
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
       i = snasc%r_con%rcv(i_r)%pr
       snasc%r_con%rcv(i_r)%pv = ppvol(i,1)
       pv = volume(snasc%r_con%rcv(i_r)%dbx)
       snasc%r_con%rcv(i_r)%av = snasc%r_con%rcv(i_r)%pv + pv
       ppvol(i,1) = ppvol(i,1) + pv
    end do
    !
    ! Pack Sends maintaining original ordering
    !
    do i_s = 1, cnt_s
       i = snasc%r_con%snd(i_s)%pr
       snasc%r_con%snd(i_s)%pv = ppvol(i,2)
       pv = volume(snasc%r_con%snd(i_s)%dbx)
       snasc%r_con%snd(i_s)%av = snasc%r_con%snd(i_s)%pv + pv
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

    call destroy(bpt)

  end subroutine syncassoc_build

  subroutine boxassoc_destroy(bxasc)
    type(boxassoc), intent(inout) :: bxasc
    if ( .not. built_q(bxasc) ) call bl_error("BOXASSOC_DESTROY: not built")
    deallocate(bxasc%nodal)
    deallocate(bxasc%l_con%cpy)
    deallocate(bxasc%r_con%snd)
    deallocate(bxasc%r_con%rcv)
    deallocate(bxasc%r_con%str)
    deallocate(bxasc%r_con%rtr)
    call mem_stats_dealloc(bxa_ms)
  end subroutine boxassoc_destroy

  subroutine syncassoc_destroy(snasc)
    type(syncassoc), intent(inout) :: snasc
    if ( .not. built_q(snasc) ) call bl_error("SYNCASSOC_DESTROY: not built")
    deallocate(snasc%nodal)
    deallocate(snasc%l_con%cpy)
    deallocate(snasc%r_con%snd)
    deallocate(snasc%r_con%rcv)
    deallocate(snasc%r_con%str)
    deallocate(snasc%r_con%rtr)
  end subroutine syncassoc_destroy

  subroutine boxassoc_print(bxasc, str, unit, skip)
    use bl_IO_module
    type(boxassoc), intent(in) :: bxasc
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: skip
    integer :: un
    integer :: i, j, ii
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un,fmt='("BOXASSOC")', advance='no')
    if ( present(str) ) then
       write(unit=un, fmt='(": ", A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un,fmt='(" NODAL = ", 3L2)') bxasc%nodal
    do ii = 0, parallel_nprocs()-1
       if ( ii == parallel_myproc() ) then
          call unit_skip(un, skip)
          write(unit=un, fmt='(" PROCESSOR ", i4)') ii
          call unit_skip(un, skip)
          write(unit=un, fmt='(" L_CON")')
          do i = 1, bxasc%l_con%ncpy
             call unit_skip(un, skip)
             write(unit=un, fmt='(" ",i5,":(",i4,"<-",i4,"): ")', advance = 'no') &
                  i, &
                  bxasc%l_con%cpy(i)%nd, &
                  bxasc%l_con%cpy(i)%ns
             call print(bxasc%l_con%cpy(i)%dbx, unit=un, advance = 'no')
             write(unit=un, fmt='(" <-- ")', advance = 'no')
             call print(bxasc%l_con%cpy(i)%sbx, unit=un)
          end do
          call unit_skip(un, skip)
          write(unit=un, fmt='(" R_CON")')
          call unit_skip(un, skip)
          write(unit=un, fmt='(" S BUF: Volume: ", i10)') bxasc%r_con%svol
          do i = 1, bxasc%r_con%nsp
             call unit_skip(un, skip)
             write(unit=un, fmt='(" ", i4,":",i4,":",i4,":",i4)') i, &
                  bxasc%r_con%str(i)%pr, &
                  bxasc%r_con%str(i)%pv, &
                  bxasc%r_con%str(i)%sz
          end do
          call unit_skip(un, skip)
          write(unit=un, fmt='(" R BUF: Volume: ", i10)') bxasc%r_con%rvol
          do i = 1, bxasc%r_con%nrp
             call unit_skip(un, skip)
             write(unit=un, fmt='(i4,":",i4,":",i4,":",i4)') i, &
                  bxasc%r_con%rtr(i)%pr, &
                  bxasc%r_con%rtr(i)%pv, &
                  bxasc%r_con%rtr(i)%sz
          end do
          write(unit=un, fmt='(" SND")')
          do i = 1, bxasc%r_con%nsnd
             call unit_skip(un, skip)
             write(unit=un, fmt='(" ",i5,":(",i4,"<-",i4,"):",i4,":",i4,":",i4,":",i4,":",4(1x,i4),": ")', &
                  advance = 'no') &
                  i, &
                  bxasc%r_con%snd(i)%nd, &
                  bxasc%r_con%snd(i)%ns, &
                  bxasc%r_con%snd(i)%pr, &
                  bxasc%r_con%snd(i)%pv, &
                  bxasc%r_con%snd(i)%av, &
                  bxasc%r_con%snd(i)%s1, &
                  bxasc%r_con%snd(i)%sh
             call print(bxasc%r_con%snd(i)%sbx, unit=un)
          end do
          call unit_skip(un, skip)
          write(unit=un, fmt='(" RCV")')
          do i = 1, bxasc%r_con%nrcv
             call unit_skip(un, skip)
             write(unit=un, fmt='(" ",i5,":(",i4,"<-",i4,"):",i4,":",i4,":",i4,":",i4,":",4(1x,i4),": ")', &
                  advance = 'no') &
                  i, &
                  bxasc%r_con%rcv(i)%nd, &
                  bxasc%r_con%rcv(i)%ns, &
                  bxasc%r_con%rcv(i)%pr, &
                  bxasc%r_con%rcv(i)%pv, &
                  bxasc%r_con%rcv(i)%av, &
                  bxasc%r_con%rcv(i)%s1, &
                  bxasc%r_con%rcv(i)%sh
             call print(bxasc%r_con%rcv(i)%dbx, unit=un)
          end do
       end if
       call parallel_barrier()
    end do
  end subroutine boxassoc_print

  subroutine copyassoc_build(cpasc, la_dst, la_src, nd_dst, nd_src)

    type(copyassoc),  intent(inout) :: cpasc
    type(layout),     intent(in)    :: la_src, la_dst
    logical,          intent(in)    :: nd_dst(:), nd_src(:)

    integer                        :: i, j, pv, rpv, spv, pi_r, pi_s, pcnt_r, pcnt_s
    integer                        :: sh(MAX_SPACEDIM+1)
    type(box)                      :: bx
    type(boxarray)                 :: bxa_src, bxa_dst
    integer                        :: lcnt_r, li_r, cnt_r, cnt_s, i_r, i_s
    integer                        :: lcnt_r_max, cnt_r_max, cnt_s_max
    integer, allocatable           :: pvol(:,:), ppvol(:,:), parr(:,:)
    type(local_copy_desc), pointer :: n_cpy(:) => Null()
    type(comm_dsc), pointer        :: n_snd(:) => Null(), n_rcv(:) => Null()
    integer, parameter             :: chunksize = 100
    type(bl_prof_timer), save      :: bpt

    if ( built_q(cpasc) ) call bl_error("COPYASSOC_BUILD: already built")

    call build(bpt, "copyassoc_build")

    bxa_src       =  get_boxarray(la_src)
    bxa_dst       =  get_boxarray(la_dst)
    cpasc%dim     =  bxa_src%dim
    cpasc%lap_dst => la_dst%lap
    cpasc%lap_src => la_src%lap

    allocate(cpasc%nd_dst(la_dst%lap%dim))
    allocate(cpasc%nd_src(la_src%lap%dim))
    allocate(parr(0:parallel_nprocs()-1,2))
    allocate(pvol(0:parallel_nprocs()-1,2))
    allocate(ppvol(0:parallel_nprocs()-1,2))
    allocate(cpasc%l_con%cpy(chunksize))
    allocate(cpasc%r_con%snd(chunksize))
    allocate(cpasc%r_con%rcv(chunksize))

    cpasc%nd_dst = nd_dst
    cpasc%nd_src = nd_src

    parr = 0; pvol = 0; lcnt_r = 0; cnt_r = 0; cnt_s = 0; li_r = 1; i_r = 1; i_s = 1
    !
    ! Consider all copies I <- J.
    !
    do i = 1, bxa_dst%nboxes
       do j = 1, bxa_src%nboxes
          if ( remote(la_dst,i) .and. remote(la_src,j) ) cycle
          bx = intersection(box_nodalize(get_box(bxa_dst,i),nd_dst), box_nodalize(get_box(bxa_src,j),nd_src))
          if ( empty(bx) ) cycle
          if ( local(la_dst, i) .and. local(la_src, j) ) then
             if ( li_r > size(cpasc%l_con%cpy) ) then
                allocate(n_cpy(size(cpasc%l_con%cpy) + chunksize))
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
                allocate(n_snd(size(cpasc%r_con%snd) + chunksize))
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
             cpasc%r_con%snd(i_s)%s1    = volume(bx)
             i_s                        = i_s + 1
          else if ( local(la_dst, i) ) then
             if ( i_r > size(cpasc%r_con%rcv) ) then
                allocate(n_rcv(size(cpasc%r_con%rcv) + chunksize))
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
    end do

    if ( .false. ) then
       call parallel_reduce(lcnt_r_max, lcnt_r, MPI_MAX, proc = parallel_IOProcessorNode())
       call parallel_reduce(cnt_s_max,   cnt_s, MPI_MAX, proc = parallel_IOProcessorNode())
       call parallel_reduce(cnt_r_max,   cnt_r, MPI_MAX, proc = parallel_IOProcessorNode())
       if ( parallel_IOProcessor() ) then
          print*, '*** chunksize = ', chunksize
          print*, '*** max(lcnt_r) = ', lcnt_r_max
          print*, '*** max(cnt_s) = ', cnt_s_max
          print*, '*** max(cnt_r) = ', cnt_r_max
       end if
    end if

    cpasc%l_con%ncpy = lcnt_r
    cpasc%r_con%nsnd = cnt_s
    cpasc%r_con%nrcv = cnt_r

    allocate(n_cpy(lcnt_r))
    n_cpy(1:lcnt_r) = cpasc%l_con%cpy(1:lcnt_r)
    deallocate(cpasc%l_con%cpy)
    cpasc%l_con%cpy => n_cpy

    allocate(n_snd(cnt_s))
    n_snd(1:cnt_s)  = cpasc%r_con%snd(1:cnt_s)
    deallocate(cpasc%r_con%snd)
    cpasc%r_con%snd => n_snd

    allocate(n_rcv(cnt_r))
    n_rcv(1:cnt_r)  = cpasc%r_con%rcv(1:cnt_r)
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
       cpasc%r_con%rcv(i_r)%av = cpasc%r_con%rcv(i_r)%pv + pv
       ppvol(i,1) = ppvol(i,1) + pv
    end do
    !
    ! Pack Sends maintaining original ordering
    !
    do i_s = 1, cnt_s
       i = cpasc%r_con%snd(i_s)%pr
       cpasc%r_con%snd(i_s)%pv = ppvol(i,2)
       pv = volume(cpasc%r_con%snd(i_s)%dbx)
       cpasc%r_con%snd(i_s)%av = cpasc%r_con%snd(i_s)%pv + pv
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

    call destroy(bpt)

  end subroutine copyassoc_build

  subroutine copyassoc_destroy(cpasc)
    type(copyassoc), intent(inout) :: cpasc
    if ( .not. built_q(cpasc) ) call bl_error("COPYASSOC_DESTROY: not built")
    deallocate(cpasc%nd_dst)
    deallocate(cpasc%nd_src)
    deallocate(cpasc%l_con%cpy)
    deallocate(cpasc%r_con%snd)
    deallocate(cpasc%r_con%rcv)
    deallocate(cpasc%r_con%str)
    deallocate(cpasc%r_con%rtr)
    cpasc%dim = 0
  end subroutine copyassoc_destroy

  function copyassoc_check(cpasc, la_dst, la_src, nd_dst, nd_src) result(r)
    logical                     :: r
    type(copyassoc), intent(in) :: cpasc
    type(layout), intent(in)    :: la_src, la_dst
    logical, intent(in)         :: nd_dst(:), nd_src(:)
    r =         associated(cpasc%lap_dst, la_dst%lap)
    r = r .and. associated(cpasc%lap_src, la_src%lap)
    r = r .and. all(cpasc%nd_dst .eqv. nd_dst)
    r = r .and. all(cpasc%nd_src .eqv. nd_src)
  end function copyassoc_check

  function layout_copyassoc(la_dst, la_src, nd_dst, nd_src) result(r)
    type(copyassoc)                :: r
    type(layout),    intent(inout) :: la_dst
    type(layout),    intent(in)    :: la_src
    logical,         intent(in)    :: nd_dst(:), nd_src(:)
    type(copyassoc), pointer       :: cp
    !
    ! Do we have one stored?
    !
    cp => the_copyassoc_head
    do while ( associated(cp) )
       if ( copyassoc_check(cp, la_dst, la_src, nd_dst, nd_src) ) then
          r = cp
          return
       end if
       cp => cp%next
    end do
    !
    ! Gotta build one.
    !
    allocate(cp)
    call copyassoc_build(cp, la_dst, la_src, nd_dst, nd_src)
    cp%next => the_copyassoc_head
    the_copyassoc_head => cp
    r = cp
  end function layout_copyassoc

  function copyassoc_built_q(cpasc) result(r)
    logical :: r
    type(copyassoc), intent(in) :: cpasc
    r = cpasc%dim /= 0
  end function copyassoc_built_q

  subroutine init_box_hash_bin(la, crsn)
    type(layout), intent(inout) :: la
    integer, intent(in), optional :: crsn
    type(boxarray) :: ba
    integer, dimension(MAX_SPACEDIM) :: ext, vsz
    integer :: dm, i, j, k, n
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
    bx = boxarray_bbox(ba)
    cbx = coarsen(bx, lcrsn)
    la%lap%plo = 0; la%lap%plo(1:dm) = lwb(cbx)
    la%lap%phi = 0; la%lap%phi(1:dm) = upb(cbx)
    la%lap%vshft = int_coarsen(vsz, lcrsn+1)
    allocate(la%lap%bins(la%lap%plo(1):la%lap%phi(1),la%lap%plo(2):la%lap%phi(2),la%lap%plo(3):la%lap%phi(3)))
    bins => la%lap%bins
    do k = la%lap%plo(3), la%lap%phi(3)
       do j = la%lap%plo(2), la%lap%phi(2)
          do i = la%lap%plo(1), la%lap%phi(1)
             allocate(bins(i,j,k)%iv(0))
          end do
       end do
    end do
    do n = 1, nboxes(ba)
       ext = 0; ext(1:dm) = int_coarsen(lwb(get_box(ba,n)), lcrsn)
       if ( .not. contains(cbx, ext(1:dm)) ) then
          call bl_error("BUILD_BOX_HASH_BIN: Not Contained!")
       end if
       sz = size(bins(ext(1),ext(2),ext(3))%iv)
       allocate(ipv(sz+1))
       ipv(1:sz) = bins(ext(1),ext(2),ext(3))%iv(1:sz)
       ipv(sz+1) = n
       deallocate(bins(ext(1),ext(2),ext(3))%iv)
       bins(ext(1),ext(2),ext(3))%iv => ipv
    end do
    call destroy(bpt)
  end subroutine init_box_hash_bin

  function layout_get_box_intersector(la, bx) result(bi)
    type(box_intersector), pointer :: bi(:)
    type(layout), intent(inout) :: la
    type(box), intent(in) :: bx
    type(box_hash_bin), pointer :: bins(:,:,:)
    integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
    integer :: dm
    type(box) :: bx1
    integer :: i, j, k, n
    type(boxarray) :: ba
    integer, parameter :: MAX_BI = 100
    integer :: cnt
    type(box_intersector) :: tbi(MAX_BI)
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
                do n = 1, size(bins(i,j,k)%iv)
                   bx1 = intersection(bx, ba%bxs(bins(i,j,k)%iv(n)))
                   if ( empty(bx1) ) cycle
                   cnt = cnt + 1
                   tbi(cnt)%i = bins(i,j,k)%iv(n)
                   tbi(cnt)%bx = bx1
                end do
             end do
          end do
       end do
    case (2)
       do j = max(lo(2)-la%lap%vshft(2)-1,la%lap%plo(2)), min(hi(2)+la%lap%vshft(2), la%lap%phi(2))
          do i = max(lo(1)-la%lap%vshft(1)-1,la%lap%plo(1)), min(hi(1)+la%lap%vshft(1), la%lap%phi(1))
             do n = 1, size(bins(i,j,k)%iv)
                bx1 = intersection(bx, ba%bxs(bins(i,j,k)%iv(n)))
                if ( empty(bx1) ) cycle
                cnt = cnt + 1
                tbi(cnt)%i = bins(i,j,k)%iv(n)
                tbi(cnt)%bx = bx1
             end do
          end do
       end do
    case (1)
       do i = max(lo(1)-la%lap%vshft(1)-1,la%lap%plo(1)), min(hi(1)+la%lap%vshft(1), la%lap%phi(1))
          do n = 1, size(bins(i,j,k)%iv)
             bx1 = intersection(bx, ba%bxs(bins(i,j,k)%iv(n)))
             if ( empty(bx1) ) cycle
             cnt = cnt + 1
             tbi(cnt)%i = bins(i,j,k)%iv(n)
             tbi(cnt)%bx = bx1
          end do
       end do
    end select
    allocate(bi(cnt))
    bi(1:cnt) = tbi(1:cnt)
  end function layout_get_box_intersector

  end module layout_module
