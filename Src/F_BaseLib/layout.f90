module layout_module

  use parallel
  use boxarray_module
  use knapsack_module
  use bl_mem_stat_module

  implicit none

  integer, private, parameter :: LA_UNDF = 0
  integer, private, parameter :: LA_BASE = 1
  integer, private, parameter :: LA_CRSN = 2
  integer, private, parameter :: LA_PCHD = 3
  integer, private, parameter :: LA_PPRT = 4
  integer, private, parameter :: LA_DERV = 5

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

!  type boxother
!     integer   :: no = 0
!     type(box) :: bx
!  end type boxother

!  type boxinters
!     integer :: ddbx = 0
!     type(boxother), pointer :: bxo_r(:) => Null()
!  end type boxinters

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

  type copyassoc
     integer :: dim = 0                      ! spatial dimension 1, 2, or 3
     logical, pointer :: nd_dst(:) => Null() ! dst nodal flag
     logical, pointer :: nd_src(:) => Null() ! src nodal flag
     type(local_conn)  :: l_con
     type(remote_conn) :: r_con
     type(copyassoc), pointer :: next => Null()
     type(layout_rep), pointer :: lap_dst => Null()
     type(layout_rep), pointer :: lap_src => Null()
  end type copyassoc

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
     type(copyassoc), pointer :: cpasc => Null()
     type(coarsened_layout), pointer :: crse_la => Null()
     type(pn_layout), pointer :: pn_children => Null()
     type(derived_layout), pointer :: dlay => Null()
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
  end interface

  interface get_pd
     module procedure layout_get_pd
  end interface

  private layout_next_id
  private layout_rep_build
  private layout_rep_destroy

  type(mem_stats), private, save :: bxa_ms
  type(mem_stats), private, save :: la_ms

contains

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

  subroutine layout_rep_build(lap, ba, pd, pmask)
    type(layout_rep), intent(out) :: lap
    type(boxarray), intent(in) :: ba
    logical :: pmask(:)
    type(box), intent(in) :: pd
    integer :: i, j

    call boxarray_build_copy(lap%bxa, ba)
    lap%dim = lap%bxa%dim
    lap%nboxes = lap%bxa%nboxes
    lap%id    = layout_next_id()
    lap%pd    = pd
    allocate(lap%pmask(lap%dim))
    lap%pmask = pmask

    allocate(lap%prc(lap%nboxes))
    ! call layout_roundrobin(lap%prc, ba%bxs)
    call layout_knapsack(lap%prc, ba%bxs)

  end subroutine layout_rep_build

  recursive subroutine layout_rep_destroy(lap, la_type)
    type(layout_rep), pointer :: lap
    integer, intent(in) :: la_type
    type(coarsened_layout), pointer :: clp, oclp
    type(pn_layout), pointer :: pnp, opnp
    type(derived_layout), pointer :: dla, odla
    type(boxassoc), pointer :: bxa, obxa
    type(copyassoc), pointer :: cpa, ocpa
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
    bxa => lap%bxasc
    do while ( associated(bxa) )
       obxa => bxa%next
       call boxassoc_destroy(bxa)
       deallocate(bxa)
       bxa => obxa
    end do
    cpa => lap%cpasc
    do while ( associated(cpa) )
       ocpa => cpa%next
       call copyassoc_destroy(cpa)
       deallocate(cpa)
       cpa => ocpa
    end do
    deallocate(lap)
  end subroutine layout_rep_destroy

  subroutine layout_build_ba(la, ba, pd, pmask)
    type(layout), intent(out) :: la
    type(boxarray), intent(in) :: ba
    type(box), intent(in), optional :: pd
    logical, optional :: pmask(:)
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
    call layout_rep_build(la%lap, ba, lpd, lpmask)
  end subroutine layout_build_ba

  subroutine layout_destroy(la)
    type(layout), intent(inout) :: la
    if ( la%la_type /= LA_BASE ) call bl_error("LAYOUT_DESTROY: confused")
    call layout_rep_destroy(la%lap, LA_BASE)
!   deallocate(la%lap)
  end subroutine layout_destroy

  subroutine layout_build_pn(lapn, la, ba, rr)
    type(layout), intent(out)   :: lapn
    type(layout), intent(inout) :: la
    type(boxarray), intent(in) :: ba
    integer, intent(in) :: rr(:)
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
    call layout_rep_build(pla%la%lap, ba, rpd, la%lap%pmask)

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
    integer :: un
    un = unit_stdout(unit)
    call unit_skip(unit, skip)
    write(unit=un, fmt='("LAYOUT")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(": ", A)') str
    else
       write(unit=un, fmt='()')
    end if
    if ( .not. associated(la%lap) ) then
       call unit_skip(unit, skip)
       write(unit=un, fmt='(" NOT ASSOCIATED")')
    else
       call unit_skip(unit, skip)
       write(unit=un, fmt='(" ID = ",i0)', advance = 'no') la%lap%id
       call unit_skip(unit, skip)
       write(unit=un, fmt='(" DIM     = ",i2)') la%lap%dim
       call unit_skip(unit, skip)
       write(unit=un, fmt='(" NBOXES  = ",i2)') la%lap%nboxes
       call unit_skip(unit, skip)
       write(unit=un, fmt='(" PD      = ",i2)', advance = 'no')
       call print(la%lap%pd)
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
    ! didn't find; so have to go looking for it
    allocate (bp)
    call boxassoc_build(bp, la%lap, ng, nodal, cross)
    bp%next => la%lap%bxasc
    la%lap%bxasc => bp
    r = bp
  end function layout_boxassoc

  function boxassoc_built_q(bxasc) result(r)
    logical :: r
    type(boxassoc), intent(in) :: bxasc
    r = bxasc%dim /= 0
  end function boxassoc_built_q

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

    lcross = .false.; if ( present(cross) ) lcross = cross

    call boxarray_box_boundary_n(tba, box_nodalize(b,nodal), ng)

    if ( lcross ) then
       call boxarray_box_corners(cba, box_nodalize(b,nodal), ng)
       call boxarray_diff(tba, cba)
    end if

    shfts = 0

    call box_periodic_shift(dmn, b, nodal, pmask, ng, shft, bxs, cnt)

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

  end subroutine boxarray_bndry_periodic

  subroutine boxassoc_build(bxasc, lap, ng, nodal, cross)

    integer,          intent(in)           :: ng
    logical,          intent(in), optional :: nodal(:)
    type(layout_rep), intent(in), target   :: lap
    type(boxassoc),   intent(inout)        :: bxasc
    logical,          intent(in), optional :: cross

    integer              :: i, j, ii, pv, rpv, spv, pi_r, pi_s, pcnt_r, pcnt_s
    integer              :: shft(2*3**lap%dim,lap%dim), sh(MAX_SPACEDIM+1)
    type(box)            :: abx, bx
    type(boxarray)       :: bxa, bxai
    type(layout)         :: la
    integer              :: lcnt, lcnt_r, li_r, cnt_r, cnt_s, i_r, i_s
    integer, allocatable :: pvol(:,:), ppvol(:,:), parr(:,:)

    if ( built_q(bxasc) ) call bl_error("BOXASSOC_BUILD: alread built")

    la%lap => lap

    bxa = get_boxarray(la)

    bxasc%dim    = bxa%dim
    bxasc%grwth  = ng
    bxasc%nboxes = bxa%nboxes

    allocate(bxasc%nodal(bxasc%dim))
    allocate(parr(0:parallel_nprocs()-1,2))
    allocate(pvol(0:parallel_nprocs()-1,2))
    allocate(ppvol(0:parallel_nprocs()-1,2))

    bxasc%nodal = .false.; if ( present(nodal) ) bxasc%nodal = nodal

    parr = 0; pvol = 0; lcnt_r = 0; cnt_r = 0; cnt_s = 0
    !
    ! We here consider all copies I <- J.
    !
    do i = 1, bxa%nboxes
       lcnt = 0
       call boxarray_bndry_periodic(bxai, lap%pd, bxa%bxs(i), bxasc%nodal, lap%pmask, ng, shft, cross)
       do j = 1, bxa%nboxes
          if ( remote(la,i) .and. remote(la,j) ) cycle
          bx = box_nodalize(get_box(bxa, j), nodal)
          do ii = 1, bxai%nboxes
             abx = intersection(bx, bxai%bxs(ii))
             if ( .not. empty(abx) ) then
                if ( local(la,i) .and. local(la, j) ) then
                   lcnt   = lcnt   + 1
                   lcnt_r = lcnt_r + 1
                else if ( local(la, j) ) then
                   cnt_s               = cnt_s + 1
                   parr(lap%prc(i), 2) = parr(lap%prc(i), 2) + 1
                   pvol(lap%prc(i), 2) = pvol(lap%prc(i), 2) + volume(abx)
                else if ( local(la, i) ) then
                   cnt_r               = cnt_r + 1
                   parr(lap%prc(j), 1) = parr(lap%prc(j), 1) + 1
                   pvol(lap%prc(j), 1) = pvol(lap%prc(j), 1) + volume(abx)
                end if
             end if
          end do
       end do
       call destroy(bxai)
    end do
    !
    ! Fill in the boxassoc structure.
    !
    bxasc%l_con%ncpy = lcnt_r
    bxasc%r_con%nsnd = cnt_s
    bxasc%r_con%nrcv = cnt_r
    allocate(bxasc%l_con%cpy(lcnt_r))
    allocate(bxasc%r_con%snd(cnt_s))
    allocate(bxasc%r_con%rcv(cnt_r))
    li_r = 1; i_r = 1; i_s = 1

    do i = 1, bxa%nboxes
       call boxarray_bndry_periodic(bxai, lap%pd, bxa%bxs(i), bxasc%nodal, lap%pmask, ng, shft, cross)
       do j = 1, bxa%nboxes
          if ( remote(la,i) .and. remote(la,j) ) cycle
          bx = box_nodalize(get_box(bxa, j), nodal)
          do ii = 1, bxai%nboxes
             abx = intersection(bx, bxai%bxs(ii))
             if ( .not. empty(abx) ) then
                if ( local(la,i) .and. local(la, j) ) then
                      bxasc%l_con%cpy(li_r)%nd  = i
                      bxasc%l_con%cpy(li_r)%ns  = j
                      bxasc%l_con%cpy(li_r)%sbx = abx
                      bxasc%l_con%cpy(li_r)%dbx = shift(abx,-shft(ii,:))
                      li_r                      = li_r + 1
                   else if ( local(la,j) ) then
                      bxasc%r_con%snd(i_s)%nd  = i
                      bxasc%r_con%snd(i_s)%ns  = j
                      bxasc%r_con%snd(i_s)%sbx = abx
                      bxasc%r_con%snd(i_s)%dbx = shift(abx,-shft(ii,:))
                      bxasc%r_con%snd(i_s)%pr  = get_proc(la, i)
                      bxasc%r_con%snd(i_s)%s1  = volume(abx)
                      i_s                      = i_s + 1
                   else if ( local(la,i) ) then
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
                end if
          end do
       end do
       call destroy(bxai)
    end do
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

  end subroutine boxassoc_build

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
             write(unit=un, fmt='(" ",i4,":(",i4,"<-",i4,"): ")', advance = 'no') &
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
             write(unit=un, fmt='(" ",i4,":(",i4,"<-",i4,"):",i4,":",i4,":",i4,":",i4,":",4(1x,i4),": ")', &
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
             write(unit=un, fmt='(" ",i4,":(",i4,"<-",i4,"):",i4,":",i4,":",i4,":",i4,":",4(1x,i4),": ")', &
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

    integer              :: i, j, pv, rpv, spv, pi_r, pi_s, pcnt_r, pcnt_s
    integer              :: sh(MAX_SPACEDIM+1)
    type(box)            :: bx
    type(boxarray)       :: bxa_src, bxa_dst
    integer              :: lcnt, lcnt_r, li_r, cnt_r, cnt_s, i_r, i_s
    integer, allocatable :: pvol(:,:), ppvol(:,:), parr(:,:)

    if ( built_q(cpasc) ) call bl_error("COPYASSOC_BUILD: alread built")

    bxa_src = get_boxarray(la_src)
    bxa_dst = get_boxarray(la_dst)

    cpasc%dim     = bxa_src%dim
    cpasc%lap_dst => la_dst%lap
    cpasc%lap_src => la_src%lap

    allocate(cpasc%nd_dst(la_dst%lap%dim))
    allocate(cpasc%nd_src(la_src%lap%dim))
    allocate(parr(0:parallel_nprocs()-1,2))
    allocate(pvol(0:parallel_nprocs()-1,2))
    allocate(ppvol(0:parallel_nprocs()-1,2))

    cpasc%nd_dst = nd_dst
    cpasc%nd_src = nd_src

    parr = 0; pvol = 0; lcnt_r = 0; cnt_r = 0; cnt_s = 0

    do i = 1, bxa_dst%nboxes
       lcnt = 0
       do j = 1, bxa_src%nboxes
          if ( remote(la_dst,i) .and. remote(la_src,j) ) cycle
          bx = intersection(box_nodalize(get_box(bxa_dst,i),nd_dst), box_nodalize(get_box(bxa_src,j),nd_src))
          if ( .not. empty(bx) ) then
             if ( local(la_dst,i) .and. local(la_src,j) ) then
                lcnt   = lcnt   + 1
                lcnt_r = lcnt_r + 1
             else if ( local(la_src,j) ) then
                cnt_s               = cnt_s + 1
                parr(la_dst%lap%prc(i), 2) = parr(la_dst%lap%prc(i), 2) + 1
                pvol(la_dst%lap%prc(i), 2) = pvol(la_dst%lap%prc(i), 2) + volume(bx)
             else if ( local(la_dst,i) ) then
                cnt_r               = cnt_r + 1
                parr(la_src%lap%prc(j), 1) = parr(la_src%lap%prc(j), 1) + 1
                pvol(la_src%lap%prc(j), 1) = pvol(la_src%lap%prc(j), 1) + volume(bx)
             end if
          end if
       end do
    end do
    !
    ! Fill in the copyassoc structure.
    !
    cpasc%l_con%ncpy = lcnt_r
    cpasc%r_con%nsnd = cnt_s
    cpasc%r_con%nrcv = cnt_r
    allocate(cpasc%l_con%cpy(lcnt_r))
    allocate(cpasc%r_con%snd(cnt_s))
    allocate(cpasc%r_con%rcv(cnt_r))
    li_r = 1; i_r = 1; i_s = 1

    do i = 1, bxa_dst%nboxes
       do j = 1, bxa_src%nboxes
          if ( remote(la_dst,i) .and. remote(la_src,j) ) cycle
          bx = intersection(box_nodalize(get_box(bxa_dst,i),nd_dst), box_nodalize(get_box(bxa_src,j),nd_src))
          if ( .not. empty(bx) ) then
             if ( local(la_dst,i) .and. local(la_src,j) ) then
                cpasc%l_con%cpy(li_r)%nd  = i
                cpasc%l_con%cpy(li_r)%ns  = j
                cpasc%l_con%cpy(li_r)%sbx = bx
                cpasc%l_con%cpy(li_r)%dbx = bx
                li_r                      = li_r + 1
             else if ( local(la_src,j) ) then
                cpasc%r_con%snd(i_s)%nd  = i
                cpasc%r_con%snd(i_s)%ns  = j
                cpasc%r_con%snd(i_s)%sbx = bx
                cpasc%r_con%snd(i_s)%dbx = bx
                cpasc%r_con%snd(i_s)%pr  = get_proc(la_dst,i)
                cpasc%r_con%snd(i_s)%s1  = volume(bx)
                i_s                      = i_s + 1
             else if ( local(la_dst,i) ) then
                cpasc%r_con%rcv(i_r)%nd  = i
                cpasc%r_con%rcv(i_r)%ns  = j
                cpasc%r_con%rcv(i_r)%sbx = bx
                cpasc%r_con%rcv(i_r)%dbx = bx
                cpasc%r_con%rcv(i_r)%pr  = get_proc(la_src,j)
                sh                       = 1
                sh(1:cpasc%dim)          = extent(bx)
                cpasc%r_con%rcv(i_r)%sh  = sh
                i_r                      = i_r + 1
             end if
          end if
       end do
    end do
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
    ! Do we have one stored in the "dst" layout?
    !
    cp => la_dst%lap%cpasc
    do while ( associated(cp) )
       if ( copyassoc_check(cp, la_dst, la_src, nd_dst, nd_src) ) then
          r = cp
          return
       end if
       cp => cp%next
    end do
    !
    ! Gotta build one then store in the "dst" layout.
    !
    allocate (cp)
    call copyassoc_build(cp, la_dst, la_src, nd_dst, nd_src)
    cp%next => la_dst%lap%cpasc
    la_dst%lap%cpasc => cp
    r = cp
  end function layout_copyassoc

  function copyassoc_built_q(cpasc) result(r)
    logical :: r
    type(copyassoc), intent(in) :: cpasc
    r = cpasc%dim /= 0
  end function copyassoc_built_q

end module layout_module
