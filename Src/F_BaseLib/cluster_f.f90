module cluster_module

  use bl_types
  use box_module
  use list_box_module
  use boxarray_module
  use multifab_module

  implicit none

  private

  logical, private, save :: verbose = .false.
  integer, private, save :: cut_threshold = 2
  real(dp_t), private, save :: beta = 1.05_dp_t

  integer        , parameter, private :: minwidth_default        = 8
  integer        , parameter, private :: blocking_factor_default = 8
  integer        , parameter, private :: ref_ratio_default       = 2
  real(kind=dp_t), parameter, private :: min_eff_default  = 0.7_dp_t

  integer        , save :: minwidth  = minwidth_default
  integer        , save :: ref_ratio = ref_ratio_default
  integer        , save :: blocking_factor = blocking_factor_default
  real(kind=dp_t), save :: min_eff = min_eff_default

  public :: cluster_set_verbose
  public :: cluster_get_verbose
  public :: cluster_set_cut_threshold
  public :: cluster_get_cut_threshold
  public :: cluster_set_beta
  public :: cluster_get_beta
  public :: cluster_set_minwidth
  public :: cluster_set_blocking_factor
  public :: cluster_set_min_eff
  public :: cluster

  public :: tagboxes_coarsen

  interface cluster
     module procedure cls_3d_mf
  end interface

contains

  subroutine cluster_set_minwidth(val)
    integer, intent(in) :: val
    minwidth = val
  end subroutine cluster_set_minwidth

  subroutine cluster_set_blocking_factor(val)
    integer, intent(in) :: val
    blocking_factor = val
  end subroutine cluster_set_blocking_factor

  subroutine cluster_set_min_eff(val)
    real(dp_t), intent(in) :: val
    min_eff = val
  end subroutine cluster_set_min_eff

  subroutine cluster_set_verbose(val)
    logical, intent(in) :: val
    verbose = val
  end subroutine cluster_set_verbose
  function cluster_get_verbose() result(r)
    logical :: r
    r = verbose
  end function cluster_get_verbose

  subroutine cluster_set_cut_threshold(val)
    integer, intent(in) :: val
    cut_threshold= val
  end subroutine cluster_set_cut_threshold
  function cluster_get_cut_threshold() result(r)
    integer :: r
    r = cut_threshold
  end function cluster_get_cut_threshold

  subroutine cluster_set_beta(val)
    real(dp_t), intent(in) :: val
    beta = val
  end subroutine cluster_set_beta
  function cluster_get_beta() result(r)
    real(dp_t)  :: r
    r = beta
  end function cluster_get_beta

  subroutine cls_3d_mf(boxes, tagboxes, buf_wid, overall_eff)
    use bl_error_module
    use bl_prof_module

    type(boxarray),  intent(out)           :: boxes
    type(lmultifab), intent(in )           :: tagboxes
    integer,         intent(in )           :: buf_wid
    real(dp_t),      intent(out), optional :: overall_eff

    type(layout)                 :: la, cla, la_buf
    type(list_box)               :: lboxes
    type(lmultifab)              :: cbuf, buf, lbuf
    logical, pointer             :: mt(:,:,:,:), mb(:,:,:,:)
    integer                      :: b, dm, num_flag, i, k, bxcnt
    integer                      :: ratio, lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
    integer, allocatable         :: iprocs(:), bxs(:)
    integer(kind=ll_t)           :: bboxinte
    type(box)                    :: bx,buf_pd
    real(dp_t)                   :: bx_eff
    type(list_box_node), pointer :: bln

    type(bl_prof_timer), save :: bpt

    call build(bpt, "cluster")

    if ( nghost(tagboxes) /= 0 ) &
       call bl_error("CLUSTER: tagboxes must have NG = 0")

    if (mod(minwidth,blocking_factor) .ne. 0) &
         call bl_error("CLUSTER: minwidth must be an integer multiple of blocking_factor")

    if ( min_eff < 0 .or. min_eff > 1.0_dp_t ) &
       call bl_error("CLUSTER: min_eff must be >= 0 and <= 1: ", min_eff)

    if ( buf_wid < 0 ) &
       call bl_error("CLUSTER: buf_wid must be >= 0: ", buf_wid)

    dm       = get_dim(tagboxes)
    num_flag = lmultifab_count(tagboxes)

    if ( num_flag == 0 ) then
       if ( present(overall_eff) ) overall_eff = 0
       call destroy(bpt)
       return
    end if

    la_buf = get_layout(tagboxes)

    call build(buf, la_buf, nc = 1, ng = buf_wid)
    call setval(buf, .false., all = .true.)
    !
    ! Buffer the cells.
    !
    do b = 1, nfabs(buf)
       mt => dataptr(tagboxes, b, get_box(tagboxes,b))
       mb => dataptr(buf, b)
       select case (dm)
       case (1)
          call buffer_1d(mb(:,1,1,1), mt(:,1,1,:), buf_wid)
       case (2)
          call buffer_2d(mb(:,:,1,1), mt(:,:,1,:), buf_wid)
       case (3)
          call buffer_3d(mb(:,:,:,1), mt(:,:,:,:), buf_wid)
       end select
    end do

    ! modify the valid region using logical or with ghost cells that
    ! share the same physical space
    call lmultifab_sum_boundary(buf)

    ! make ghost cells consistent with valid region
    call lmultifab_fill_boundary(buf)

    ! Remove any tags outside the problem domain.
    buf_pd = get_pd(la_buf)
    call pd_mask(buf,buf_pd)

    ratio = max(blocking_factor / ref_ratio, 1)

    call tagboxes_coarsen(buf, cbuf, ratio)
    !
    ! Remove any coarsened tags outside the problem domain.
    !
    call pd_mask(cbuf,coarsen(buf_pd,ratio))

    call destroy(buf)

    cla = get_layout(cbuf)

    bx = coarsen(get_pd(la_buf),ratio)
    !
    ! The cluster algorithm is inherently serial.
    ! We'll set up the problem to do on the IO processor.
    !
    allocate(iprocs(nboxes(cla)))

    iprocs = parallel_IOProcessorNode()

    call build(la, get_boxarray(cla), bx, get_pmask(cla), explicit_mapping = iprocs)

    call build(lbuf, la, nc = 1)

    call copy(lbuf, cbuf)  ! This is a parallel copy.

    call destroy(cbuf)
    call destroy(cla)

    if ( parallel_IOProcessor() ) then
       call owner_mask(lbuf)
       call cluster_mf_private_recursive(lboxes, bx, bx_eff, lbuf)
       bxcnt = size(lboxes)

       if (ratio > 1) then
          bln => begin(lboxes)
          do while ( associated(bln) )
             bln = refine(value(bln),ratio)
             bln => next(bln)
          end do
       end if
    end if

    call destroy(lbuf)
    call destroy(la)
    !
    ! We now must broadcast the info in "lboxes" back to all CPUs.
    !
    ! First broadcast the number of boxes to expect.
    !
    call parallel_bcast(bxcnt, parallel_IOProcessorNode())
    !
    ! We'll pass 2*MAX_SPACEDIM integers for each box we need to send.
    !
    allocate(bxs(2*MAX_SPACEDIM*bxcnt))

    if ( parallel_IOProcessor() ) then
       i   =  1
       bln => begin(lboxes)
       do while (associated(bln))
          bx = value(bln)
          bxs(i:i+MAX_SPACEDIM-1) = bx%lo
          i = i + MAX_SPACEDIM
          bxs(i:i+MAX_SPACEDIM-1) = bx%hi
          i = i + MAX_SPACEDIM
          bln => next(bln)
       end do
    end if

    call parallel_bcast(bxs, parallel_IOProcessorNode())

    if ( .not. parallel_IOProcessor() ) then
       i = 1
       do k = 1, bxcnt
          lo = bxs(i:i+MAX_SPACEDIM-1)
          i  = i + MAX_SPACEDIM
          hi = bxs(i:i+MAX_SPACEDIM-1)
          i  = i + MAX_SPACEDIM
          call push_back(lboxes, make_box(lo(1:dm), hi(1:dm)))
       end do
    end if

    call build(boxes, lboxes)

    call destroy(lboxes)

    if ( present(overall_eff) ) then
       bboxinte    = volume(boxes)
       overall_eff = real(count(buf,all=.true.),dp_t) / real(bboxinte, dp_t)
    end if

    call destroy(bpt)

  contains

    subroutine pd_mask(mask,pd)

      type(lmultifab), intent(inout) :: mask
      type(box)      , intent(in   ) :: pd
      integer :: i, j
      type(box) :: bxi
      type(boxarray) :: ba
      logical, pointer :: lp(:,:,:,:)

      do i = 1, nfabs(mask)
         bxi = get_pbox(mask, i)
         call boxarray_box_diff(ba, bxi, pd)
         do j = 1, nboxes(ba)
            lp => dataptr(mask, i, get_box(ba,j))
            lp = .false.
         end do
         call destroy(ba)
      end do
    end subroutine pd_mask

    subroutine owner_mask(mask)
      type(lmultifab), intent(inout) :: mask
      integer :: i, j
      type(box) :: bxi, bxj, bxij
      logical, pointer :: lp(:,:,:,:)

      do i = 1, nfabs(mask)
         bxi = get_pbox(mask, i)
         do j = 1, i-1
            bxj = get_pbox(mask, j)
            bxij = intersection(bxi, bxj)
            if ( empty(bxij) ) cycle
            lp => dataptr(mask, j, bxij)
            lp = .false.
         end do
      end do
    end subroutine owner_mask

    subroutine buffer_1d(bb, tt, ng)
      integer, intent(in) :: ng
      logical, intent(in) :: tt(:,:)
      logical, intent(out) :: bb(1-ng:)
      integer :: i
      do i = 1, size(tt,1)
         if ( any(tt(i,:)) ) bb(i-ng:i+ng) = .true.
      end do
    end subroutine buffer_1d

    subroutine buffer_2d(bb, tt, ng)
      integer, intent(in) :: ng
      logical, intent(in) :: tt(:,:,:)
      logical, intent(out) :: bb(1-ng:,1-ng:)
      integer :: i, j
      do j = 1, size(tt,2); do i = 1, size(tt,1)
         if ( any(tt(i,j,:)) ) bb(i-ng:i+ng,j-ng:j+ng) = .true.
      end do; end do
    end subroutine buffer_2d

    subroutine buffer_3d(bb, tt, ng)
      integer, intent(in) :: ng
      logical, intent(in) :: tt(:,:,:,:)
      logical, intent(out) :: bb(1-ng:,1-ng:,1-ng:)
      integer :: i, j, k
      do k = 1, size(tt,3); do j = 1, size(tt,2); do i = 1, size(tt,1)
         if ( any(tt(i,j,k,:)) ) bb(i-ng:i+ng,j-ng:j+ng,k-ng:k+ng) = .true.
      end do; end do; end do
    end subroutine buffer_3d

  end subroutine cls_3d_mf

  subroutine map_periodic(mask)

    type(lmultifab), intent(inout) :: mask

    integer                        :: i, j, ii, jj, cnt, shft(3**get_dim(mask), get_dim(mask)), li, lj
    type(box)                      :: pd, bxs(3**get_dim(mask)), bx_from, bx_to, ibx
    logical                        :: pmask(get_dim(mask))
    logical, pointer               :: ap(:,:,:,:), bp(:,:,:,:)
    logical, allocatable           :: pt(:,:,:,:)
    integer, parameter             :: tag = 2121
    type(box_intersector), pointer :: bi(:)
    type(layout)                   :: la
    type(bl_prof_timer),   save    :: bpt

    call build(bpt, "map_periodic")

    la    = get_layout(mask)
    pd    = get_pd(la)
    pmask = get_pmask(la)

    if ( .not. any(pmask) ) then
       call destroy(bpt)
       return
    end if

    do i = 1, nboxes(mask%la)

       ibx = box_nodalize(get_box(mask%la,i),mask%nodal)

       call box_periodic_shift(pd, ibx, nodal_flags(mask), pmask, nghost(mask), shft, cnt, bxs)

       do jj = 1, cnt
          bi => layout_get_box_intersector(la, bxs(jj))

          do ii = 1, size(bi)
             j = bi(ii)%i

             if ( remote(mask%la,i) .and. remote(mask%la,j) ) cycle

             bx_to   = bi(ii)%bx
             bx_from = shift(bx_to, -shft(jj,:))

             if ( local(mask%la,i) .and. local(mask%la,j) ) then
                li =  local_index(mask,i)
                lj =  local_index(mask,j)
                ap => dataptr(mask,lj,bx_to)
                bp => dataptr(mask,li,bx_from)
                ap = ap .or. bp
             else if ( local(mask%la,i) ) then
                !
                ! We own index i.  Got to send it to processor owning index j.
                !
                li =  local_index(mask,i)
                bp => dataptr(mask,li,bx_from)
                call parallel_send(bp, get_proc(get_layout(mask),j), tag)
             else
                !
                ! We own index j.  Got to get index i data from processor owning it.
                !
                lj =  local_index(mask,j)
                ap => dataptr(mask,lj,bx_to)
                allocate(pt(1:size(ap,1),1:size(ap,2),1:size(ap,3),1))
                call parallel_recv(pt, get_proc(get_layout(mask),i), tag)
                ap = ap .or. pt
                deallocate(pt)
             end if
          end do

          deallocate(bi)
       end do
    end do

    call destroy(bpt)

  end subroutine map_periodic

  function cluster_tally(tagboxes, bx, tx, ty, tz, ll, hh) result(r)
    type(box) :: r
    type(lmultifab), intent(in) :: tagboxes
    type(box), intent(in) :: bx
    integer, intent(in) :: ll(:), hh(:)
    integer, intent(out) :: tx(ll(1):), ty(ll(2):), tz(ll(3):)
    integer :: ll1(3), hh1(3)
    integer :: n, ii, jj, kk, i, j, k
    logical, pointer :: tp(:,:,:,:)
    integer :: dm
    type(box) :: bx1,pbx
    integer :: llo(3), hho(3)

    dm = get_dim(tagboxes)

    tx = 0
    ty = 0
    tz = 0

    do n = 1, nfabs(tagboxes)
       pbx = grow(box_nodalize(get_box(tagboxes%la,n),tagboxes%nodal),tagboxes%ng)
       bx1 = intersection(pbx, bx)
       if ( empty(bx1) ) cycle
       tp => dataptr(tagboxes, n, bx1)

       ll1 = 1; ll1(1:dm) = lwb(bx1)
       hh1 = 1; hh1(1:dm) = extent(bx1)
       do i = 1, hh1(1)
          ii = ll1(1) + i - 1
          tx(ii) = tx(ii) + count(tp(i,:,:,1))
       end do

       do j = 1, hh1(2)
          jj = ll1(2) + j - 1
          ty(jj) = ty(jj) + count(tp(:,j,:,1))
       end do

       do k = 1, hh1(3)
          kk = ll1(3) + k - 1
          tz(kk) = tz(kk) + count(tp(:,:,k,1))
       end do
    end do

    call bracket(tx, ll(1), hh(1), llo(1), hho(1))
    call bracket(ty, ll(2), hh(2), llo(2), hho(2))
    call bracket(tz, ll(3), hh(3), llo(3), hho(3))

    if ( all(llo <= hho) ) then
       r = make_box(llo(1:dm), hho(1:dm))
    end if

  contains

    subroutine bracket(sig, ll, hh, llo, hho)
      integer, intent(in) :: ll, hh
      integer, intent(out) :: llo, hho
      integer, intent(in)    :: sig(ll:hh)

      llo = ll
      hho = hh
      do while ( sig(llo) == 0  .and.  hho-llo >= minwidth/blocking_factor )
         llo = llo + 1
      end do

      do while ( sig(hho) ==  0  .and. hho-llo >= minwidth/blocking_factor )
         hho = hho - 1
      end do

    end subroutine bracket

  end function cluster_tally

  recursive subroutine cluster_mf_private_recursive(boxes, bx, bx_eff, tagboxes)
    type(list_box), intent(out) :: boxes
    type(box), intent(in) :: bx
    type(lmultifab), intent(in) ::  tagboxes
    real(dp_t), intent(out) ::  bx_eff ! Only meaningfull for the case of 1 box returned
    type(box) :: bbx
    integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
    type(list_box) :: C1, C2
    type(box) :: b1, b2
    integer, allocatable :: sx(:), sy(:), sz(:)
    real(dp_t) :: b1_eff, b2_eff, c_eff
    integer :: ll(3), hh(3)
    integer :: dm

    dm  = get_dim(tagboxes)

    if ( verbose ) call print(bx, 'in bx')

    if ( empty(bx) ) return

    ll = 1; ll(1:dm) = lwb(bx)
    hh = 1; hh(1:dm) = upb(bx)

    allocate(sx(ll(1):hh(1)), sy(ll(2):hh(2)), sz(ll(3):hh(3)))

    bbx = cluster_tally(tagboxes, bx, sx, sy, sz, ll, hh)

    if ( verbose ) call print(bbx, '  bbx')

    if ( empty(bbx) ) return

    bx_eff = box_eff_mf(tagboxes, bbx)

    if ( verbose ) then
       print *, 'bx_eff ', bx_eff
       print *, 'min_eff ', min_eff
    end if

    if ( bx_eff == 0.0_dp_t) return

    if ( bx_eff >= min_eff ) then
       call push_back(boxes, bbx)
    else
       lo = 1; lo(1:dm) = lwb(bbx); hi = 1; hi(1:dm) = upb(bbx)
       call find_split(bbx, b1, b2, &
            sx(lo(1):hi(1)), sy(lo(2):hi(2)), sz(lo(3):hi(3)), lo, hi)
       if ( verbose ) then
          call print(bbx, 'bbx;find_split')
          call print(b1, 'b1')
          call print(b2, 'b2')
       end if
       call cluster_mf_private_recursive(c1, b1, b1_eff, tagboxes)
       call cluster_mf_private_recursive(c2, b2, b2_eff, tagboxes)
       if ( verbose ) then
          print *, 'size(c1) = ', size(c1), 'size(c2) = ', size(c2)
       end if
       if ( size(c1) > 0 .or. size(c2) > 0 ) then
          if ( size(c1) == 1 .and. size(c2) == 1 ) then
             b1 = front(c1)
             b2 = front(c2)
             c_eff = ( dvolume(b1)*b1_eff + dvolume(b2)*b2_eff )/(dvolume(b1) + dvolume(b2))
             if ( verbose ) then
                print *, 'c_eff = ', c_eff, bx_eff, c_eff-beta*bx_eff>0.0_dp_t
             end if
             if ( c_eff > beta*bx_eff ) then
                call push_back(boxes, b1)
                call push_back(boxes, b2)
             else
                call push_back(boxes, bbx)
             end if
             call destroy(c1)
             call destroy(c2)
          else
             boxes = c1
             call splice(boxes, c2)
          end if
       else
          call push_back(boxes, bbx)
       end if
    end if

    contains

      function box_eff_mf(tagboxes, bx) result(r)
        real(dp_t) :: r
        type(box), intent(in) :: bx
        type(lmultifab), intent(in) :: tagboxes
        logical, pointer :: tp(:,:,:,:)
        integer :: n
        type(box) :: bx1, pbx
        r = 0.0_dp_t
        do n = 1, nboxes(tagboxes%la)
           pbx = grow(box_nodalize(get_box(tagboxes%la,n),tagboxes%nodal),tagboxes%ng)
           bx1 = intersection(pbx,bx)
           if ( empty(bx1) ) cycle
           tp => dataptr(tagboxes, n, bx1)
           r = r + real(count(tp),dp_t)
        end do
        r = r/dvolume(bx)
      end function box_eff_mf

  end subroutine cluster_mf_private_recursive

  subroutine find_split(bx, b1, b2, sx, sy, sz, ll, hh)
    type(box), intent(in) ::  bx
    type(box), intent(out) :: b1, b2
    integer, intent(in) :: ll(:), hh(:)
    integer, intent(in) :: sx(ll(1):), sy(ll(2):), sz(ll(3):)
    integer, allocatable :: lx(:), ly(:), lz(:)
    integer :: i, j, k
    integer :: hi(3), ip(3)

    if ( holes(bx, b1, b2, sx, ll(1), hh(1), 1) ) return
    if ( holes(bx, b1, b2, sy, ll(2), hh(2), 2) ) return
    if ( holes(bx, b1, b2, sz, ll(3), hh(3), 3) ) return

    allocate(lx(ll(1):hh(1)),ly(ll(2):hh(2)),lz(ll(3):hh(3)))

    lx = 0
    ly = 0
    lz = 0

    do i = ll(1)+1, hh(1)-1
       lx(i) = sx(i+1)-2*sx(i)+sx(i-1)
    end do

    do j = ll(2)+1, hh(2)-1
       ly(j) = sy(j+1)-2*sy(j)+sy(j-1)
    end do

    do k = ll(3)+1, hh(3)-1
       lz(k) = sz(k+1)-2*sz(k)+sz(k-1)
    end do

    ! Find the inflection points in each direction
    call inflection(lx, ll(1), hh(1), hi(1), ip(1))
    call inflection(ly, ll(2), hh(2), hi(2), ip(2))
    call inflection(lz, ll(3), hh(3), hi(3), ip(3))

    if ( any(ip /= -1) ) then
       if ( maxval(hi) >= CUT_THRESHOLD ) then
          i = maxloc(hi, dim=1)
          call box_chop(bx, b1, b2, i, ip(i))
          if ( verbose ) then
             call print(bx, 'bx;inflection')
             call print(b1, 'b1')
             call print(b2, 'b2')
             print *, 'extents b1,b2 ', extent(b1), extent(b2)
             print *, 'i = ', i
             print *, 'ip ', ip
             print *, 'hi ', hi
             print *, 'lx ', lx
             print *, 'ly ', ly
             print *, 'lz ', lz
          end if
          return
       end if
    end if

    ! final fall back:

    return

    i = maxloc(extent(bx),dim=1)
    if ( extent(bx,i) >= 2*minwidth/blocking_factor ) then
       call box_chop(bx, b1, b2, i, lwb(bx,i) + extent(bx,i)/2)
       if ( verbose ) then
          call print(bx, 'bx;fallback')
          call print(b1, 'b1')
          call print(b2, 'b2')
       end if
    end if

  contains

    function holes(bx, b1, b2, sig, ll, hh, dim) result(r)
      logical :: r
      integer, intent(in) :: ll, hh, dim
      integer, intent(in)    :: sig(ll:hh)
      integer :: i
      type(box), intent(in) :: bx
      type(box), intent(out) :: b1, b2

      r  = .false.
      do i = ll + minwidth/blocking_factor, hh - minwidth/blocking_factor
         if ( sig(i)  ==  0 ) then
            call box_chop(bx, b1, b2, dim, i)
            b2 = grow(b2, -1, dim, -1)
            r  = .true.
            if ( verbose ) then
               call print(bx, 'bx;holes')
               call print(b1, 'b1')
               call print(b2, 'b2')
            end if
            return
         end if
      end do

    end function holes

    subroutine inflection(lp, ll, hh, hiv, infp)
      integer, intent(in) :: ll, hh
      integer, intent(in) :: lp(ll:hh)
      integer, intent(out) :: hiv, infp
      integer :: i, tmp

      hiv  = 0
      infp = -1
      do i = ll + minwidth/blocking_factor, hh - max(minwidth/blocking_factor,1)
         if( (lp(i) > 0 .and. lp(i+1) < 0)  .or. (lp(i) < 0 .and. lp(i+1) > 0) ) then
            tmp = abs(lp(i)-lp(i+1))
            if ( tmp  >  hiv ) then
               hiv = tmp
               infp  = i+1
            end if
         end if
      end do

    end subroutine inflection

  end subroutine find_split

  subroutine tagboxes_coarsen(tagboxes,ctagboxes,ratio)

    type(lmultifab), intent(in   ) :: tagboxes
    type(lmultifab), intent(inout) :: ctagboxes
    integer,         intent(in   ) :: ratio

    integer          :: ii, i, j, k, ic, jc, kc, dm
    integer          :: flo(get_dim(tagboxes)), fhi(get_dim(tagboxes))
    type(layout)     :: cla, cla_grow
    type(boxarray)   :: cba
    logical, pointer :: fp(:,:,:,:), cp(:,:,:,:)
    type(lmultifab)  :: ctagboxes_grow
    type(box)        :: pbf, pbc

    call bl_assert(ncomp(tagboxes) == 1, 'tagboxes should only have one component')
    !
    ! ctagboxes should be an lmultifab on which build() has not been called.
    ! ctagboxes will be built on the appropriately grown & coarsen'd boxarray
    ! and will have no grow cells itself.  We want all the coarsen'd values
    ! to be in ctagboxes valid region.  Callers of this routine need to
    ! destroy both ctagboxes and its layout.
    !
    call boxarray_build_copy(cba, get_boxarray(tagboxes))

    call boxarray_grow(cba, nghost(tagboxes))

    call boxarray_coarsen(cba, ratio)

    dm = get_dim(tagboxes)
    !
    ! I'm playing a little fast & loose here.
    ! I'm assuming all we really need to get right is the mapping.
    !
    call build(cla, cba, boxarray_bbox(cba), explicit_mapping = get_proc(get_layout(tagboxes)))

    call destroy(cba)

    call build(ctagboxes, cla, 1, 0)

    call setval(ctagboxes, .false., all = .true.)

    do ii = 1, nfabs(tagboxes)

       fp  => dataptr(tagboxes,  ii)
       cp  => dataptr(ctagboxes, ii)

       flo = lwb(get_pbox(tagboxes, ii))
       fhi = upb(get_pbox(tagboxes, ii))

       select case (dm)
       case (2)
          do j = flo(2), fhi(2)
             jc = int_coarsen(j,ratio)
             do i = flo(1), fhi(1)
                ic = int_coarsen(i,ratio)
                if ( fp(i,j,1,1) ) cp(ic,jc,1,1) = .true.
             end do
          end do
       case  (3)
          do k = flo(3), fhi(3)
             kc = int_coarsen(k,ratio)
             do j = flo(2), fhi(2)
                jc = int_coarsen(j,ratio)
                do i = flo(1), fhi(1)
                   ic = int_coarsen(i,ratio)
                   if ( fp(i,j,k,1) ) cp(ic,jc,kc,1) = .true.
                end do
             end do
          end do
       end select
    end do

    call boxarray_build_copy(cba, get_boxarray(tagboxes))

    call boxarray_coarsen(cba, ratio)

    call build(cla_grow, cba, boxarray_bbox(cba), explicit_mapping = get_proc(get_layout(tagboxes)))

    call destroy(cba)

    call build(ctagboxes_grow, cla_grow, 1, nghost(tagboxes)/ratio + 1)

    call setval(ctagboxes_grow, .false., all = .true.)

    do ii = 1, nfabs(ctagboxes_grow)

       pbf = get_pbox(ctagboxes, ii)
       pbc = get_pbox(ctagboxes_grow, ii)

       call bl_assert(box_contains(pbc,pbf), 'pbc does not contain pbf')

       fp  => dataptr(ctagboxes_grow,  ii, pbf)
       cp  => dataptr(ctagboxes,       ii, pbf)

       fp = cp
    end do

    call lmultifab_sum_boundary(ctagboxes_grow)

    call lmultifab_fill_boundary(ctagboxes_grow)

    do ii = 1, nfabs(ctagboxes_grow)

       pbf = get_pbox(ctagboxes, ii)
       pbc = get_pbox(ctagboxes_grow, ii)

       fp  => dataptr(ctagboxes_grow,  ii, pbf)
       cp  => dataptr(ctagboxes,       ii, pbf)

       cp = fp
    end do

    call destroy(ctagboxes_grow)
    call destroy(cla_grow)

  end subroutine tagboxes_coarsen

end module cluster_module
