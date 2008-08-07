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

  integer        , parameter, private :: minwidth_default = 4
  real(kind=dp_t), parameter, private :: min_eff_default  = .7

  integer        , save :: minwidth = minwidth_default
  real(kind=dp_t), save :: min_eff = min_eff_default

  public :: cluster_set_verbose
  public :: cluster_get_verbose
  public :: cluster_set_cut_threshold
  public :: cluster_get_cut_threshold
  public :: cluster_set_beta
  public :: cluster_get_beta
  public :: cluster_set_minwidth
  public :: cluster_set_min_eff
  public :: cluster

  interface cluster
     module procedure cls_3d_mf
  end interface

contains

  subroutine cluster_set_minwidth(val)
    integer, intent(in) :: val
    minwidth = val
  end subroutine cluster_set_minwidth

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

  subroutine filter_lor(out, in)
    logical, intent(inout) :: out(:,:,:,:)
    logical, intent(in)    ::  in(:,:,:,:)
    out = out .or. in
  end subroutine filter_lor

  subroutine cls_3d_mf(boxes, tagboxes, buf_wid, overall_eff, blocking_factor)
    use bl_error_module
    use bl_prof_module
    type(boxarray), intent(out) :: boxes
    type(lmultifab), intent(in) :: tagboxes
    integer, intent(in) :: buf_wid
    real(dp_t), intent(out), optional :: overall_eff
    integer, intent(in), optional :: blocking_factor

    type(layout) :: la, lag, la_buf
    type(list_box) :: lboxes
    type(lmultifab) :: buf, lbuf, bufg
    type(boxarray) :: bag
    logical, pointer :: mt(:,:,:,:), mb(:,:,:,:)
    integer :: b, dm, num_flag, i, k, bxcnt, lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
    integer, allocatable :: iprocs(:), bxs(:)
    integer(kind=ll_t) :: bboxinte
    type(box) :: bx
    real(dp_t) :: bx_eff
    type(list_box_node), pointer :: bln
    type(bl_prof_timer), save :: bpt

    call build(bpt, "cluster")

    if ( present(blocking_factor) ) &
       call bl_error("CLUSTER: blocking factor not implemented yet")

    if ( nghost(tagboxes) /= 0 ) &
       call bl_error("CLUSTER: tagboxes must have NG = 0")

    if ( minwidth < 1 ) &
       call bl_error("CLUSTER: minwidth  must be > 0; ", minwidth)

    if ( min_eff < 0 .or. min_eff > 1.0_dp_t ) &
       call bl_error("CLUSTER: min_eff must be >= 0 and <= 1: ", min_eff)

    if ( buf_wid < 0 ) &
       call bl_error("CLUSTER: buf_wid must be >= 0: ", buf_wid)

    dm       = tagboxes%dim
    num_flag = lmultifab_count(tagboxes)

    if ( num_flag == 0 ) then
!      call bl_warn("No cells are flagged/No boxes returned")
       if ( present(overall_eff) ) overall_eff = 0
       return
    end if

    la_buf = get_layout(tagboxes)

    call build(buf, la_buf, nc = 1, ng = buf_wid)
    call setval(buf, .false., all = .true.)
    !
    ! Buffer the cells.
    !
    do b = 1, buf%nboxes; if ( remote(buf, b) ) cycle
       mt => dataptr(tagboxes, b, get_box(tagboxes,b))
       mb => dataptr(buf, b)
       select case (tagboxes%dim)
       case (1)
          call buffer_1d(mb(:,1,1,1), mt(:,1,1,:), buf_wid)
       case (2)
          call buffer_2d(mb(:,:,1,1), mt(:,:,1,:), buf_wid)
       case (3)
          call buffer_3d(mb(:,:,:,1), mt(:,:,:,:), buf_wid)
       end select
    end do
    !
    ! Set any valid region that can be covered by a tagged periodically-shifted ghost cell.
    !
    call map_periodic(buf)
    !
    ! Remove any tags outside the problem domain.
    !
    call pd_mask(buf)
    !
    ! Make sure that all tagged cells in are properly replicated in the high indexed boxes.
    !
    call internal_sync(buf, all = .true., filter = filter_lor)

    bx = get_pd(la_buf)

    if ( parallel_nprocs() > 1 ) then
       !
       ! The cluster algorithm is inherently serial.
       ! We'll set up the problem to do on the IO processor.
       !
       call boxarray_build_copy(bag, get_boxarray(buf))
       call boxarray_grow(bag, buf_wid)
       call build(lag, bag, bx, get_pmask(la_buf), explicit_mapping = get_proc(la_buf))
       call destroy(bag)
       call build(bufg, lag, nc = 1, ng = 0)

       do i = 1, bufg%nboxes
          if ( remote(bufg, i) ) cycle
          mt => dataptr(buf,  i)
          mb => dataptr(bufg, i)
          mb = mt
       end do

       call destroy(buf)

       allocate(iprocs(nboxes(get_layout(buf))))

       iprocs = parallel_IOProcessorNode()

       call build(la, get_boxarray(buf%la), bx, get_pmask(buf%la), explicit_mapping = iprocs)

       call build(lbuf, la, nc = 1)

       call copy(lbuf, bufg)  ! This is a parallel copy.

       call destroy(bufg)
       call destroy(lag)

       if ( parallel_IOProcessor() ) then
          call owner_mask(lbuf)
          call cluster_mf_private_recursive(lboxes, bx, bx_eff, lbuf)
          bxcnt = size(lboxes)
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

    else
       call owner_mask(buf)
       call cluster_mf_private_recursive(lboxes, bx, bx_eff, buf)
       call destroy(buf)
    end if

    call build(boxes, lboxes)

    call destroy(lboxes)
    !
    ! TODO - remove this check after it's clear the parallel implementation is OK.
    !
    if ( .not. boxarray_clean(boxes%bxs) ) then
       call bl_error('cls_3d_mf: boxes are NOT disjoint')
    end if

    if ( present(overall_eff) ) then
       bboxinte    = volume(boxes)
       overall_eff = real(count(buf,all=.true.),dp_t) / real(bboxinte, dp_t)
    end if

    call destroy(bpt)

  contains

    subroutine pd_mask(mask)
      type(lmultifab), intent(inout) :: mask
      integer :: i, j
      type(box) :: bxi, pd
      type(boxarray) :: ba

      pd = get_pd(get_layout(mask))
      do i = 1, mask%nboxes; if ( remote(mask, i) ) cycle
         bxi = get_pbox(mask, i)
         call boxarray_box_diff(ba, bxi, pd)
         do j = 1, nboxes(ba)
            call setval(mask%fbs(i), .false., get_box(ba,j))
         end do
         call destroy(ba)
      end do
    end subroutine pd_mask

    subroutine owner_mask(mask)
      type(lmultifab), intent(inout) :: mask
      integer :: i, j
      type(box) :: bxi, bxj, bxij

      do i = 1, mask%nboxes; if ( remote(mask, i) ) cycle
         bxi = get_pbox(mask, i)
         do j = 1, i-1
            bxj = get_pbox(mask, j)
            bxij = intersection(bxi, bxj)
            if ( empty(bxij) ) cycle
            call setval(mask%fbs(j), .false., bxij)
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

    integer                        :: i, j, ii, jj, cnt, shft(3**mask%dim,mask%dim)
    type(box)                      :: pd, bxs(3**mask%dim), bx_from, bx_to
    logical                        :: pmask(mask%dim)
    logical, pointer               :: ap(:,:,:,:), bp(:,:,:,:)
    logical, allocatable           :: pt(:,:,:,:)
    integer, parameter             :: tag = 2121
    type(box_intersector), pointer :: bi(:)
    type(bl_prof_timer),   save    :: bpt

    call build(bpt, "map_periodic")

    pd    = get_pd(mask%la)
    pmask = layout_get_pmask(mask%la)

    if ( .not. any(pmask) ) return

    do i = 1, nboxes(mask)
       call box_periodic_shift(pd, get_ibox(mask,i), mask%nodal, pmask, mask%ng, shft, cnt, bxs)

       do jj = 1, cnt
          bi => layout_get_box_intersector(mask%la, bxs(jj))

          do ii = 1, size(bi)
             j = bi(ii)%i

             if ( remote(mask,i) .and. remote(mask,j) ) cycle

             bx_to   = bi(ii)%bx
             bx_from = shift(bx_to, -shft(jj,:))

             if ( local(mask,i) .and. local(mask,j) ) then
                ap => dataptr(mask,j,bx_to)
                bp => dataptr(mask,i,bx_from)
                ap = ap .or. bp
             else if ( local(mask,i) ) then
                !
                ! We own index i.  Got to send it to processor owning index j.
                !
                bp => dataptr(mask,i,bx_from)
                call parallel_send(bp, get_proc(mask%la,j), tag)
             else
                !
                ! We own index j.  Got to get index i data from processor owning it.
                !
                ap => dataptr(mask,j,bx_to)
                allocate(pt(1:size(ap,1),1:size(ap,2),1:size(ap,3),1))
                call parallel_recv(pt, get_proc(mask%la,i), tag)
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
    type(box) :: bx1
    integer :: llo(3), hho(3)

    dm = tagboxes%dim

    tx = 0
    ty = 0
    tz = 0

    do n = 1, tagboxes%nboxes; if ( remote(tagboxes, n) ) cycle
       bx1 = intersection(get_pbox(tagboxes, n), bx)
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
      do while ( sig(llo) == 0  .and.  hho-llo >= minwidth )
         llo = llo + 1
      end do

      do while ( sig(hho) ==  0  .and. hho-llo >= minwidth )
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

    dm  = tagboxes%dim

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
        type(box) :: bx1
        r = 0
        do n = 1, tagboxes%nboxes;
           bx1 = intersection(get_pbox(tagboxes, n), bx)
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
    if ( extent(bx,i) >= 2*minwidth ) then
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
      do i = ll + minwidth, hh - minwidth
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
      do i = ll + minwidth, hh - max(minwidth,1)
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

end module cluster_module
