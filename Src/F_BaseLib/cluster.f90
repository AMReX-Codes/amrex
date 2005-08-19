module cluster_module

  use bl_types
  use bl_error_module
  use box_module
  use list_box_module
  use boxarray_module
  use multifab_module

  implicit none

  private

  logical, private, save :: verbose = .false.
  integer, private, save :: cut_threshold = 2

  public :: cluster_set_verbose
  public :: cluster_set_cut_threshold
  public :: cluster

  interface cluster
     module procedure cls_3d_mf
  end interface

contains

  subroutine cluster_set_verbose(val)
    logical, intent(in) :: val
    verbose = val
  end subroutine cluster_set_verbose

  subroutine cluster_set_cut_threshold(val)
    integer, intent(in) :: val
    cut_threshold= val
  end subroutine cluster_set_cut_threshold

  subroutine filter_lor(out, in)
    logical, intent(inout) :: out(:,:,:,:)
    logical, intent(in)    ::  in(:,:,:,:)
    out = out .or. in
  end subroutine filter_lor

  subroutine cls_3d_mf(boxes, tagboxes, minwidth, buf_wid, min_eff, overall_eff, blocking_factor)
    type(boxarray), intent(out) :: boxes
    type(lmultifab), intent(in) :: tagboxes
    integer, intent(in) :: minwidth
    integer, intent(in) :: buf_wid
    real(dp_t), intent(in) :: min_eff
    real(dp_t), intent(out), optional :: overall_eff
    integer, intent(in), optional :: blocking_factor

    type(list_box) :: lboxes
    type(lmultifab) :: buf
    logical, pointer :: mt(:,:,:,:)
    logical, pointer :: mb(:,:,:,:)
    integer :: b
    integer :: num_flag
    integer :: dm, i
    integer :: bboxinte
    integer :: lblocking_factor
    type(box) :: bx

    lblocking_factor = 1
    if ( present(blocking_factor) ) then
       call bl_error("CLUSTER: blocking factor not implemented yet")
    end if

    if ( nghost(tagboxes) /= 0 ) then
       call bl_error("CLUSTER: tagboxes must have NG = 0")
    end if

    dm = tagboxes%dim

    call build(lboxes)
    call build(buf, get_layout(tagboxes), nc = 1, ng = buf_wid)

    num_flag = lmultifab_count(tagboxes)

    if ( num_flag == 0 ) then
       call bl_warn("No cells are flagged/No boxes returned")
       if ( present(overall_eff) ) then
          overall_eff = 0
       end if
       return
    end if

    ! buffer the cells
    do b = 1, buf%nboxes; if ( remote(buf, b) ) cycle
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

print *, 'HERE'
print *, 'count1 ', count(buf, all = .true.)
    ! remove any tags outside the problem domain.
    call pd_mask(buf)
print *, 'count2 ', count(buf, all = .true.)
    ! make sure that all tagged cells in are replicated in the high
    ! indexed boxes
    call internal_sync(buf, all = .true., filter = filter_lor)
print *, 'count3 ', count(buf, all = .true.)
    ! remove all tags from the low index fabs that overlap high index fabs
    call owner_mask(buf)
print *, 'count4 ', count(buf, all = .true.)

    bx = get_pd(get_layout(buf))
    call cluster_mf_private_recursive(lboxes, bx, buf, minwidth, min_eff)

    call build(boxes, lboxes)
    boxes%dim = dm

    call destroy(lboxes)

    if ( present(overall_eff) ) then
       bboxinte = volume(boxes)
       overall_eff = real(count(buf,all=.true.),dp_t) / real(bboxinte, dp_t)
       do i = 1, nboxes(boxes)
          print *, i, box_eff_mf(buf, get_box(boxes,i))
       end do
    end if

    call destroy(buf)

  contains

    subroutine pd_mask(mask)
      type(lmultifab), intent(inout) :: mask
      integer :: i, j
      type(box) :: bxi, bxj, bxij, pd
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
      integer :: i, j, k, l, m, n

      do i = 1, size(tt,1)
         if ( any(tt(i,:)) ) bb(i-ng:i+ng) = .true.
      end do

    end subroutine buffer_1d

    subroutine buffer_2d(bb, tt, ng)
      integer, intent(in) :: ng
      logical, intent(in) :: tt(:,:,:)
      logical, intent(out) :: bb(1-ng:,1-ng:)
      integer :: i, j, k, l, m, n

      do j = 1, size(tt,2); do i = 1, size(tt,1)
         if ( any(tt(i,j,:)) ) bb(i-ng:i+ng,j-ng:j+ng) = .true.
      end do; end do

    end subroutine buffer_2d

    subroutine buffer_3d(bb, tt, ng)
      integer, intent(in) :: ng
      logical, intent(in) :: tt(:,:,:,:)
      logical, intent(out) :: bb(1-ng:,1-ng:,1-ng:)
      integer :: i, j, k, l, m, n

      do k = 1, size(tt,3); do j = 1, size(tt,2); do i = 1, size(tt,1)
         if ( any(tt(i,j,k,:)) ) bb(i-ng:i+ng,j-ng:j+ng,k-ng:k+ng) = .true.
      end do; end do; end do

    end subroutine buffer_3d

  end subroutine cls_3d_mf

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
      do while ( sig(llo) == 0  .and.  llo < hh )
         llo = llo + 1
      end do

      hho = hh
      do while ( sig(hho) ==  0  .and.  hho > llo )
         hho = hho - 1
      end do

    end subroutine bracket

  end function cluster_tally

  recursive subroutine cluster_mf_private_recursive(boxes, bx, tagboxes, minwidth, min_eff)
    type(list_box), intent(out) :: boxes
    type(box), intent(in) :: bx
    type(lmultifab), intent(in) ::  tagboxes
    integer, intent(in) ::  minwidth
    real(dp_t), intent(in) ::   min_eff
    real(dp_t) bx_eff
    type(box) :: bbx
    integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
    type(list_box) :: C1, C2
    type(box) :: b1, b2
    integer, allocatable :: sx(:), sy(:), sz(:)
    integer :: ll(3), hh(3)
    integer :: dm

    dm  = tagboxes%dim

call print(bx, 'in bx')
    if ( empty(bx) ) then
       return
    endif

    ll = 1; ll(1:dm) = lwb(bx)
    hh = 1; hh(1:dm) = upb(bx)

    allocate(sx(ll(1):hh(1)), sy(ll(2):hh(2)), sz(ll(3):hh(3)))

    bbx = cluster_tally(tagboxes, bx, sx, sy, sz, ll, hh)

call print(bbx, '  bbx')

    if ( empty(bbx) ) then
       return
    endif

    bx_eff = box_eff_mf(tagboxes, bbx)

print *, 'bx_eff ', bx_eff
print *, 'min_eff ', min_eff

    if ( bx_eff >= min_eff ) then
       call push_back(boxes, bbx)
    else
       lo = 1; lo(1:dm) = lwb(bbx); hi = 1; hi(1:dm) = upb(bbx)
       call find_split(bbx, b1, b2, minwidth, &
            sx(lo(1):hi(1)), sy(lo(2):hi(2)), sz(lo(3):hi(3)), lo, hi)
call print(b1, 'b1')
call print(b2, 'b2')
       call cluster_mf_private_recursive(c1, b1, tagboxes, minwidth, min_eff)
       call cluster_mf_private_recursive(c2, b2, tagboxes, minwidth, min_eff)
print *, 'size(c1) = ', size(c1)
print *, 'size(c2) = ', size(c2)
       if ( size(c1) > 0 .or. size(c2) > 0 ) then
          boxes = c1
          call splice(boxes, c2)
       else
          call push_back(boxes, bbx)
       end if
    end if

  end subroutine cluster_mf_private_recursive

  function box_eff_mf(tagboxes, bx) result(r)
    real(dp_t) :: r, r1
    type(box), intent(in) :: bx
    type(lmultifab), intent(in) :: tagboxes
    logical, pointer :: tp(:,:,:,:)
    integer :: n
    type(box) :: bx1

    r1 = 0
    do n = 1, tagboxes%nboxes; if ( remote(tagboxes, n) ) cycle
       bx1 = intersection(get_pbox(tagboxes, n), bx)
       if ( empty(bx1) ) cycle
       tp => dataptr(tagboxes, n, bx1)
       r1 = r1 + real(count(tp),dp_t)
    end do
    call parallel_reduce(r,r1, MPI_SUM)
    r = r/dvolume(bx)
  end function box_eff_mf

  subroutine find_split(bx, b1, b2, minwidth, sx, sy, sz, ll, hh)
    type(box), intent(in) ::  bx
    type(box), intent(out) :: b1, b2
    integer, intent(in) :: minwidth
    integer, intent(in) :: ll(:), hh(:)
    integer, intent(in) :: sx(ll(1):), sy(ll(2):), sz(ll(3):)
    integer, allocatable :: lx(:), ly(:), lz(:)
    integer :: dm, n, i, j, k, ii, jj, kk
    type(box) :: bx1
    integer :: hi(3), ip(3)

    dm = bx%dim

    if ( holes(bx, b1, b2, sx, ll(1), hh(1), minwidth, 1) ) return
    if ( holes(bx, b1, b2, sy, ll(2), hh(2), minwidth, 2) ) return
    if ( holes(bx, b1, b2, sz, ll(3), hh(3), minwidth, 3) ) return

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
    call inflection(lx, ll(1), hh(1), minwidth, hi(1), ip(1))
    call inflection(ly, ll(2), hh(2), minwidth, hi(2), ip(2))
    call inflection(lz, ll(3), hh(3), minwidth, hi(3), ip(3))

    if ( any(ip /= -1) ) then
       if ( maxval(hi) >= CUT_THRESHOLD ) then
          i = maxloc(hi, dim=1)
          call box_chop(bx, b1, b2, i, ip(i))
print *, 'extents b1,b2 ', extent(b1), extent(b2)
          return
       end if
    end if

    ! final fall back:

print *, 'fall back'
print *, 'extent ', extent(bx)
    i = maxloc(extent(bx),dim=1)
call print(bx, 'bx')
print *, 'i = ', i
    if ( extent(bx,i) >= 2*minwidth ) then
       call box_chop(bx, b1, b2, i, lwb(bx,i) + extent(bx,i)/2)
call print(b1, 'b1')
call print(b2, 'b2')
print *, 'extents b1,b2 ', extent(b1), extent(b2)
    end if

  contains

    function holes(bx, b1, b2, sig, ll, hh, minwidth, dim) result(r)
      logical :: r
      integer, intent(in) :: minwidth
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
print *, 'extents b1,b2 ', extent(b1), extent(b2)
            return
         end if
      end do

    end function holes

    subroutine inflection(lp, ll, hh, minwidth, hiv, infp)
      logical :: r
      integer, intent(in) :: minwidth
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
