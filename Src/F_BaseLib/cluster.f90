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
     module procedure cls_3d
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

    ! remove any tags outside the problem domain.
    call pd_mask(buf)
    ! make sure that all tagged cells in are replicated in the high
    ! indexed boxes
    call internal_sync(buf, all = .true., filter = filter_lor)
    ! remove all tags from the low index fabs that overlap high index fabs
    call owner_mask(buf)

    call cluster_mf_private(lboxes, buf, minwidth, min_eff)

    call build(boxes, lboxes)
    do i = 1, nboxes(boxes)
       if ( dm < 3 ) boxes%bxs(i) = reduce(get_box(boxes,i),dim=3)
       if ( dm < 2 ) boxes%bxs(i) = reduce(get_box(boxes,i),dim=2)
    end do
    boxes%dim = dm

    call destroy(lboxes)

    if ( present(overall_eff) ) then
       bboxinte = volume(boxes)
       overall_eff = real(count(buf,all=.true.),dp_t) / real(bboxinte, dp_t)
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

  subroutine cls_3d(boxes, tagbox, minwidth, buf_wid, min_eff, overall_eff, dim)
    type(boxarray), intent(out) :: boxes
    logical, intent(in) :: tagbox(:,:,:)
    integer, intent(in) :: minwidth
    integer, intent(in) :: buf_wid
    real(dp_t), intent(in) :: min_eff
    real(dp_t), intent(out), optional :: overall_eff
    integer, optional, intent(in) :: dim

    type(list_box_node), pointer :: bn

    logical, allocatable :: buf_box(:,:,:)
    integer :: nx, ny, nz
    type(list_box) :: lboxes

    integer :: i, j, k, l, m, n 
    integer :: num_flag
    integer :: bboxinte

    integer :: dm

    call build(lboxes)

    dm = 3; if ( present(dim) ) dm = dim
    nx = size(tagbox,1); ny = size(tagbox,2); nz = size(tagbox,3)

    allocate(buf_box(nx,ny,nz))

    num_flag = count(tagbox)

    if ( num_flag  ==  0 ) then
       call bl_warn('No cells are flagged/No boxes returned')
       if ( present(overall_eff) ) then
          overall_eff = 0
       end if
       return
    end if

    buf_box = .False.

    do i = 1, nx; do j = 1, ny; do k = 1, nz
       if( tagbox(i,j,k) ) then
          do l = max(1,i-buf_wid), min(i+buf_wid,nx)
             do m = max(1,j-buf_wid), min(j+buf_wid,ny)
                do n = max(1,k-buf_wid), min(k+buf_wid,nz)
                   buf_box(l,m,n) = .True.
                end do
             end do
          end do
       end if
    end do; end do; end do

    call cluster_l_private(lboxes, buf_box, minwidth, min_eff)

    if ( present(overall_eff) ) then
       bboxinte = 0
       bn => begin(lboxes)
       do while ( associated(bn) )
          bboxinte = bboxinte + volume(value(bn))
          bn => next(bn)
       end do
       overall_eff = real(count(buf_box),dp_t) / real(bboxinte, dp_t)
    end if

    call build(boxes, lboxes)
    do i = 1, nboxes(boxes)
       if ( dm < 3 ) boxes%bxs(i) = reduce(get_box(boxes,i),dim=3)
       if ( dm < 2 ) boxes%bxs(i) = reduce(get_box(boxes,i),dim=2)
    end do
    boxes%dim = dm

    call destroy(lboxes)

  end subroutine cls_3d

  subroutine cluster_mf_private(boxes, tagboxes, minwidth, min_eff)
    type(list_box), intent(inout) :: boxes
    type(lmultifab), intent(in) ::  tagboxes
    integer, intent(in) ::  minwidth
    real(dp_t), intent(in) ::   min_eff
    real(dp_t) cur_box_eff
    type(list_box_node), pointer :: bn
    integer, allocatable :: sigx(:), sigy(:), sigz(:)
    integer, allocatable :: lplx(:), lply(:), lplz(:)
    logical :: flag
    type(box) :: bbx
    integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)


    lo = 0 ; hi = 0
    bbx = bbox(get_boxarray(tagboxes))
    lo(1:bbx%dim) = lwb(bbx); hi(1:bbx%dim) = upb(bbx)
    bbx = make_box(lo, hi)
    allocate(sigx(lo(1):hi(1)), sigy(lo(2):hi(2)), sigz(lo(3):hi(3)))
    allocate(lplx(lo(1):hi(1)), lply(lo(2):hi(2)), lplz(lo(3):hi(3)))
    call push_back(boxes, bbx)
    bn => begin(boxes)
    do while ( associated(bn) )
       cur_box_eff = box_eff_mf(tagboxes, value(bn))
       if ( verbose  .and. parallel_IOProcessor() ) then
          print *, 'cur_box_eff', cur_box_eff
       end if
       if ( cur_box_eff  <  min_eff ) then
          call sigma_laplace_mf(tagboxes, value(bn), sigx, sigy, sigz, lplx, lply, lplz, lo)
          flag = find_split(boxes, bn, minwidth, sigx, sigy, sigz, lplx, lply, lplz)
          if ( .not. flag ) then
             bn => next(bn)
          end if
       else
          bn => next(bn)
       end if
    end do

  end subroutine cluster_mf_private

  subroutine cluster_l_private(boxes, tagbox, minwidth, min_eff)
    type(list_box), intent(inout) :: boxes
    logical, intent(in) ::  tagbox(0:,0:,0:)
    integer, intent(in) ::  minwidth
    real(dp_t), intent(in) ::   min_eff

    type(list_box_node), pointer :: bn
    integer :: sigx(size(tagbox,1)), sigy(size(tagbox,2)), sigz(size(tagbox,3))
    integer :: lplx(size(tagbox,1)), lply(size(tagbox,2)), lplz(size(tagbox,3))
    logical :: flag
    integer :: nx, ny, nz
    real(dp_t) :: cur_box_eff

    nx = size(tagbox,1); ny = size(tagbox,2); nz = size(tagbox,3)

    call push_back(boxes, make_box((/0, 0, 0/), (/nx-1, ny-1, nz-1/)))
    bn => begin(boxes)
    do while ( associated(bn) )
       cur_box_eff = box_eff(tagbox, value(bn))
       if ( verbose  .and. parallel_IOProcessor() ) then
          print *, 'cur_box_eff', cur_box_eff
       end if
       if ( cur_box_eff  <  min_eff ) then
          call sigma_laplace(tagbox, value(bn), sigx, sigy, sigz, lplx, lply, lplz)
          flag = find_split(boxes, bn, minwidth, sigx, sigy, sigz, lplx, lply, lplz)
          if ( .not. flag ) then
             bn => next(bn)
          end if
       else
          bn => next(bn)
       end if
    end do

  end subroutine cluster_l_private

  subroutine sigma_laplace_mf(tagboxes, bx, sigx, sigy, sigz, lplx, lply, lplz, lo)
    type(lmultifab), intent(in) :: tagboxes
    integer, intent(in):: lo(:)
    type(box), intent(in) :: bx
    integer, intent(out) :: sigx(lo(1):),sigy(lo(2):),sigz(lo(3):)
    integer, intent(out) :: lplx(lo(1):),lply(lo(2):),lplz(lo(3):)
    logical, pointer :: tp(:,:,:,:)
    !   logical, pointer :: mp(:,:,:,:)
    integer :: n
    type(box) :: bx1

    integer, allocatable :: tsigx(:),tsigy(:),tsigz(:)
    integer :: i, j, k, lx, ly, lz, hx, hy, hz
    integer :: hh(4), dm, ll(3), ii, jj, kk

    dm = tagboxes%dim

    sigx = 0
    sigy = 0
    sigz = 0

    lx = lwb(bx,1); ly = 1; lz = 1
    hx = upb(bx,1); hy = 1; hz = 1
    if ( tagboxes%dim > 1 ) then
       ly = lwb(bx,2)
       hy = upb(bx,2)
    else if ( tagboxes%dim > 2 ) then
       lz = lwb(bx,3)
       hz = upb(bx,3)
    end if

    allocate(tsigx(lx:hx), tsigy(ly:hy), tsigz(lz:hz))
    tsigx = 0
    tsigy = 0
    tsigz = 0

    do n = 1, tagboxes%nboxes; if ( remote(tagboxes, n) ) cycle
       bx1 = intersection(get_pbox(tagboxes, n), bx)
       if ( empty(bx1) ) cycle
       tp => dataptr(tagboxes, n, bx1)

       ll = 1; ll(1:dm) = lwb(bx1)
       hh = 1; hh(1:dm) = extent(bx1)
       do i = 1, hh(1)
          ii = ll(1) + i - 1
          tsigx(ii) = tsigx(ii) + count(tp(i,:,:,1))
       end do

       do j = 1, hh(2)
          jj = ll(2) + j - 1
          tsigy(jj) = tsigy(jj) + count(tp(:,j,:,1))
       end do

       do k = 1, hh(3)
          kk = ll(3) + k - 1
          tsigz(kk) = tsigz(kk) + count(tp(:,:,k,1))
       end do
    end do

    call parallel_reduce(sigx(lx:hx), tsigx, MPI_SUM)
    call parallel_reduce(sigy(ly:hy), tsigy, MPI_SUM)
    call parallel_reduce(sigz(lz:hz), tsigz, MPI_SUM)

    !! Note: only one of berger/rigotsis schemes is here used.
    !! Note: a fill boundary needs to have b???
    lplx = 0
    lply = 0
    lplz = 0

    do i = lx+1, hx-1
       lplx(i) = sigx(i+1)-2*sigx(i)+sigx(i-1)
    end do

    do j = ly+1, hy-1
       lply(j) = sigy(j+1)-2*sigy(j)+sigy(j-1)
    end do

    do k = lz+1, hz-1
       lplz(k) = sigz(k+1)-2*sigz(k)+sigz(k-1)
    end do

    if ( verbose .and. parallel_IOProcessor() ) then
       print '(a,1x,20(i3,1x))', 'sigx', sigx
       print '(a,1x,20(i3,1x))', 'sigy', sigy
       print '(a,1x,20(i3,1x))', 'sigz', sigz
       print '(a,1x,20(i3,1x))', 'lplx', lplx
       print '(a,1x,20(i3,1x))', 'lply', lply
       print '(a,1x,20(i3,1x))', 'lplz', lplz
    end if

  end subroutine sigma_laplace_mf

  subroutine sigma_laplace(tagbox, bx, sigx, sigy, sigz, lplx, lply, lplz)
    logical, intent(in) :: tagbox(0:,0:,0:)
    type(box), intent(in) :: bx
    integer, intent(out) :: sigx(0:),sigy(0:),sigz(0:)
    integer, intent(out) :: lplx(0:),lply(0:),lplz(0:)

    integer :: i, j, k, lx, ly, lz, hx, hy, hz

    sigx = 0
    sigy = 0
    sigz = 0

    lx = lwb(bx,1); ly = lwb(bx,2); lz = lwb(bx,3)
    hx = upb(bx,1); hy = upb(bx,2); hz = upb(bx,3)

    do i = lx, hx
       sigx(i) = count(tagbox(i,ly:hy,lz:hz))
    end do

    do j = ly, hy
       sigy(j) = count(tagbox(lx:hx,j,lz:hz))
    end do

    do k = lz, hz
       sigz(k) = count(tagbox(lx:hx,ly:hy,k))
    end do

    lplx = 0
    lply = 0
    lplz = 0

    if ( .true. ) then

       do i = lx+1, hx-1
          lplx(i) = sigx(i+1)-2*sigx(i)+sigx(i-1)
       end do

       do j = ly+1, hy-1
          lply(j) = sigy(j+1)-2*sigy(j)+sigy(j-1)
       end do

       do k = lz+1, hz-1
          lplz(k) = sigz(k+1)-2*sigz(k)+sigz(k-1)
       end do

    else 

       do i = lx, hx-1
          lplx(i) = sum(absdiff(tagbox(i,ly:hy,lz:hz),tagbox(i+1,ly:hy,lz:hz)))
       end do

       do j = ly, hy-1
          lply(j) = sum(absdiff(tagbox(lx:hx,j,lz:hz),tagbox(lx:hx,j+1,lz:hz)))
       end do
       do k = lz, hz-1
          lplz(k) = sum(absdiff(tagbox(lx:hx,ly:hy,k),tagbox(lx:hx,ly:hy,k+1)))
       end do

       do i = hx-1, lx+1, -1
          lplx(i) = lplx(i)-lplx(i-1)
       end do

       do j = hy-1, ly+1, -1
          lply(j) = lply(j)-lply(j-1)
       end do

       do k = hz-1, lz+1, -1
          lplz(k) = lplz(k)-lplz(k-1)
       end do

       lplx(lx) = 0; lply(ly) = 0; lplz(lz) = 0
       lplx(hx) = 0; lply(hy) = 0; lplz(hz) = 0


    end if

    if ( verbose  .and. parallel_IOProcessor() ) then
       print '(a,1x,20(i3,1x))', 'sigx', sigx
       print '(a,1x,20(i3,1x))', 'sigy', sigy
       print '(a,1x,20(i3,1x))', 'sigz', sigz
       print '(a,1x,20(i3,1x))', 'lplx', lplx
       print '(a,1x,20(i3,1x))', 'lply', lply
       print '(a,1x,20(i3,1x))', 'lplz', lplz
    end if
  contains
    elemental function absdiff(L1, L2) result(i)
      logical, intent(in) :: l1, l2
      integer :: i, i1, i2
      i1 = 0; if ( l1 ) i1 = 1
      i2 = 0; if ( l2 ) i2 = 1
      i = abs(i1-i2)
    end function absdiff

  end subroutine sigma_laplace

  function box_eff_mf(tagboxes, bx) result(r)
    real(dp_t) :: r
    type(box), intent(in) :: bx
    type(lmultifab), intent(in) :: tagboxes
    logical, pointer :: tp(:,:,:,:)
    integer :: n
    type(box) :: bx1

    r = 0
    do n = 1, tagboxes%nboxes; if ( remote(tagboxes, n) ) cycle
       bx1 = intersection(get_pbox(tagboxes, n), bx)
       if ( empty(bx1) ) cycle
       tp => dataptr(tagboxes, n, bx1)
       r = r + real(count(tp),dp_t)
    end do
    r = r/dvolume(bx)
  end function box_eff_mf

  function box_eff(tagbox, bx) result(r)
    real(dp_t) :: r

    type(box), intent(in) :: bx
    logical, intent(in) :: tagbox(0:,0:,0:)

    integer :: lx, ly, lz, hx, hy, hz

    lx = lwb(bx,1); ly = lwb(bx,2); lz = lwb(bx,3)
    hx = upb(bx,1); hy = upb(bx,2); hz = upb(bx,3)

    r = real(count(tagbox(lx:hx, ly:hy, lz:hz)),dp_t)/dvolume(bx)

  end function box_eff

  function find_split(boxes, bn, minwidth, sigx, sigy, sigz, lplx, lply, lplz ) result(r)
    logical :: r, rr
    type(list_box), intent(inout) ::  boxes
    type(list_box_node), pointer :: bn
    integer, intent(in) :: minwidth
    integer, intent(in) :: sigx(0:), sigy(0:), sigz(0:)
    integer, intent(in) :: lplx(0:), lply(0:), lplz(0:)


    rr = find_holes(boxes, bn, minwidth, sigx, sigy, sigz)
    if ( verbose  .and. parallel_IOProcessor() ) then
       print *, 'FIND_SPLIT(1) r = ', rr
    end if
    if ( .not. rr ) then
       rr = find_inflx(boxes, bn, minwidth, lplx, lply, lplz)
    end if
    if ( verbose  .and. parallel_IOProcessor() ) then
       print *, 'FIND_SPLIT(2) r = ', rr
    end if
    r = rr

  end function find_split

  function find_split0(boxes, bn, minwidth, sigx, sigy, sigz, lplx, lply, lplz ) result(r)
    logical :: r
    type(list_box), intent(inout) ::  boxes
    type(list_box_node), pointer :: bn
    integer, intent(in) :: minwidth
    integer, intent(in) :: sigx(0:), sigy(0:), sigz(0:)
    integer, intent(in) :: lplx(0:), lply(0:), lplz(0:)


    r = find_holes(boxes, bn, minwidth, sigx, sigy, sigz)
    if ( verbose  .and. parallel_IOProcessor() ) then
       print *, 'FIND_SPLIT(1) r = ', r
    end if
    if ( .not. r ) then
       r = find_inflx(boxes, bn, minwidth, lplx, lply, lplz)
    end if
    if ( verbose  .and. parallel_IOProcessor() ) then
       print *, 'FIND_SPLIT(2) r = ', r
    end if

  end function find_split0

  function find_holes(boxes, bn, minwidth, sigx, sigy, sigz) result(r)
    logical :: r
    type(list_box), intent(inout) ::  boxes
    type(list_box_node), pointer :: bn
    integer, intent(in) :: minwidth
    integer, intent(in) :: sigx(0:), sigy(0:), sigz(0:)

    integer :: lx, hx, ly, hy, lz, hz

    if ( verbose  .and. parallel_IOProcessor() ) then
       call print(value(bn), 'find_holes')
    end if

    r  = .false.

    lx = lwb(value(bn),1); ly = lwb(value(bn),2); lz = lwb(value(bn),3)
    hx = upb(value(bn),1); hy = upb(value(bn),2); hz = upb(value(bn),3)

    call bracket(sigx, lx, hx)
    call bracket(sigy, ly, hy)
    call bracket(sigz, lz, hz)

    call set(bn, make_box((/lx, ly, lz/), (/hx, hy, hz/)))

    if ( verbose  .and. parallel_IOProcessor() ) then
       call print(value(bn), 'find_holes ==>')
    end if

    if ( holes(sigx, lx, hx, 1) ) then
       r = .true.
    else if ( holes(sigy, ly, hy, 2) ) then
       r = .true.
    else if( holes(sigz, lz, hz, 3) ) then
       r = .true.
    end if

    if ( verbose  .and. parallel_IOProcessor() ) then
       print *, 'FIND_HOLES => ', r
    end if

  contains

    function holes(sig, ll, hh, dim) result(r)
      logical :: r
      integer, intent(in) :: ll, hh, dim
      integer, intent(in)    :: sig(0:)
      type(box) :: b1, b2
      integer :: i
      r  = .false.
      do i = ll+minwidth, hh-minwidth
         if ( sig(i)  ==  0 ) then
            call box_chop(value(bn), b1, b2, dim, i)
            b2 = grow(b2, -1, dim, -1)
            call boxStack(boxes, bn, b1, b2)
            r  = .true.
            return
         end if
      end do
    end function holes

    subroutine bracket(sig, ll, hh)
      integer, intent(inout) :: ll, hh
      integer, intent(in)    :: sig(0:)

      integer :: i

      i = ll
      do while ( sig(i)  ==  0  .and.  (hh-ll)  >=  minwidth)
         ll = i + 1
         r = .true.
         i = i + 1
      end do

      i = hh
      do while ( sig(i)  ==  0  .and.  (hh-ll)  >=  minwidth)
         hh = i - 1
         r = .true.
         i = i - 1
      end do

    end subroutine bracket

  end function find_holes

  subroutine boxStack(boxes, bn, bx1, bx2)
    type(list_box), intent(inout) :: boxes
    type(list_box_node), pointer :: bn
    type(box), intent(in) :: bx1, bx2

    call set(bn, bx1)
    call push_back(boxes, bx2)

  end subroutine boxStack

  function find_inflx(boxes, bn, minwidth, lplx, lply, lplz) result(r)
    logical :: r
    type(list_box), intent(inout) ::  boxes
    type(list_box_node), pointer :: bn
    integer, intent(in) :: minwidth
    integer, intent(in) :: lplx(0:), lply(0:), lplz(0:)

    type(box) :: b1, b2
    integer :: lx, ly, lz, hx, hy, hz, wx, wy, wz
    integer :: hix, ipx
    integer :: hiy, ipy
    integer :: hiz, ipz

    r = .false.

    lx = lwb(value(bn),1); ly = lwb(value(bn),2); lz = lwb(value(bn),3)
    hx = upb(value(bn),1); hy = upb(value(bn),2); hz = upb(value(bn),3)

    call inflection(lplx, lx, hx, hix, ipx)
    call inflection(lply, ly, hy, hiy, ipy)
    call inflection(lplz, lz, hz, hiz, ipz)

    if ( verbose  .and. parallel_IOProcessor() ) then
       call print(value(bn), 'FIND_INFLX')
       print '(a,1x,20(i3,1x))', 'lplx', lplx
       print *, '      l', lx, 'h', hx, 'hiv', hix, 'inf', ipx
       print '(a,1x,20(i3,1x))', 'lply', lply
       print *, '      l', ly, 'h', hy, 'hiv', hiy, 'inf', ipy
       print '(a,1x,20(i3,1x))', 'lplz', lplz
       print *, '      l', lz, 'h', hz, 'hiv', hiz, 'inf', ipz
    end if


    if (ipx /= -1 .or. ipy /= -1 .or. ipz /= -1) then
       if ( max(hix, hiy, hiz) >= CUT_THRESHOLD ) then
          if ( hix >= hiy .and. hix >= hiz ) then
             call box_chop(value(bn), b1, b2, 1, ipx)
          else if ( hiy >= hix .and. hiy >= hiz ) then
             call box_chop(value(bn), b1, b2, 2, ipy)
          else
             call box_chop(value(bn), b1, b2, 3, ipz)
          end if
          call boxStack(boxes, bn, b1, b2)
          r  = .true.
       end if
    else if ( .false. ) then
       wx = hx - lx + 1
       wy = hy - ly + 1
       wz = hz - lz + 1
       if ( wx >= wy .and. wx >= wz ) then
          wx = lx + wx/2
          call box_chop(value(bn), b1, b2, 1, wx)
       else if ( wy >= wx .and. wy >= wz ) then
          wy = ly + wy/2
          call box_chop(value(bn), b1, b2, 2, wy)
       else
          wz = lz + wz/2
          call box_chop(value(bn), b1, b2, 3, wz)
       end if
       call boxStack(boxes, bn, b1, b2)
       r = .true.
    end if

    if ( verbose  .and. parallel_IOProcessor() ) then
       print *, 'FIND_INFLX => ', r
    end if

  contains

    subroutine inflection(lpl, ll, hh, hiv, infp)
      integer, intent(in) :: lpl(0:), ll, hh
      integer, intent(out) :: hiv, infp
      integer :: i, tmp
      hiv  = 0
      infp = -1
      do i = ll + minwidth, hh - minwidth
         if( (lpl(i) > 0 .and. lpl(i+1) < 0)  .or. (lpl(i) < 0 .and. lpl(i+1) > 0) ) then
            tmp = abs(lpl(i)-lpl(i+1))
            if ( tmp  >  hiv ) then
               hiv = tmp
               infp  = i+1
            end if
         end if
      end do
    end subroutine inflection

  end function find_inflx

  subroutine lmultifab_owner_mask(mask, all)
    type(lmultifab), intent(inout) :: mask
    logical, intent(in), optional :: all
    integer :: i, j
    type(box) :: bxi, bxj, bxij
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    call setval(mask, val = .true.)
    do i = 1, mask%nboxes; if ( remote(mask, i) ) cycle
       if ( lall ) then
          bxi = get_pbox(mask, i)
       else
          bxi = get_ibox(mask, i)
       end if
       do j = 1, i-1
          if ( lall ) then
             bxj = get_pbox(mask, j)
          else
             bxj = get_ibox(mask, j)
          end if
          bxij = intersection(bxi, bxj)
          if ( empty(bxij) ) cycle
          call setval(mask%fbs(i), .false., bxij)
       end do
    end do
  end subroutine lmultifab_owner_mask

end module cluster_module
