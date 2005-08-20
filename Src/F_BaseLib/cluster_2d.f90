module cluster2d_module

  use bl_types
  use bl_error_module
  use box_module
  use list_box_module
  use multifab_module

  implicit none

  logical, private, save :: verbose = .false.
  integer, private, save :: cut_threshold = 2

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

  subroutine bit_print(tagbox)
    logical, intent(in) :: tagbox(0:,0:)
    integer :: i, j
    integer :: nx, ny
    character c
    nx = size(tagbox,1)
    ny = size(tagbox,2)
    do j = ny-1, 0, -1
       write(*,fmt='(i3," : ")', advance='no') j
       do i = 0, nx-1
          c = ' '; if ( tagbox(i,j) ) c = '*'
          write(*,fmt='(a1)', advance = 'no') c
       end do
       write(*,fmt='()')
    end do
  end subroutine bit_print

  subroutine cls_2d_mf(boxes, tagboxes, minwidth, buf_wid, min_eff, crse_domain)
    type(list_box), intent(out) :: boxes
    type(lmultifab), intent(in) :: tagboxes
    integer, intent(in) :: minwidth
    integer, intent(in) :: buf_wid
    real(dp_t), intent(in) :: min_eff
    type(box), intent(in) :: crse_domain

    type(lmultifab) :: buf
    logical, pointer :: mt(:,:,:,:)
    logical, pointer :: mb(:,:,:,:)
    integer :: b
    integer :: num_flag
    integer :: lo(2),lo_d(2),hi_d(2)
    
    call build(boxes)
    call build(buf, get_layout(tagboxes), nc = 1, ng = buf_wid)

    num_flag = lmultifab_count(tagboxes)

    if ( num_flag == 0 ) then
       call bl_warn("No cells are flagged/No boxes returned")
       return
    end if

    lo_d = lwb(crse_domain)
    hi_d = upb(crse_domain)

    do b = 1, buf%nboxes; if ( remote(buf, b) ) cycle
       mt => dataptr(tagboxes, b)
       mb => dataptr(buf, b)
       lo = lwb(get_box(tagboxes,b))
       call frob(mb(:,:,1,1), mt(:,:,1,:), lo, buf_wid, lo_d, hi_d)
       print *,'TAG BOX ',lwb(get_pbox(tagboxes,b)),upb(get_pbox(tagboxes,b))
       print *,'BUF BOX ',lwb(get_pbox(buf,b)),upb(get_pbox(buf,b))
    end do

    call internal_sync(buf, all = .true., filter = filter_lor)
    call owner_mask(buf)

    call cluster_mf(boxes, buf, minwidth, min_eff)

    call destroy(buf)

  contains
    
    subroutine owner_mask(mask)
      type(lmultifab), intent(inout) :: mask
      integer :: i, j
      type(box) :: bxi, bxj, bxij

!     call setval(mask, val = .true.)
      do i = 1, mask%nboxes; if ( remote(mask, i) ) cycle
         bxi = get_pbox(mask, i)
         do j = 1, i-1
            bxj = get_pbox(mask, j)
            bxij = intersection(bxi, bxj)
            if ( empty(bxij) ) cycle
            call setval(mask%fbs(i), .false., bxij)
         end do
      end do
    end subroutine owner_mask

    subroutine frob(bb, tt, lo, ng, lo_d, hi_d)
      integer, intent(in   ) :: lo(:),ng
      logical, intent(in   ) :: tt(lo(1)   :,lo(2)   :,:)
      logical, intent(  out) :: bb(lo(1)-ng:,lo(2)-ng:)
      integer, intent(in   ) :: lo_d(:),hi_d(:)

      integer :: i, j, l, m, nx,ny

      ny = size(tt,2); nx = size(tt,1)
      print *,'SIZE OF BB ',nx,ny
      print *,'SIZE OF TT ',size(bb,1),size(bb,2)

      do j = lo(2),lo(2)+ny-1; 
      do i = lo(1),lo(1)+nx-1; 
         if ( any(tt(i,j,:)) ) then
            do l = max(1,i-ng), min(i+ng,nx)
               do m = max(1,j-ng), min(j+ng,ny)
                  bb(l,m) = .true.
               end do
            end do
         end if
      end do
      end do

      ! Dont allow grids with only one crse cell between fine and the
      !   physical boundary
      do j = lo(2),lo(2)+ny-1
        if (bb(lo(1)+   1,j)) bb(lo(1)     ,j) = .true.
        if (bb(lo(1)+nx-2,j)) bb(lo(1)+nx-1,j) = .true.
      end do
      do i = lo(1),lo(1)+nx-1
        if (bb(i,lo(2)+   1)) bb(i,lo(2)     ) = .true.
        if (bb(i,lo(2)+ny-2)) bb(i,lo(2)+ny-1) = .true.
      end do

      ! Dont allow tagging outside the physical domain
      do j = lo(2),lo(2)+ny-1; 
        bb(lo(1)-ng:lo(1)-1      ,j) = .false.
        bb(lo(1)+nx:lo(1)+nx+ng-1,j) = .false.
      end do
      do i = lo(1),lo(1)+nx-1; 
        bb(i,lo(2)-ng:lo(2)-1      ) = .false.
        bb(i,lo(2)+ny:lo(2)+ny+ng-1) = .false.
      end do

    end subroutine frob

  end subroutine cls_2d_mf

  subroutine cluster_mf(boxes, tagboxes, minwidth, min_eff)
    type(list_box), intent(inout) :: boxes
    type(lmultifab), intent(in) ::  tagboxes
    integer, intent(in) ::  minwidth
    real(dp_t), intent(in) ::   min_eff
    real(dp_t) cur_box_eff
    type(list_box_node), pointer :: bn
    integer, allocatable :: sigx(:), sigy(:)
    integer, allocatable :: lplx(:), lply(:)
    logical :: flag
    type(box) :: bbx
    integer :: lo(2), hi(2)

    lo = 0 ; hi = 0
    bbx = bbox(get_boxarray(tagboxes))
    print *,'BBX ',bbx
    lo(1:bbx%dim) = lwb(bbx); hi(1:bbx%dim) = upb(bbx)
    allocate(sigx(lo(1):hi(1)), sigy(lo(2):hi(2)))
    allocate(lplx(lo(1):hi(1)), lply(lo(2):hi(2)))
    call push_back(boxes, bbx)
    bn => begin(boxes)
    print *,'BN ',value(bn)
    do while ( associated(bn) )
    print *,'DOING EFF '
       cur_box_eff = box_eff_mf(tagboxes, value(bn))
    print *,'EFF IS ',cur_box_eff
       if ( verbose  .and. parallel_IOProcessor() ) then
          print *, 'cur_box_eff', cur_box_eff
       end if
       if ( cur_box_eff  <  min_eff ) then
          call sigma_laplace_mf_2d(tagboxes, value(bn), sigx, sigy, lplx, lply, lo)
          flag = find_split_2d(boxes, bn, minwidth, sigx, sigy, lplx, lply)
          if ( .not. flag ) then
             bn => next(bn)
          end if
       else
          bn => next(bn)
       end if
    end do

  end subroutine cluster_mf

  subroutine sigma_laplace_mf_2d(tagboxes, bx, sigx, sigy, lplx, lply, lo)
    type(lmultifab), intent(in) :: tagboxes
    integer, intent(in):: lo(:)
    type(box), intent(in) :: bx
    integer, intent(out) :: sigx(lo(1):),sigy(lo(2):)
    integer, intent(out) :: lplx(lo(1):),lply(lo(2):)
    logical, pointer :: tp(:,:,:,:)
!   logical, pointer :: mp(:,:,:,:)
    integer :: n

    integer, allocatable :: tsigx(:),tsigy(:)
    integer :: i, j, k, lx, ly, hx, hy

    sigx = 0
    sigy = 0

    lplx = 0
    lply = 0

    lx = lwb(bx,1); ly = 1
    hx = upb(bx,1); hy = 1
    if ( tagboxes%dim > 1 ) then
       ly = lwb(bx,2)
       hy = upb(bx,2)
    end if

    allocate(tsigx(lx:hx), tsigy(ly:hy))
    tsigx = 0
    tsigy = 0

    do n = 1, tagboxes%nboxes; if ( remote(tagboxes, n) ) cycle
       if ( .not. intersects(get_box(tagboxes,n), bx) ) cycle
       tp => dataptr(tagboxes, n)
       do i = lx, hx
          tsigx(i) = tsigx(i) + count(tp(i,ly:hy,1,1))
       end do

       do j = ly, hy
          ! tsigy(j) = tsigy(j) + count(tp(lx:hx,j,1,1) .and. mp(lx:hx,j,1,1))
          tsigy(j) = tsigy(j) + count(tp(lx:hx,j,1,1))
       end do
    end do

    call parallel_reduce(sigx(lx:hx), tsigx, MPI_SUM)
    call parallel_reduce(sigy(ly:hy), tsigy, MPI_SUM)

       !! Note: only one of berger/rigotsis schemes is here used.
       !! Note: a fill boundary needs to have b???
    do i = lx+1, hx-1
       lplx(i) = lplx(i) + sigx(i+1)-2*sigx(i)+sigx(i-1)
    end do

    do j = ly+1, hy-1
       lply(j) = lply(j) + sigy(j+1)-2*sigy(j)+sigy(j-1)
    end do

!   if ( verbose .and. parallel_IOProcessor() ) then
!      print '(a,1x,20(i3,1x))', 'sigx', sigx
!      print '(a,1x,20(i3,1x))', 'sigy', sigy
!      print '(a,1x,20(i3,1x))', 'lplx', lplx
!      print '(a,1x,20(i3,1x))', 'lply', lply
!   end if

  end subroutine sigma_laplace_mf_2d

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
    print *,'R / VOL IN BOX_EFF ',r,dvolume(bx)
    r = r/dvolume(bx)
  end function box_eff_mf

  function box_eff(tagbox, bx) result(r)
    real(dp_t) :: r

    type(box), intent(in) :: bx
    logical, intent(in) :: tagbox(0:,0:)

    integer :: lx, ly, hx, hy

    lx = lwb(bx,1); ly = lwb(bx,2)
    hx = upb(bx,1); hy = upb(bx,2)

    r = real(count(tagbox(lx:hx, ly:hy)),dp_t)/dvolume(bx)

  end function box_eff

  function find_split_2d(boxes, bn, minwidth, sigx, sigy, lplx, lply) result(r)
    logical :: r, rr
    type(list_box), intent(inout) ::  boxes
    type(list_box_node), pointer :: bn
    integer, intent(in) :: minwidth
    integer, intent(in) :: sigx(0:), sigy(0:)
    integer, intent(in) :: lplx(0:), lply(0:)

    rr = find_holes(boxes, bn, minwidth, sigx, sigy)
    if ( verbose  .and. parallel_IOProcessor() ) then
       print *, 'FIND_SPLIT(1) r = ', rr
    end if
    if ( .not. rr ) then
       rr = find_inflx(boxes, bn, minwidth, lplx, lply)
    end if
    if ( verbose  .and. parallel_IOProcessor() ) then
       print *, 'FIND_SPLIT(2) r = ', rr
    end if
    r = rr

  end function find_split_2d

  function find_holes(boxes, bn, minwidth, sigx, sigy) result(r)
    logical :: r
    type(list_box), intent(inout) ::  boxes
    type(list_box_node), pointer :: bn
    integer, intent(in) :: minwidth
    integer, intent(in) :: sigx(0:), sigy(0:)

    type(box) :: bx
    integer :: lx, hx, ly, hy

    if ( verbose  .and. parallel_IOProcessor() ) then
       call print(value(bn), 'find_holes')
    end if

    r  = .false.

    bx = value(bn)

    lx = lwb(value(bn),1); ly = lwb(value(bn),2)
    hx = upb(value(bn),1); hy = upb(value(bn),2)

    call bracket(sigx, lx, hx)
    call bracket(sigy, ly, hy)

    call set(bn, make_box((/lx, ly/), (/hx, hy/)))

    if ( holes(sigx, lx, hx, 1) ) then
       r = .true.
    else if ( holes(sigy, ly, hy, 2) ) then
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

  function find_inflx(boxes, bn, minwidth, lplx, lply) result(r)
    logical :: r
    type(list_box), intent(inout) ::  boxes
    type(list_box_node), pointer :: bn
    integer, intent(in) :: minwidth
    integer, intent(in) :: lplx(0:), lply(0:)

    type(box) :: b1, b2
    integer :: lx, ly, hx, hy, wx, wy
    integer :: hix, ipx
    integer :: hiy, ipy

    r = .false.

    lx = lwb(value(bn),1); ly = lwb(value(bn),2)
    hx = upb(value(bn),1); hy = upb(value(bn),2)

    call inflection(lplx, lx, hx, hix, ipx)
    call inflection(lply, ly, hy, hiy, ipy)

    if ( verbose  .and. parallel_IOProcessor() ) then
       call print(value(bn), 'FIND_INFLX')
       print '(a,1x,20(i3,1x))', 'lplx', lplx
       print *, '      l', lx, 'h', hx, 'hiv', hix, 'inf', ipx
       print '(a,1x,20(i3,1x))', 'lply', lply
       print *, '      l', ly, 'h', hy, 'hiv', hiy, 'inf', ipy
    end if


    if (ipx /= -1 .or. ipy /= -1) then
       if ( max(hix, hiy) >= CUT_THRESHOLD ) then
          if ( hix >= hiy ) then
             call box_chop(value(bn), b1, b2, 1, ipx)
          else 
             call box_chop(value(bn), b1, b2, 2, ipy)
          end if
          call boxStack(boxes, bn, b1, b2)
          r  = .true.
       end if
    else if ( .false. ) then
       wx = hx - lx + 1
       wy = hy - ly + 1
       if ( wx >= wy ) then
          wx = lx + wx/2
          call box_chop(value(bn), b1, b2, 1, wx)
       else 
          wy = ly + wy/2
          call box_chop(value(bn), b1, b2, 2, wy)
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

end module cluster2d_module
