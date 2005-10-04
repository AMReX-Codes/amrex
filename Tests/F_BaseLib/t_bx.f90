subroutine t_bx
  use box_module
  implicit none
  type(box) :: bx
  bx = allbox(2)
  print *, bx
  bx = allbox_grc(bx, grow = 2)
  print *, bx
contains
  function allbox_grc(bx, grow, refine, coarsen) result(r)
    type(box) :: r
    type(box), intent(in) :: bx
    integer, intent(in), optional :: grow, refine, coarsen
    integer :: dm
    r = bx
    dm = r%dim
    if ( present(grow) ) then
    else if ( present(refine) ) then
    else if ( present(coarsen) ) then
    end if
  end function allbox_grc
end subroutine t_bx

subroutine t_ba_self_intersection
  use ml_boxarray_module
  use box_util_module
  use bl_prof_module
  implicit none
  type(boxarray) :: ba
  integer :: dm
  integer :: ng
  integer :: pd_xyz(MAX_SPACEDIM), ext(MAX_SPACEDIM), pext(MAX_SPACEDIM)
  integer :: ba_maxsize
  logical pmask(MAX_SPACEDIM)
  type(ml_boxarray) :: mba
  character(len=64) :: test_set
  integer :: i, f, n, j, k, sz
  type(box) :: bx
  type(bl_prof_timer), save :: bpt, bpt_r, bpt_s
  integer :: cnt, cnt1, cnt2
  integer, pointer :: ipv(:)
  integer :: tsz, mxsz

  integer :: crsn

  type bin
     integer, pointer :: iv(:) => Null()
  end type bin

  type(bin), allocatable, dimension(:,:,:) :: bins
  
  call build(bpt, "t_ba_self_intersection")

  ng = 1
  test_set = "grids.5034"

  call build(bpt_r, "ba_read")
  call read_a_mglib_grid(mba, test_set)
  call destroy(bpt_r)

  ba = mba%bas(1)

  cnt = 0
  cnt1 = 0
  cnt2 = 0
  call build(bpt_s, "ba_s")

  crsn = 32

  bx = boxarray_bbox(ba)
  print *, 'pd ', mba%pd(1)
  print *, 'cpd ', coarsen(mba%pd(1),crsn)

  pext = extent(coarsen(mba%pd(1),crsn))
  print *, 'prod(ext)', product(pext)

  allocate(bins(0:pext(1)-1,0:pext(2)-1,0:pext(3)-1))

  do k = 0, pext(3)-1; do j = 0, pext(2)-1; do i = 0, pext(1)-1
     allocate(bins(i,j,k)%iv(0))
  end do;end do; end do

  sz = 0
  do k = 0, pext(3)-1; do j = 0, pext(2)-1; do i = 0, pext(1)-1
     sz = sz + size(bins(i,j,k)%iv)
  end do;end do; end do
  print *, 'tot bins ', sz

  print *, 'size(bins)', size(bins)

  print *, 'bbox(ba) ', bx
  print *, 'extents(bx)', extent(bx)

  do n = 1, nboxes(ba)
     ext = lwb(coarsen(get_box(ba,n),crsn))
     sz = size(bins(ext(1),ext(2),ext(3))%iv)
     allocate(ipv(sz+1))
     ipv(1:sz) = bins(ext(1),ext(2),ext(3))%iv(1:sz)
     ipv(sz+1) = n
     if ( sz + 1 /= size(ipv) ) stop 'die'
     deallocate(bins(ext(1),ext(2),ext(3))%iv)
     bins(ext(1),ext(2),ext(3))%iv => ipv
     if ( sz + 1 /= size(bins(ext(1),ext(2),ext(3))%iv) ) stop 'die1'
  end do

  mxsz = -Huge(1)
  do k = 0, pext(3)-1; do j = 0, pext(2)-1; do i = 0, pext(1)-1
     mxsz = max(mxsz, size(bins(i,j,k)%iv))
  end do;end do; end do
  print *, 'max bin sz ', mxsz

  sz = Huge(1)
  do k = 0, pext(3)-1; do j = 0, pext(2)-1; do i = 0, pext(1)-1
     sz = min(sz, size(bins(i,j,k)%iv))
  end do;end do; end do
  print *, 'min bin sz ', sz

  sz = 0
  do k = 0, pext(3)-1; do j = 0, pext(2)-1; do i = 0, pext(1)-1
     sz = sz + size(bins(i,j,k)%iv)
  end do;end do; end do
  print *, 'tot bins ', sz

  do i = 1, nboxes(ba)
     bx = grow(get_box(ba,i),1)
     call self_intersection(bx, ba)
     call self_intersection_1(bx, ba)
     call chk_box(bx)
  end do
  call destroy(bpt_s)
  print *, 'cnt = ', cnt
  print *, 'cnt1 = ', cnt1
  print *, 'cnt2 = ', cnt2

  do k = 0, pext(3)-1; do j = 0, pext(2)-1; do i = 0, pext(1)-1
     deallocate(bins(i,j,k)%iv)
  end do;end do; end do

  call destroy(mba)
  call destroy(bpt)
contains

  subroutine chk_box(bx)
    type(box), intent(in) :: bx
    type(box) :: bx1
    integer :: lo(bx%dim), hi(bx%dim)
    integer :: i, j, k, n
    type(bl_prof_timer), save :: bpt_h
    call build(bpt_h, "ba_h")
    bx1 = coarsen(bx,crsn)
    lo = lwb(bx1)
    hi = upb(bx1)
    do k = max(lo(3)-1,0), min(hi(3), pext(3)-1)
       do j = max(lo(2)-1,0), min(hi(2), pext(2)-1)
          do i = max(lo(1)-1,0), min(hi(1), pext(1)-1)
             do n = 1, size(bins(i,j,k)%iv)
                bx1 = intersection(bx, ba%bxs(bins(i,j,k)%iv(n)))
                if ( empty(bx1) ) cycle
                cnt2 = cnt2 + 1
             end do
          end do
       end do
    end do
    call destroy(bpt_h)
  end subroutine chk_box

  subroutine self_intersection(bx, ba)
    type(box), intent(in) :: bx
    type(boxarray), intent(in) :: ba
    integer :: i
    type(bl_prof_timer), save :: bpt_i
    type(box) :: bx1
    call build(bpt_i, "ba_i")
    do i = 1, nboxes(ba)
       bx1 = intersection(bx, ba%bxs(i))
       if ( empty(bx1) ) cycle
       cnt = cnt + 1
    end do
    call destroy(bpt_i)

  end subroutine self_intersection

  subroutine self_intersection_1(bx, ba)
    type(box), intent(in) :: bx
    type(boxarray), intent(in) :: ba
    integer :: i
    type(bl_prof_timer), save :: bpt_i
    type(box) :: bx1(size(ba%bxs))
    logical   :: is_empty(size(ba%bxs))
    call build(bpt_i, "ba_i1")
    call box_intersection_and_empty(bx1, is_empty, bx, ba%bxs)
    call destroy(bpt_i)
    cnt1 = cnt1 + count(.not.is_empty)

  end subroutine self_intersection_1

end subroutine t_ba_self_intersection

