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
  use bl_error_module
  use ml_boxarray_module
  use box_util_module
  use layout_module
  use bl_prof_module
  implicit none
  type(boxarray) :: ba
  integer :: dm
  integer :: ng
  integer, allocatable, dimension(:) :: ext, plo, phi, vsz, crsn, vrng
  type(ml_boxarray) :: mba
  character(len=64) :: test_set
  integer :: i, n, j, k, sz
  type(box) :: bx, cbx
  type(bl_prof_timer), save :: bpt, bpt_r, bpt_s, bpt_b
  integer :: cnt, cnt1, cnt2, cnt3
  integer(ll_t) :: vol, vol1, vol2, vol3
  integer, pointer :: ipv(:)
  integer :: mxsz
  logical :: verbose
  type(layout) :: la

  type bin
     integer, pointer :: iv(:) => Null()
  end type bin

  type(bin), allocatable, dimension(:,:,:) :: bins

  call build(bpt, "t_ba_self_intersection")

  verbose = .true.
  ng = 1
  test_set = "grids.5034"
  test_set = "grids.1071"
  test_set = "grids.213"
  test_set = "grids.60"

  call build(bpt_r, "ba_read")
  call read_a_mglib_grid(mba, test_set)
  call destroy(bpt_r)

  ba = mba%bas(1)
  call build(la, ba, pd = mba%pd(1))

  dm = get_dim(ba)

  allocate(ext(dm), plo(dm), phi(dm), vsz(dm), crsn(dm), vrng(dm))

  cnt = 0; cnt1 = 0; cnt2 = 0; cnt3 = 0
  vol = 0; vol1 = 0; vol2 = 0; vol3 = 0
  crsn = 32

  call init_box_hash_bin(la)

  call build(bpt_b, "build hash")
  bx = boxarray_bbox(ba)
  cbx = coarsen(bx,crsn)
  plo = lwb(cbx)
  phi = upb(cbx)

  vsz = -Huge(1)
  do n = 1, nboxes(ba)
     vsz = max(vsz,extent(get_box(ba,n)))
  end do
  print *, 'max extent', vsz
  print *, 'crsn max extent', int_coarsen(vsz,crsn+1)
  vrng = int_coarsen(vsz,crsn+1)
  print *, 'vrng = ', vrng

  allocate(bins(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3)))
  do k = plo(3), phi(3); do j = plo(2), phi(2); do i = plo(1), phi(1)
     allocate(bins(i,j,k)%iv(0))
  end do;end do; end do

  do n = 1, nboxes(ba)
     ext = int_coarsen(lwb(get_box(ba,n)),crsn)
     if ( .not. contains(cbx, ext) ) then
        call bl_error("Not Contained!")
     end if
     sz = size(bins(ext(1),ext(2),ext(3))%iv)
     allocate(ipv(sz+1))
     ipv(1:sz) = bins(ext(1),ext(2),ext(3))%iv(1:sz)
     ipv(sz+1) = n
     deallocate(bins(ext(1),ext(2),ext(3))%iv)
     bins(ext(1),ext(2),ext(3))%iv => ipv
  end do
  call destroy(bpt_b)

  if ( verbose ) then
     call print(bx, 'bbox(ba) ')
     call print(cbx, 'coarsen(bbox(ba),crsn) ')
     print *, 'extents(bx)', extent(bx)
     print *, 'extents(cbx)', extent(cbx)
     print *, 'plo ', plo
     print *, 'phi ', phi
     mxsz = -Huge(1)
     do k = plo(3), phi(3); do j = plo(2), phi(2); do i = plo(1), phi(1)
        mxsz = max(mxsz, size(bins(i,j,k)%iv))
     end do;end do; end do
     print *, 'max bin sz ', mxsz

     sz = Huge(1)
     do k = plo(3), phi(3); do j = plo(2), phi(2); do i = plo(1), phi(1)
        sz = min(sz, size(bins(i,j,k)%iv))
     end do;end do; end do
     print *, 'min bin sz ', sz

     sz = 0
     do k = plo(3), phi(3); do j = plo(2), phi(2); do i = plo(1), phi(1)
        sz = sz + size(bins(i,j,k)%iv)
     end do;end do; end do
     print *, 'tot bins ', sz

     if ( sz /= nboxes(ba) ) then
        call bl_error("sz /= nboxes(ba): ", sz)
     end if
  end if

  call build(bpt_s, "ba_s")
  do i = 1, nboxes(ba)
     bx = grow(get_box(ba,i), ng)
     call self_intersection(bx, ba)
     call self_intersection_1(bx, ba)
     call la_chk_box(la, bx)
  end do
  call destroy(bpt_s)

  ! Just a check of the result
  print *, 'cnt = ', cnt
  print *, 'cnt1 = ', cnt1
  print *, 'cnt2 = ', cnt2
  print *, 'cnt3 = ', cnt3
  print *, 'vol = ', vol
  print *, 'vol1 = ', vol1
  print *, 'vol2 = ', vol2
  print *, 'vol3 = ', vol3

  do k = plo(3), phi(3); do j = plo(2), phi(2); do i = plo(1), phi(1)
     deallocate(bins(i,j,k)%iv)
  end do;end do; end do

  call destroy(la)
  call destroy(mba)
  call destroy(bpt)
contains

  subroutine la_chk_box(la, bx)
    type(layout), intent(inout) :: la
    type(box), intent(in) :: bx
    type(box_intersector), pointer :: bi(:)
    integer :: i
    type(bl_prof_timer), save :: bpt_h1
    call build(bpt_h1, "ba_h1")
    bi => layout_get_box_intersector(la, bx)
    do i = 1, size(bi)
       cnt3 = cnt3 + 1
       vol3 = vol3 + volume(bi(i)%bx)
    end do
    deallocate(bi)
    call destroy(bpt_h1)
  end subroutine la_chk_box

  subroutine self_intersection(bx, ba)
    type(box), intent(in) :: bx
    type(boxarray), intent(in) :: ba
    integer :: i
    type(bl_prof_timer), save :: bpt_i
    type(box) :: bx1
    call build(bpt_i, "ba_i")
    do i = 1, nboxes(ba)
       bx1 = intersection(bx, get_box(ba,i))
       if ( empty(bx1) ) cycle
       cnt = cnt + 1
       vol = vol + volume(bx1)
    end do
    call destroy(bpt_i)

  end subroutine self_intersection

  subroutine self_intersection_1(bx, ba)
    type(box), intent(in) :: bx
    type(boxarray), intent(in) :: ba
    integer :: i
    type(bl_prof_timer), save :: bpt_i
    type(box) :: bx1(nboxes(ba)),bx2(nboxes(ba))
    logical   :: is_empty(nboxes(ba))
    call build(bpt_i, "ba_i1")
    do i = 1, nboxes(ba)
       bx2(i) = get_box(ba,i)
    end do
    call box_intersection_and_empty(bx1, is_empty, bx, bx2)
    do i = 1, size(bx1)
       if ( is_empty(i) ) cycle
       vol1 = vol1 + volume(bx1(i))
    end do
    cnt1 = cnt1 + count(.not.is_empty)
    call destroy(bpt_i)
  end subroutine self_intersection_1

end subroutine t_ba_self_intersection

function log2(vin) result(r)
  implicit none
  integer :: r
  integer, intent(in) :: vin
  integer :: v
  r = 0
  v = vin
  do
     v = ishft(v, -1)
     if ( v == 0 ) exit
     r = r + 1
  end do
end function log2
