subroutine t_cluster
  use f2kcli
  use bl_types
  use cluster_module
  integer :: buf_wid, minwidth
  real(dp_t) :: min_eff

  min_eff = .7
  buf_wid = 1
  minwidth = 1

  call cluster_set_verbose(.true.)
  call t_cls_mf(buf_wid, minwidth, min_eff)

end subroutine t_cluster

subroutine t_cls_mf(buf_wid, minwidth, min_eff)
  use cluster_module
  use multifab_module
  implicit none
  integer, intent(in) :: minwidth
  integer, intent(in) :: buf_wid
  real(dp_t), intent(in) :: min_eff
  integer, parameter :: n = 16
  type(lmultifab) :: tagbox
  type(boxarray) :: boxes, ba
  type(box) :: bx(2)
  type(layout) :: la
  real(dp_t) :: overall_eff
  real(dp_t) :: d, wid1, wid2, cen
  integer :: i, j, k, dm
  logical, pointer :: lp(:,:,:,:)
  integer :: lo(3), hi(3), ng

  dm = 2

  ng = 2
  select case (ng)
  case (1)
     bx(1) = make_box((/0,0/), (/n-1,n-1/))
  case (2)
     bx(1) = make_box((/0,0/), (/n/2-1,n-1/))
     bx(2) = make_box((/n/2,0/), (/n-1,n-1/))
  end select
  call build(ba, bx(1:ng))
  call build(la, ba)
  call build(tagbox, la, nc=1, ng=0)

  wid1 = real(n,dp_t)/3
  wid2 = real(n,dp_t)/4
  cen  = real(n,dp_t)/2

  call setval(tagbox, .false.)

  do ng = 1, nboxes(tagbox)
     lo = 0; lo(1:dm) = lwb(get_box(tagbox, ng))
     hi = 0; hi(1:dm) = upb(get_box(tagbox, ng))
     lp => dataptr(tagbox, ng)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           if ( .true. ) then
              if (  i == 7 .and. j == 7 ) then
                 lp(i,j,1,1) = .true.
              end if
              if (  i == 8 .and. j == 8 ) then
                 lp(i,j,1,1) = .true.
              end if
           else
              d = sqrt((i-cen)**2 + (j-cen)**2)
              if ( d <= wid1 .and. d >= wid2 ) then
                 lp(i,j,1,1) = .true.
              end if
           endif
        end do
     end do
  end do

  call cluster(boxes, tagbox, minwidth, buf_wid, min_eff, overall_eff)

  print *, 'number of boxes ', nboxes(boxes)
  do i = 1, nboxes(boxes)
     call print(get_box(boxes,i))
  end do
  print *, 'overall_eff', overall_eff

  call destroy(tagbox)
  call destroy(la)
  call destroy(ba)

  call destroy(boxes)

end subroutine t_cls_mf

subroutine t_cls(buf_wid, minwidth, min_eff)
  use cluster_module
  use multifab_module
  implicit none
  integer, intent(in) :: minwidth
  integer, intent(in) :: buf_wid
  real(dp_t), intent(in) :: min_eff
  integer, parameter :: n = 16
  logical :: tagbox(0:n-1, 0:n-1, 0:0)
  type(boxarray) :: boxes
  type(list_box_node), pointer :: bn
  real(dp_t) :: overall_eff
  real(dp_t) :: d, wid1, wid2, cen
  integer :: i, j, k, dm

  dm = 2
  wid1 = real(n,dp_t)/3
  wid2 = real(n,dp_t)/4
  cen  = real(n,dp_t)/2

  tagbox = .false.

  k = 0
  do j = 0, n-1
     do i = 0, n-1
        d = sqrt((i-cen)**2 + (j-cen)**2)
        if ( d <= wid1 .and. d >= wid2 ) then
           tagbox(i,j,k) = .true.
        end if
     end do
  end do

  call cluster(boxes, tagbox, minwidth, buf_wid, min_eff, overall_eff, dim = dm)

  print *, 'number of boxes ', nboxes(boxes)
  do i = 1, nboxes(boxes)
     call print(get_box(boxes,i))
  end do
  print *, 'overall_eff', overall_eff

  call destroy(boxes)

end subroutine t_cls
