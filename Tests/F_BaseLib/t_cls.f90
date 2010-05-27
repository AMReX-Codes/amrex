subroutine t_cluster
  use f2kcli
  use bl_types
  use cluster_module
  integer :: buf_wid, minwidth
  real(dp_t) :: min_eff

  min_eff = .8_dp_t
  min_eff = .5_dp_t
  min_eff = .9_dp_t
  buf_wid = 2
  buf_wid = 0
  minwidth = 2

  call cluster_set_verbose(.false.)
  call cluster_set_verbose(.true.)
  call cluster_set_beta(1.0_dp_t)
  call t_cls_mf(buf_wid, minwidth, min_eff)

end subroutine t_cluster

subroutine t_cls_mf(buf_wid, minwidth, min_eff)
  use cluster_module
  use multifab_module
  use fabio_module
  use bl_constants_module
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
  type(multifab) :: mf

  dm = 2

  ng = 2
  select case (dm)
  case (2)
     select case (ng)
     case (1)
        bx(1) = make_box((/0,0/), (/n-1,n-1/))
     case (2)
        bx(1) = make_box((/0,0/), (/n/2-1,n-1/))
        bx(2) = make_box((/n/2,0/), (/n-1,n-1/))
     end select
  case (3)
     select case (ng)
     case (1)
        bx(1) = make_box((/0,0,0/), (/n-1,n-1,n-1/))
     case (2)
        bx(1) = make_box((/0,0,0/), (/n-1,n-1,n/2-1/))
        bx(2) = make_box((/0,0,n/2/), (/n-1,n-1,n-1/))
     end select
  case default
     call bl_error("NO OTHER CASES")
  end select
  call build(ba, bx(1:ng))
  call build(la, ba)
  call build(tagbox, la, nc=1, ng=0)
  call build(mf, la)

  wid1 = real(n,dp_t)/3
  wid2 = real(n,dp_t)/4
  cen  = real(n,dp_t)/2

  call setval(tagbox, .false.)

  do ng = 1, nboxes(tagbox)
     lo = 0; lo(1:dm) = lwb(get_box(tagbox, ng))
     hi = 0; hi(1:dm) = upb(get_box(tagbox, ng))
     lp => dataptr(tagbox, ng)
     select case (dm)
     case (1)
        do i = lo(1), hi(1)
           d = sqrt((i-cen)**2)
           if ( d <= wid1 .and. d >= wid2 ) then
              lp(i,1,1,1) = .true.
           end if
        end do
     case (2)
        do j = lo(2), hi(2); do i = lo(1), hi(1)
           d = sqrt((i-cen)**2 + (j-cen)**2)
           if ( d <= wid1 .and. d >= wid2 ) then
              lp(i,j,1,1) = .true.
           end if
        end do; end do
     case(3)
        do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
           d = sqrt((i-cen)**2 + (j-cen)**2 + (k-cen)**2)
           if ( d <= wid1 .and. d >= wid2 ) then
              lp(i,j,k,1) = .true.
           end if
        end do; end do; end do
     end select
  end do

  call setval_mask(mf, tagbox, ONE, ZERO)
  call fabio_multifab_write_d(mf, "tdir1", "tags1")

!  call cluster(boxes, tagbox, minwidth, buf_wid, min_eff, overall_eff)

  call cluster(boxes, tagbox, buf_wid, overall_eff)

  print *, 'number of boxes ', nboxes(boxes)
  do i = 1, nboxes(boxes)
     write(6,fmt='(i5,1x)', advance = 'no') i
     call print(get_box(boxes,i))
  end do
  print *, 'overall_eff', overall_eff

  call destroy(mf)
  call destroy(tagbox)
  call destroy(la)
  call destroy(ba)

  call destroy(boxes)

contains

  subroutine setval_mask(mf, tags, on, off)
    type(multifab), intent(inout) :: mf
    type(lmultifab), intent(in) :: tags
    real(dp_t), intent(in) :: on, off
    integer :: i
    logical, pointer :: tp(:,:,:,:)
    real(dp_t), pointer :: mp(:,:,:,:)
    integer :: lo(3)

    call setval(mf, off)
    do i = 1, tags%nboxes; if ( remote(tags, i) ) cycle
       tp => dataptr(tags, i)
       mp => dataptr(mf, i)
       where(tp) mp = on
    end do
  end subroutine setval_mask

end subroutine t_cls_mf
