module cc_rhs_module

    use BoxLib
    use ml_layout_module
    use multifab_module

    implicit none

contains

  subroutine cc_rhs(mla, rh)

    type(ml_layout), intent(inout) :: mla
    type(box      )                :: pd
    type( multifab), intent(inout) :: rh(:)

    integer                        :: n, dm, nlevs

    dm = mla%dim

    nlevs = mla%nlevel

    do n = nlevs, 1, -1
       call setval(rh(n), val = 0.d0, all=.true.)
    end do

    call mf_init_2(rh(nlevs))

  end subroutine cc_rhs

  subroutine mf_init(mf,pd)

    type(multifab), intent(inout) :: mf
    type(box)     , intent(in   ) :: pd
    type(box) :: bx

    bx = pd
    bx%lo(1:bx%dim) = (bx%hi(1:bx%dim) + bx%lo(1:bx%dim))/2
    bx%hi(1:bx%dim) = bx%lo(1:bx%dim)
    call setval(mf, 1.d0, bx)

  end subroutine mf_init

  subroutine mf_init_2(mf)
    type(multifab), intent(inout) :: mf
    integer i
    type(box) bx
    do i = 1, nfabs(mf)

       bx = get_box(mf,i)
       bx%lo(1:bx%dim) = (bx%hi(1:bx%dim) + bx%lo(1:bx%dim))/2
       bx%hi(1:bx%dim) = bx%lo(1:bx%dim)
       call setval(mf%fbs(i), 1.0_dp_t, bx)
!       print *,'SETTING RHS TO 1.d0 ON ',bx%lo(1:bx%dim)

!      Single point of non-zero RHS: use this to make system solvable
       bx = get_box(mf,i)
       bx%lo(1       ) = (bx%hi(1       ) + bx%lo(1       ))/2 + 1
       bx%lo(2:bx%dim) = (bx%hi(2:bx%dim) + bx%lo(2:bx%dim))/2
       bx%hi(1:bx%dim) = bx%lo(1:bx%dim)
       call setval(mf%fbs(i), -1.0_dp_t, bx)
!       print *,'SETTING RHS TO -1.d0 ON ',bx%lo(1:bx%dim)

!      1-d Strip: Variation in x-direction
!      bx%lo(1) = (bx%hi(1) + bx%lo(1))/2
!      bx%hi(1) = bx%lo(1)+1

!      1-d Strip: Variation in y-direction
!      bx%lo(2) = (bx%hi(2) + bx%lo(2))/2
!      bx%hi(2) = bx%lo(2)+1

!      1-d Strip: Variation in z-direction
!      bx%lo(3) = (bx%hi(3) + bx%lo(3))/2
!      bx%hi(3) = bx%lo(3)+1

    end do
  end subroutine mf_init_2

  subroutine mf_init_1(mf)
    type(multifab), intent(inout) :: mf
    integer i
    type(box) bx
    type(box) rhs_box, rhs_intersect_box

    rhs_box%dim = mf%dim
    rhs_box%lo(1:rhs_box%dim) = 7
    rhs_box%hi(1:rhs_box%dim) = 8

    do i = 1, nfabs(mf)
       bx = get_ibox(mf,i)
       rhs_intersect_box = box_intersection(bx,rhs_box)
       if (.not. empty(rhs_intersect_box)) then
         bx%lo(1:bx%dim) = lwb(rhs_intersect_box)
         bx%hi(1:bx%dim) = upb(rhs_intersect_box)
!         print *,'SETTING RHS IN BOX ',i,' : ', bx%lo(1:bx%dim),bx%hi(1:bx%dim)
         call setval(mf%fbs(i), 1.d0, bx)
       end if
    end do

  end subroutine mf_init_1

end module cc_rhs_module
