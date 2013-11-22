module nodal_rhs_module

    use BoxLib
    use ml_layout_module
    use multifab_module
    use mg_module

    implicit none

contains

    subroutine nodal_rhs(mla, rh)

      type(ml_layout), intent(inout) :: mla
      type( multifab), intent(inout) :: rh(:)

      integer                        :: n, dm, nlevs

      dm = mla%dim

      nlevs = mla%nlevel

      do n = nlevs, 1, -1
         call setval(rh(n), val = ZERO, all=.true.)
      end do

      call mf_init(rh(nlevs))

    end subroutine nodal_rhs

  subroutine mf_init(mf)
    type(multifab), intent(inout) :: mf
    integer i
    type(box) bx
    do i = 1, nfabs(mf)

!      Single point of non-zero RHS
       bx = get_ibox(mf,i)
       bx%lo(1:bx%dim) = (bx%hi(1:bx%dim) + bx%lo(1:bx%dim))/2
       bx%hi(1:bx%dim) = bx%lo(1:bx%dim)
       call setval(mf%fbs(i), ONE, bx)

!      Single point of non-zero RHS: use this to make system solvable
       bx = get_ibox(mf,i)
       bx%lo(1       ) = (bx%hi(1       ) + bx%lo(1       ))/2 + 1
       bx%lo(2:bx%dim) = (bx%hi(2:bx%dim) + bx%lo(2:bx%dim))/2
       bx%hi(1:bx%dim) = bx%lo(1:bx%dim)
       call setval(mf%fbs(i), -ONE, bx)

    end do
  end subroutine mf_init

  subroutine mf_init1(mf)
    type(multifab), intent(inout) :: mf
    integer i
    type(box) bx
    type(box) rhs_box, rhs_intersect_box

    rhs_box%dim = mf%dim
    rhs_box%lo(1:rhs_box%dim) = 8
    rhs_box%hi(1:rhs_box%dim) = 8

    do i = 1, nfabs(mf)
       bx = get_ibox(mf,i)
       rhs_intersect_box = box_intersection(bx,rhs_box)
       if (.not. empty(rhs_intersect_box)) then
         bx%lo(1:bx%dim) = lwb(rhs_intersect_box)
         bx%hi(1:bx%dim) = upb(rhs_intersect_box)
         call setval(mf%fbs(i), ONE, bx)
       end if
    end do

  end subroutine mf_init1

end module nodal_rhs_module
