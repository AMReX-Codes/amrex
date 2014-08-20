module edge_restriction_module

  use bl_constants_module
  use bl_error_module
  use bl_types
  use bc_functions_module

  implicit none

contains


  subroutine edge_restriction_1d(cc, loc, ff, lof, lo, hi, ir)
    integer,     intent(in)    :: loc(:)
    integer,     intent(in)    :: lof(:)
    integer,     intent(in)    :: lo(:), hi(:)
    real (dp_t), intent(inout) :: cc(loc(1):)
    real (dp_t), intent(in)    :: ff(lof(1):)
    integer,     intent(in)    :: ir(:)

    integer :: i

    do i = lo(1), hi(1)
       cc(i) = ff(ir(1)*i)
    end do

  end subroutine edge_restriction_1d

  subroutine edge_restriction_2d(cc, loc, ff, lof, lo, hi, ir, face)
    integer,     intent(in)    :: loc(:)
    integer,     intent(in)    :: lof(:)
    integer,     intent(in)    :: lo(:), hi(:)
    integer,     intent(in)    :: face
    real (dp_t), intent(inout) :: cc(loc(1):,loc(2):)
    real (dp_t), intent(in)    :: ff(lof(1):,lof(2):)
    integer,     intent(in)    :: ir(:)

    real (dp_t) :: fac
    integer     :: i, j, l, m

    fac = one

    if (face .eq. 1) then
       fac = fac/real(ir(2),kind=dp_t)
    else if (face .eq. 2) then
       fac = fac/real(ir(1),kind=dp_t)
    else
       call bl_error('edge_restriction_2d: face must be 1 or 2')
    end if

    if ( face .eq. 1 ) then
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             cc(i,j) = zero
             do m = 0, ir(2)-1
                cc(i,j) = cc(i,j) + ff(ir(1)*i,ir(2)*j+m)
             end do
             cc(i,j) = cc(i,j)*fac
          end do
       end do
    else if ( face .eq. 2 ) then
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             cc(i,j) = zero
             do l = 0, ir(1)-1
                cc(i,j) = cc(i,j) + ff(ir(1)*i+l,ir(2)*j)
             end do
             cc(i,j) = cc(i,j)*fac
          end do
       end do
    end if

  end subroutine edge_restriction_2d

  subroutine edge_restriction_3d(cc, loc, ff, lof, lo, hi, ir, face)
    integer,     intent(in)    :: loc(:)
    integer,     intent(in)    :: lof(:)
    integer,     intent(in)    :: lo(:),hi(:)
    integer,     intent(in)    :: face
    real (dp_t), intent(inout) :: cc(loc(1):,loc(2):,loc(3):)
    real (dp_t), intent(in)    :: ff(lof(1):,lof(2):,lof(3):)
    integer,     intent(in)    :: ir(:)

    real (dp_t) :: fac
    integer     :: i, j, k, l, m, n

    fac = one

    if ( face .eq. 1 ) then
       fac = fac/real(ir(2)*ir(3),kind=dp_t)
    else if (face .eq. 2) then
       fac = fac/real(ir(1)*ir(3),kind=dp_t)
    else if (face .eq. 3) then
       fac = fac/real(ir(1)*ir(2),kind=dp_t)
    else
       call bl_error('edge_restriction_3d: face must be 1, 2 or 3')
    end if

    if ( face .eq. 1 ) then
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                cc(i,j,k) = zero
                do n = 0, ir(3)-1
                   do m = 0, ir(2)-1
                      cc(i,j,k) = cc(i,j,k) + ff(ir(1)*i,ir(2)*j+m,ir(3)*k+n)
                   end do
                end do
                cc(i,j,k) = cc(i,j,k)*fac
             end do
          end do
       end do
    else if ( face .eq. 2 ) then
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                cc(i,j,k) = zero
                do n = 0, ir(3)-1
                   do l = 0, ir(1)-1
                      cc(i,j,k) = cc(i,j,k) + ff(ir(1)*i+l,ir(2)*j,ir(3)*k+n)
                   end do
                end do
                cc(i,j,k) = cc(i,j,k)*fac
             end do
          end do
       end do
    else if ( face .eq. 3 ) then
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                cc(i,j,k) = zero
                do m = 0, ir(2)-1
                   do l = 0, ir(1)-1
                      cc(i,j,k) = cc(i,j,k) + ff(ir(1)*i+l,ir(2)*j+m,ir(3)*k)
                   end do
                end do
                cc(i,j,k) = cc(i,j,k)*fac
             end do
          end do
       end do
    end if

  end subroutine edge_restriction_3d

end module edge_restriction_module
