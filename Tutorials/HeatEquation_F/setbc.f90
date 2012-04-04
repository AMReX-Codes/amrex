module setbc_module

  use bl_types
  use bl_error_module

  implicit none

  private

  public :: setbc_2d, setbc_3d

contains

  subroutine setbc_2d(s,lo,ng,bc)

    use bl_constants_module
    use bc_module

    integer        , intent(in   ) :: lo(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)

    ! Local variables
    integer :: i,j,hi(2)

    hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)

    if (bc(1,1) .eq. EXT_DIR) then
       call bl_error('BC(1,1) = EXT_DIR NOT YET SUPPORTED')
    else if (bc(1,1) .eq. FOEXTRAP) then
       do j = lo(2)-ng, hi(2)+ng
          s(lo(1)-ng:lo(1)-1,j) = s(lo(1),j)
       end do
    else if (bc(1,1) .eq. HOEXTRAP) then
       do j = lo(2)-ng, hi(2)+ng
          s(lo(1)-ng:lo(1)-1,j) = EIGHTH* &
               (15.0_dp_t*s(lo(1),j) - 10.0_dp_t*s(lo(1)+1,j) + 3.0_dp_t*s(lo(1)+2,j))
       end do
    else if (bc(1,1) .eq. REFLECT_EVEN) then
       do j = lo(2)-ng, hi(2)+ng
          do i = 1, ng
             s(lo(1)-i,j) = s(lo(1)+i-1,j)
          end do
       end do
    else if (bc(1,1) .eq. REFLECT_ODD) then
       do j = lo(2)-ng, hi(2)+ng
          do i = 1, ng
             s(lo(1)-i,j) = -s(lo(1)+i-1,j)
          end do
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(1,1) = ',bc(1,1)
       call bl_error('BC(1,1) = NOT YET SUPPORTED')
    end if

    if (bc(1,2) .eq. EXT_DIR) then
       call bl_error('BC(1,1) = EXT_DIR NOT YET SUPPORTED')
    else if (bc(1,2) .eq. FOEXTRAP) then
       do j = lo(2)-ng, hi(2)+ng
          s(hi(1)+1:hi(1)+ng,j) = s(hi(1),j)
       end do
    else if (bc(1,2) .eq. HOEXTRAP) then
       do j = lo(2)-ng, hi(2)+ng
          s(hi(1)+1:hi(1)+ng,j) = EIGHTH* &
               (15.0_dp_t * s(hi(1)  ,j) - 10.0_dp_t*s(hi(1)-1,j) + 3.0_dp_t*s(hi(1)-2,j))
       end do
    else if (bc(1,2) .eq. REFLECT_EVEN) then
       do j = lo(2)-ng, hi(2)+ng
          do i = 1, ng
             s(hi(1)+i,j) = s(hi(1)-i+1,j)
          end do
       end do
    else if (bc(1,2) .eq. REFLECT_ODD) then
       do j = lo(2)-ng, hi(2)+ng
          do i = 1, ng
             s(hi(1)+i,j) = -s(hi(1)-i+1,j)
          end do
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(1,2) = ',bc(1,2)
       call bl_error('BC(1,2) = NOT YET SUPPORTED')
    end if

    if (bc(2,1) .eq. EXT_DIR) then
       call bl_error('BC(1,1) = EXT_DIR NOT YET SUPPORTED')
    else if (bc(2,1) .eq. FOEXTRAP) then
       do i = lo(1)-ng, hi(1)+ng
          s(i,lo(2)-ng:lo(2)-1) = s(i,lo(2))
       end do
    else if (bc(2,1) .eq. HOEXTRAP) then
       do i = lo(1)-ng, hi(1)+ng
          s(i,lo(2)-ng:lo(2)-1) = EIGHTH* &
               (15.0_dp_t*s(i,lo(2)) - 10.0_dp_t*s(i,lo(2)+1) + 3.0_dp_t*s(i,lo(2)+2))
       end do
    else if (bc(2,1) .eq. REFLECT_EVEN) then
       do i = lo(1)-ng, hi(1)+ng
          do j = 1, ng
             s(i,lo(2)-j) = s(i,lo(2)+j-1)
          end do
       end do
    else if (bc(2,1) .eq. REFLECT_ODD) then
       do i = lo(1)-ng, hi(1)+ng
          do j = 1, ng
             s(i,lo(2)-j) = -s(i,lo(2)+j-1)
          end do
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(2,1) = ',bc(2,1)
       call bl_error('BC(2,1) = NOT YET SUPPORTED')
    end if

    if (bc(2,2) .eq. EXT_DIR) then
       call bl_error('BC(1,1) = EXT_DIR NOT YET SUPPORTED')
    else if (bc(2,2) .eq. FOEXTRAP) then
       do i = lo(1)-ng, hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = s(i,hi(2))
       end do
    else if (bc(2,2) .eq. HOEXTRAP) then
       do i = lo(1)-ng, hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = EIGHTH* &
               (15.0_dp_t*s(i,hi(2)) - 10.0_dp_t*s(i,hi(2)-1) + 3.0_dp_t*s(i,hi(2)-2))
       end do
    else if (bc(2,2) .eq. REFLECT_EVEN) then
       do i = lo(1)-ng, hi(1)+ng
          do j = 1, ng
             s(i,hi(2)+j) = s(i,hi(2)-j+1)
          end do
       end do
    else if (bc(2,2) .eq. REFLECT_ODD) then
       do i = lo(1)-ng, hi(1)+ng
          do j = 1, ng
             s(i,hi(2)+j) = -s(i,hi(2)-j+1)
          end do
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(2,2) = ',bc(2,2)
       call bl_error('BC(2,2) = NOT YET SUPPORTED')
    end if

  end subroutine setbc_2d

  subroutine setbc_3d(s,lo,ng,bc)

    use bl_constants_module
    use bc_module

    integer        , intent(in   ) :: lo(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:, lo(2)-ng:, lo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)

    ! Local variables
    integer :: i,j,k,hi(3)

    hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)
    hi(3) = lo(3) + size(s,dim=3) - (2*ng+1)

    if (bc(1,1) .eq. EXT_DIR) then

    else if (bc(1,2) .eq. FOEXTRAP) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             s(hi(1)+1:hi(1)+ng,j,k) = s(hi(1),j,k)
          end do
       end do
    else if (bc(1,2) .eq. HOEXTRAP) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             s(hi(1)+1:hi(1)+ng,j,k) = &
                  ( 15.0_dp_t * s(hi(1)  ,j,k) &
                  -10.0_dp_t * s(hi(1)-1,j,k) &
                  + 3.0_dp_t * s(hi(1)-2,j,k) ) * EIGHTH
          end do
       end do
    else if (bc(1,2) .eq. REFLECT_EVEN) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             do i = 1,ng
                s(hi(1)+i,j,k) = s(hi(1)-i+1,j,k)
             end do
          end do
       end do
    else if (bc(1,2) .eq. REFLECT_ODD) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             do i = 1,ng
                s(hi(1)+i,j,k) = -s(hi(1)-i+1,j,k)
             end do
          end do
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(1,2) = ',bc(1,2)
       call bl_error('BC(1,2) = NOT YET SUPPORTED')
    end if

    if (bc(2,1) .eq. EXT_DIR) then

    else if (bc(2,1) .eq. FOEXTRAP .or. bc(2,1) .eq. REFLECT_EVEN) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,lo(2)-ng:lo(2)-1,k) = s(i,lo(2),k)
          end do
       end do
    else if (bc(2,1) .eq. HOEXTRAP) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,lo(2)-ng:lo(2)-1,k) = &
                  ( 15.0_dp_t * s(i,lo(2)  ,k) &
                  -10.0_dp_t * s(i,lo(2)+1,k) &
                  + 3.0_dp_t * s(i,lo(2)+2,k) ) * EIGHTH
          end do
       end do
    else if (bc(2,1) .eq. REFLECT_EVEN) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             do j = 1,ng
                s(i,lo(2)-j,k) = s(i,lo(2)+j-1,k)
             end do
          end do
       end do
    else if (bc(2,1) .eq. REFLECT_ODD) then
       do k = lo(3)-1,hi(3)+1
          do i = lo(1)-ng,hi(1)+ng
             do j = 1,ng
                s(i,lo(2)-j,k) = -s(i,lo(2)+j-1,k)
             end do
          end do
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(2,1) = ',bc(2,1)
       call bl_error('BC(2,1) = NOT YET SUPPORTED')
    end if

    if (bc(2,2) .eq. EXT_DIR) then
    else if (bc(2,2) .eq. FOEXTRAP .or. bc(2,2) .eq. REFLECT_EVEN) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,hi(2)+1:hi(2)+ng,k) = s(i,hi(2),k)
          end do
       end do
    else if (bc(2,2) .eq. HOEXTRAP) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,hi(2)+1:hi(2)+ng,k) = &
                  ( 15.0_dp_t * s(i,hi(2)  ,k) &
                  -10.0_dp_t * s(i,hi(2)-1,k) &
                  + 3.0_dp_t * s(i,hi(2)-2,k) ) * EIGHTH
          end do
       end do
    else if (bc(2,2) .eq. REFLECT_EVEN) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             do j = 1,ng
                s(i,hi(2)+j,k) = s(i,hi(2)-j+1,k)
             end do
          end do
       end do
    else if (bc(2,2) .eq. REFLECT_ODD) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             do j = 1,ng
                s(i,hi(2)+j,k) = -s(i,hi(2)-j+1,k)
             end do
          end do
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(2,2) = ',bc(2,2)
       call bl_error('BC(2,2) = NOT YET SUPPORTED')
    end if

    if (bc(3,1) .eq. EXT_DIR) then
    else if (bc(3,1) .eq. FOEXTRAP .or. bc(3,1) .eq. REFLECT_EVEN) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,lo(3)-ng:lo(3)-1) = s(i,j,lo(3))
          end do
       end do
    else if (bc(3,1) .eq. HOEXTRAP) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,lo(3)-ng:lo(3)-1) = &
                  ( 15.0_dp_t * s(i,j,lo(3)  ) &
                  -10.0_dp_t * s(i,j,lo(3)+1) &
                  + 3.0_dp_t * s(i,j,lo(3)+2) ) * EIGHTH
          end do
       end do
    else if (bc(3,1) .eq. REFLECT_EVEN) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             do k = 1,ng
                s(i,j,lo(3)-k) = s(i,j,lo(3)+k-1)
             end do
          end do
       end do
    else if (bc(3,1) .eq. REFLECT_ODD) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             do k = 1,ng
                s(i,j,lo(3)-k) = -s(i,j,lo(3)+k-1)
             end do
          end do
       end do
    else if (bc(3,1) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(3,1) = ',bc(3,1)
       call bl_error('BC(3,1) = NOT YET SUPPORTED')
    end if

    if (bc(3,2) .eq. EXT_DIR) then
    else if (bc(3,2) .eq. FOEXTRAP .or. bc(3,2) .eq. REFLECT_EVEN) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,hi(3)+1:hi(3)+ng) = s(i,j,hi(3))
          end do
       end do
    else if (bc(3,2) .eq. HOEXTRAP) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,hi(3)+1:hi(3)+ng) = &
                  ( 15.0_dp_t * s(i,j,hi(3)  ) &
                  -10.0_dp_t * s(i,j,hi(3)-1) &
                  + 3.0_dp_t * s(i,j,hi(3)-2) ) * EIGHTH
          end do
       end do
    else if (bc(3,2) .eq. REFLECT_EVEN) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             do k = 1,ng
                s(i,j,hi(3)+k) = s(i,j,hi(3)-k+1)
             end do
          end do
       end do
    else if (bc(3,2) .eq. REFLECT_ODD) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             do k = 1,ng
                s(i,j,hi(3)+k) = -s(i,j,hi(3)-k+1)
             end do
          end do
       end do
    else if (bc(3,2) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(3,2) = ',bc(3,2)
       call bl_error('BC(3,2) = NOT YET SUPPORTED')
    end if

  end subroutine setbc_3d

end module setbc_module
