
module amrex_ebcellflag_module
  implicit none
  private
  public :: is_regular, is_single_valued, is_multi_valued, is_covered, get_neighbors
  
  integer, parameter :: w_type      = 2
  integer, parameter :: w_numvofs   = 3
  integer, parameter :: pos_numvofs = 2;
  integer, parameter :: pos_ngbr(-1:1,-1:1,-1:1) = &
       reshape((/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26/),&
               [3,3,3]) + (w_type+w_numvofs)

  integer, parameter :: regular       = 0
  integer, parameter :: single_valued = 1
  integer, parameter :: multi_valued  = 2
  integer, parameter :: covered       = 3

contains
  
  pure logical function is_regular (flag)
    integer, intent(in) :: flag
    is_regular = ibits(flag,0,w_type) .eq. regular
  end function is_regular

  pure logical function is_single_valued (flag)
    integer, intent(in) :: flag
    is_single_valued = ibits(flag,0,w_type) .eq. single_valued
  end function is_single_valued

  pure logical function is_multi_valued (flag)
    integer, intent(in) :: flag
    is_multi_valued = ibits(flag,0,w_type) .eq. multi_valued
  end function is_multi_valued

  pure logical function is_covered (flag)
    integer, intent(in) :: flag
    is_covered = ibits(flag,0,w_type) .eq. covered
  end function is_covered

  pure subroutine get_neighbors (flag, ngbr)
    integer, intent(in) :: flag
    integer, intent(inout) :: ngbr(-1:1,-1:1,-1:1)
    integer :: i, j, k
    do k = -1,1
       do j = -1,1
          do i = -1,1
             if (btest(flag,pos_ngbr(i,j,k))) then
                ngbr(i,j,k) = 1
             else
                ngbr(i,j,k) = 0
             end if
          end do
       end do
    end do
  end subroutine get_neighbors

end module amrex_ebcellflag_module
