module poisson_module

  use amrex_fort_module, only : rt=>amrex_real
  
  implicit none

  private
  
contains

  subroutine init_data (lo,hi,q,qlo,qhi) bind(c,name='init_data')

    integer , intent(in   ) :: lo(3), hi(3), qlo(3), qhi(3)
    real(rt), intent(inout) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))

    integer :: i,j,k

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if (i .eq. 70 .and. j .eq. 70) then
                q(i,j,k) = 1.d0
             else
                q(i,j,k) = 0.d0                
             end if
          end do
       end do
    end do
    
  end subroutine init_data

end module poisson_module
