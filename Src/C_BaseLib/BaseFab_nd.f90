module basefab_nd_module

  use bl_fort_module, only : c_real

  implicit none

contains

  ! This function adds scaled floating point numbers from one array to another.
  subroutine fort_saxpy(lo, hi, dst, dlo, dhi, a, src, slo, shi, sblo, ncomp) bind(c,name='fort_saxpy')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp
    real(c_real), intent(in   ) :: a
    real(c_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(c_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    
    integer :: i,j,k,n,off(3)

    off = lo - sblo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = dst(i,j,k,n) + a * src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do

  end subroutine fort_saxpy

end module basefab_nd_module
