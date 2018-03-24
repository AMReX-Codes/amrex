module amrex_mlnodelap_nd_module

  use amrex_fort_module, only : amrex_real
  implicit none

  ! external dirichlet at physical boundary or internal dirichlet at crse/fine boundary
  integer, parameter :: dirichlet = 1

  integer, parameter :: fine_node = 2

  private
  public :: amrex_mlndlap_copy_fine_node, amrex_mlndlap_zero_dirichlet_node

contains

  subroutine amrex_mlndlap_copy_fine_node (lo, hi, d, dlo, dhi, s, slo, shi, m, mlo, mhi) &
       bind(c, name='amrex_mlndlap_copy_fine_node')
    integer, dimension(3), intent(in) :: lo, hi, dlo, dhi, slo, shi, mlo, mhi
    real(amrex_real), intent(inout) :: d(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
    real(amrex_real), intent(in   ) :: s(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    integer         , intent(in   ) :: m(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (m(i,j,k) .eq. fine_node) d(i,j,k) = s(i,j,k)
          end do
       end do
    end do
  end subroutine amrex_mlndlap_copy_fine_node

  subroutine amrex_mlndlap_zero_dirichlet_node (lo, hi, d, dlo, dhi, m, mlo, mhi) &
       bind(c, name='amrex_mlndlap_zero_dirichlet_node')
    integer, dimension(3), intent(in) :: lo, hi, dlo, dhi, mlo, mhi
    real(amrex_real), intent(inout) :: d(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
    integer         , intent(in   ) :: m(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (m(i,j,k) .eq. dirichlet) d(i,j,k) = 0.0
          end do
       end do
    end do    
  end subroutine amrex_mlndlap_zero_dirichlet_node

end module amrex_mlnodelap_nd_module
