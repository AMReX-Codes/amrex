module amrex_mlnodelap_3d_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlndlap_sigma_cctoedge

contains

  subroutine amrex_mlndlap_sigma_cctoedge (xlo, xhi, ylo, yhi, zlo, zhi, &
       sigx, sxlo, sxhi, sigy, sylo, syhi, sigz, szlo, szhi, &
       sigcc, clo, chi) bind(c, name='amrex_mlndlap_sigma_cctoedge')
    integer, dimension(3), intent(in) :: xlo, xhi, ylo, yhi, zlo, zhi, sxlo, sxhi, sylo, syhi, szlo, szhi, clo, chi
    real(amrex_real), intent(inout) :: sigx (sxlo(1):sxhi(1),sxlo(2):sxhi(2),sxlo(3):sxhi(3))
    real(amrex_real), intent(inout) :: sigy (sylo(1):syhi(1),sylo(2):syhi(2),sylo(3):syhi(3))
    real(amrex_real), intent(inout) :: sigz (szlo(1):szhi(1),szlo(2):szhi(2),szlo(3):szhi(3))
    real(amrex_real), intent(in   ) :: sigcc( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3))

    integer :: i, j, k

    do       k = xlo(3), xhi(3)
       do    j = xlo(2), xhi(2)
          do i = xlo(1), xhi(1)
             sigx(i,j,k) = 0.25d0*(sigcc(i,j-1,k-1) &
                  &              + sigcc(i,j  ,k-1) &
                  &              + sigcc(i,j-1,k  ) &
                  &              + sigcc(i,j  ,k  ))
          end do
       end do
    end do

    do       k = ylo(3), yhi(3)
       do    j = ylo(2), yhi(2)
          do i = ylo(1), yhi(1)
             sigy(i,j,k) = 0.25d0*(sigcc(i-1,j,k-1) &
                  &              + sigcc(i  ,j,k-1) &
                  &              + sigcc(i-1,j,k  ) &
                  &              + sigcc(i  ,j,k  ))
          end do
       end do
    end do

    do       k = zlo(3), zhi(3)
       do    j = zlo(2), zhi(2)
          do i = zlo(1), zhi(1)
             sigz(i,j,k) = 0.25d0*(sigcc(i-1,j-1,k) &
                  &              + sigcc(i  ,j-1,k) &
                  &              + sigcc(i-1,j  ,k) &
                  &              + sigcc(i  ,j  ,k))
          end do
       end do
    end do

  end subroutine amrex_mlndlap_sigma_cctoedge

end module amrex_mlnodelap_3d_module
