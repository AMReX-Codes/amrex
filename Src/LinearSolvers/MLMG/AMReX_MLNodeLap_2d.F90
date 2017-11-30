module amrex_mlnodelap_2d_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlndlap_sigma_cctoedge

contains

  subroutine amrex_mlndlap_sigma_cctoedge (xlo, xhi, ylo, yhi, sigx, sxlo, sxhi, sigy, sylo, syhi, &
       sigcc, clo, chi) bind(c, name='amrex_mlndlap_sigma_cctoedge')
    integer, dimension(2), intent(in) :: xlo, xhi, ylo, yhi, sxlo, sxhi, sylo, syhi, clo, chi
    real(amrex_real), intent(inout) :: sigx (sxlo(1):sxhi(1),sxlo(2):sxhi(2))
    real(amrex_real), intent(inout) :: sigy (sylo(1):syhi(1),sylo(2):syhi(2))
    real(amrex_real), intent(in   ) :: sigcc( clo(1): chi(1), clo(2): chi(2))

    integer :: i, j

    do    j = xlo(2), xhi(2)
       do i = xlo(1), xhi(1)
          sigx(i,j) = 0.5d0*(sigcc(i,j-1)+sigcc(i,j))
       end do
    end do

    do    j = ylo(2), yhi(2)
       do i = ylo(1), yhi(1)
          sigy(i,j) = 0.5d0*(sigcc(i-1,j)+sigcc(i,j))
       end do
    end do

  end subroutine amrex_mlndlap_sigma_cctoedge

end module amrex_mlnodelap_2d_module
