module lapl_simple_module

! since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
! for e.g., #if (BL_SPACEDIM == 1) statements.

  implicit none

  public

contains

  subroutine lapl_simple(lph, lph_lo, lph_hi, phi,phi_lo,phi_hi,reglo,reghi,dx) &
       bind(C, name="fort_lapl_simple")

    use amrex_fort_module, only : amrex_spacedim, c_real=>amrex_real

    implicit none

    integer      :: phi_lo(3),phi_hi(3),i,j,k
    integer      :: lph_lo(3),lph_hi(3)
    integer      :: reglo(3), reghi(3)
    real(c_real) :: dx
    real(c_real) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    real(c_real) :: lph(lph_lo(1):lph_hi(1),lph_lo(2):lph_hi(2),lph_lo(3):lph_hi(3))

#if (BL_SPACEDIM == 3)
      do k = reglo(3), reghi(3)
         do j = reglo(2), reghi(2)
            do i = reglo(1), reghi(1)
               lph(i,j,k) = &
               (phi(i+1,j  ,k  ) + phi(i-1,j  ,k  ) - 2.d0*phi(i,j,k))/dx/dx + &
               (phi(i  ,j+1,k  ) + phi(i  ,j-1,k  ) - 2.d0*phi(i,j,k))/dx/dx + &
               (phi(i  ,j  ,k+1) + phi(i  ,j  ,k-1) - 2.d0*phi(i,j,k))/dx/dx  
            enddo
         enddo
      enddo
#else
      do k = reglo(3), reghi(3)
         do j = reglo(2), reghi(2)
            do i = reglo(1), reghi(1)
               lph(i,j,k) = &
               (phi(i+1,j  ,k  ) + phi(i-1,j  ,k  ) - 2.d0*phi(i,j,k))/dx/dx + &
               (phi(i  ,j+1,k  ) + phi(i  ,j-1,k  ) - 2.d0*phi(i,j,k))/dx/dx 
            enddo
         enddo
      enddo

#endif

  end subroutine lapl_simple
  
end module lapl_simple_module
