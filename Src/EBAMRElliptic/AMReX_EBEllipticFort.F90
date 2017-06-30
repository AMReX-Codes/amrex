#include "AMReX_CONSTANTS.H"

module ebefnd_module

  !     since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
  !     for e.g., #if (BL_SPACEDIM == 1) statements.

  implicit none

  public

contains
  
  integer function identitymat(i, j)
    implicit none
    integer i, j, retval
    retval = 0
    if (i.eq.j) then
       retval = 1
    endif

    return retval
  end function identitymat

    
  subroutine ebefnd_decrinvrelcoefebco( &
       relco, relco_lo, relco_hi,  relco_nco, &
       bcoef, bcoef_lo, bcoef_hi,  bcoef_nco, &
       gridlo, gridhi, &
       beta, dx, idir)
       bind(C, name="ebefnd_decrinvrelcoefebco")

    use amrex_fort_module, only : amrex_spacedim, c_real=>amrex_real

    implicit none

    integer      :: i, j, k, relco_nco, bcoef_nco, idir
    integer      :: relco_lo(0:2),relco_hi(0:2)
    integer      :: bcoef_lo(0:2),bcoef_hi(0:2)
    integer      :: gridlo(0:2), gridhi(0:2), ii, jj, kk
    real(c_real) :: beta, dx
    real(c_real) :: bcoef(bcoef_lo(0):bcoef_hi(0),bcoef_lo(1):bcoef_hi(1),bcoef_lo(2):bcoef_hi(2), 0:bcoef_nco-1)
    real(c_real) :: relco(relco_lo(0):relco_hi(0),relco_lo(1):relco_hi(1),relco_lo(2):relco_hi(2), 0:relco_nco-1)

    ii = identitymat(idir, 0)
    jj = identitymat(idir, 1)
    kk = identitymat(idir, 2)

    do k = gridlo(2), gridhi(2)
       do j = gridlo(1), gridhi(1)
          do i = gridlo(0), gridhi(0)

             relco(i,j,k,0) = relco(i,j,k,0) &
                  - beta*( &
                  bcoef(i+ii,j+jj,k+kk,0) + &
                  bcoef(i   ,j   ,k   ,0))/(dx*dx)

          enddo
       enddo
    enddo

    return 
  end subroutine ebefnd_decrinvrelcoefebco


  subroutine ebefnd_invertlambdaebco( &
       relco, relco_lo, relco_hi,  relco_nco, &
       gridlo, gridhi, safety)
       bind(C, name="ebefnd_invertlambdaebco")

    use amrex_fort_module, only : amrex_spacedim, c_real=>amrex_real

    implicit none

    integer      :: i, j, k, relco_nco,  idir
    integer      :: relco_lo(0:2),relco_hi(0:2)
    integer      :: gridlo(0:2), gridhi(0:2)
    real(c_real) :: safety


    do k = gridlo(2), gridhi(2)
       do j = gridlo(1), gridhi(1)
          do i = gridlo(0), gridhi(0)

             relco(i,j,k,0) = safety/relco(i,j,k,0)

          enddo
       enddo
    enddo

    return 
  end subroutine ebefnd_decrinvrelcoefebco

end module ebefnd_module

