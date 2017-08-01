module nbrsTest_nd_module

  use amrex_fort_module, only : amrex_real, dim=>bl_spacedim
  use nbrs_test_module, only : nbr_sten

  implicit none

contains

  subroutine fill_redist_stencil(lo, hi, sten, Nsten, mask, mask_lo, mask_hi, &
       vf, vf_lo, vf_hi) bind(C,name="fill_redist_stencil")

    integer,          intent(in) :: lo(0:2), hi(0:2)
    integer,          intent(in) :: Nsten
    type(nbr_sten),intent(inout) :: sten(0:Nsten-1)
    integer,          intent(in) :: mask_lo(3), mask_hi(3)
    integer,          intent(in) :: vf_lo(3), vf_hi(3)
    integer,          intent(in) :: mask(mask_lo(1):mask_hi(1),mask_lo(2):mask_hi(2),mask_lo(3):mask_hi(3))
    real(amrex_real), intent(in) :: vf(vf_lo(1):vf_hi(1),vf_lo(2):vf_hi(2),vf_lo(3):vf_hi(3))
    integer :: i,j,k,n,ii,jj,kk,jjlo,jjhi,kklo,kkhi,iii,jjj,kkk
    real(amrex_real) :: kappa_tot

    ! Note: mask is such that 1: reg, 0:sv, -1: covered, 2: mv, -2: outside domain, -3: crse

    if (dim .ge. 2) then
       jjlo = -1
       jjhi = +1
    else
       jjlo = 0
       jjhi = 0
    endif

    if (dim .eq. 3) then
       kklo = -1
       kkhi = +1
    else
       kklo = 0
       kkhi = 0
    endif

    do n = 0, Nsten-1
       i = sten(n) % iv(0)
       j = sten(n) % iv(1)
       k = sten(n) % iv(2)
       if (i.ge.lo(0) .and. i.le.hi(0) &
            .and. (dim.lt.2 .or. (j.ge.lo(1) .and. j.le.hi(1)))&
            .and. (dim.lt.3 .or. (k.ge.lo(2) .and. k.le.hi(2))) ) then
          kappa_tot = 0.d0
          do kk=kklo,kkhi
             kkk = k+kk
             do jj=jjlo,jjhi
                jjj = j+jj
                do ii=-1,1
                   iii = i+ii
                   if (mask(iii,jjj,kkk) .eq. 1 &
                        .or. mask(iii,jjj,kkk) .eq. 0) then
                      kappa_tot = kappa_tot + vf(iii,jjj,kkk)
                   endif
                enddo
             enddo
          enddo
          kappa_tot = 1.d0 / kappa_tot
          do kk=kklo,kkhi
             kkk = k+kk
             do jj=jjlo,jjhi
                jjj = j+jj
                do ii=-1,1
                   iii = i+ii
                   if (mask(iii,jjj,kkk) .eq. 1 &
                        .or. mask(iii,jjj,kkk) .eq. 0) then
                      sten(n) % val(ii,jj,kk) = vf(iii,jjj,kkk) * kappa_tot
                   else
                      sten(n) % val(ii,jj,kk) = 0.d0
                   endif
                enddo
             enddo
          enddo
       endif
    enddo

  end subroutine fill_redist_stencil

  subroutine apply_redist_stencil(lo, hi, sten, Nsten, vin, vin_lo, vin_hi, &
    vout, vout_lo, vout_hi) bind(C,name="apply_redist_stencil")

    integer,          intent(in   ) :: lo(0:2), hi(0:2)
    integer,          intent(in   ) :: Nsten
    type(nbr_sten),   intent(in   ) :: sten(0:Nsten-1)
    integer,          intent(in   ) ::  vin_lo(3),  vin_hi(3)
    integer,          intent(in   ) :: vout_lo(3), vout_hi(3)
    real(amrex_real), intent(in   ) ::  vin( vin_lo(1):vin_hi(1),  vin_lo(2):vin_hi(2),  vin_lo(3):vin_hi(3))
    real(amrex_real), intent(inout) :: vout(vout_lo(1):vout_hi(1),vout_lo(2):vout_hi(2),vout_lo(3):vout_hi(3))
    integer :: i,j,k,n,ii,jj,kk,jjlo,jjhi,kklo,kkhi,iii,jjj,kkk
    real(amrex_real) :: vtot

    if (dim .ge. 2) then
       jjlo = -1
       jjhi = +1
    else
       jjlo = 0
       jjhi = 0
    endif

    if (dim .eq. 3) then
       kklo = -1
       kkhi = +1
    else
       kklo = 0
       kkhi = 0
    endif

    do n = 0, Nsten-1
       i = sten(n) % iv(0)
       j = sten(n) % iv(1)
       k = sten(n) % iv(2)
       if (i.ge.lo(0) .and. i.le.hi(0) &
            .and. (dim.lt.2 .or. (j.ge.lo(1) .and. j.le.hi(1)))&
            .and. (dim.lt.3 .or. (k.ge.lo(2) .and. k.le.hi(2))) ) then
          vtot = 0.d0
          do kk=kklo,kkhi
             kkk = k+kk
             do jj=jjlo,jjhi
                jjj = j+jj
                do ii=-1,1
                   iii = i+ii
                   vtot = vtot + sten(n)%val(ii,jj,kk) * vin(iii,jjj,kkk)
                enddo
             enddo
          enddo
          vout(i,j,k) = vtot
       endif
    enddo
    
  end subroutine apply_redist_stencil

end module nbrsTest_nd_module
