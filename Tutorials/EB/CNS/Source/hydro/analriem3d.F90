module riemann_module

  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private

  public :: analriem

contains

  subroutine analriem (gamma, smallp, smallr, lo, hi, j, k, &
       rl, ul, pl, ut1l, ut2l, rr, ur, pr, ut1r, ut2r, flux, flo, fhi, iu, iut1, iut2)
    integer, intent(in) :: lo, hi, j, k, flo(3), fhi(3), iu, iut1, iut2
    real(rt), intent(in) :: gamma, smallp, smallr
    real(rt), dimension(lo:hi), intent(in) :: rl, ul, pl, ut1l, ut2l
    real(rt), dimension(lo:hi), intent(in) :: rr, ur, pr, ut1r, ut2r
    real(rt), intent(inout) :: flux(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3),5)

    real(rt) :: clsql,clsqr
    real(rt) :: wl,wr,wlsq,wrsq
    real(rt) :: dpditer,zp,zm,denom
    real(rt) :: ustnm1,ustnp1
    real(rt) :: uo,ro,rstar
    real(rt) :: po,pgdnv,rgdnv,ugdnv,dpjmp,cstar
    real(rt) :: spin,spout,sgnm,co,wo, wosq
    real(rt) :: frac
    real(rt) :: utrans1,utrans2
    integer :: i, iter
    real(rt), dimension(lo:hi) :: cleft, cright, pstar, pstnm1, ustarp, ustarm, ustar
    real(rt), parameter :: weakwv = 1.d-03
    real(rt), parameter :: small = 1.d-6
    integer, parameter :: itno = 3

    do i = lo, hi
       clsql = gamma*pl(i)*rl(i)
       clsqr = gamma*pr(i)*rr(i)
       wl = sqrt(clsql)
       wr = sqrt(clsqr)
       cleft(i) = wl/rl(i)
       cright(i) = wr/rr(i)
    
       pstar(i) = (wl*pr(i) + wr*pl(i) - wr*wl*(ur(i)-ul(i)))/(wl+wr)
       pstar(i) = max(pstar(i),smallp)
       pstnm1(i) = pstar(i)

       wlsq = (.5d0*(gamma-1.d0)*(pstar(i)+pl(i))+pstar(i))*rl(i)
       wrsq = (.5d0*(gamma-1.d0)*(pstar(i)+pr(i))+pstar(i))*rr(i)

       wl = sqrt(wlsq)
       wr = sqrt(wrsq)
       ustarp(i) = ul(i) - (pstar(i)-pl(i))/wl
       ustarm(i) = ur(i) + (pstar(i)-pr(i))/wr

       pstar(i) = (wl*pr(i) + wr*pl(i) - wr*wl*(ur(i)-ul(i)))/(wl+wr)
       pstar(i) = max(pstar(i),smallp)
    end do

    do iter = 1, itno
       do i = lo, hi  ! what changes in the iteration: ustarm, ustarp, pstnm1, pstar, ustar
          
          wlsq = (.5d0*(gamma-1.d0)*(pstar(i)+pl(i))+pstar(i))*rl(i)
          wrsq = (.5d0*(gamma-1.d0)*(pstar(i)+pr(i))+pstar(i))*rr(i)

          wl = 1.d0/sqrt(wlsq)
          wr = 1.d0/sqrt(wrsq)

          ustnm1 = ustarm(i)
          ustnp1 = ustarp(i)

          ustarm(i) = ur(i) - (pr(i) - pstar(i))*wr
          ustarp(i) = ul(i) + (pl(i) - pstar(i))*wl

          dpditer = abs(pstnm1(i)-pstar(i))
          zp = abs(ustarp(i)-ustnp1)
          if (zp-weakwv*cleft(i) .lt. 0.0d0 ) then
             zp = dpditer*wl
          endif
          zm = abs(ustarm(i)-ustnm1)
          if (zm-weakwv*cright(i) .lt. 0.0d0 ) then
             zm = dpditer*wr
          endif

          denom = dpditer/max(zp+zm,small*(cleft(i)+cright(i)))
          pstnm1(i) = pstar(i)
          pstar(i) = pstar(i) - denom*(ustarm(i)-ustarp(i))
          pstar(i) = max(pstar(i),smallp)
          ustar(i) = 0.5d0*(ustarm(i)+ustarp(i))

       end do
    end do

    do i = lo, hi
       if(ustar(i) .gt. 0d0)then
          ro = rl(i)
          uo = ul(i)
          po = pl(i)
          sgnm = 1.d0
          utrans1 = ut1l(i)
          utrans2 = ut2l(i)
       elseif (ustar(i) .lt. 0.d0)then
          ro = rr(i)
          uo = ur(i)
          po = pr(i)
          sgnm = -1.d0
          utrans1 = ut1r(i)
          utrans2 = ut2r(i)
       else
          uo = 0.5d0*(ur(i)+ul(i))
          po = 0.5d0*(pr(i)+pl(i))
          ro = 2.d0*(rl(i)*rr(i))/(rl(i)+rr(i))
          sgnm = 1.d0
          utrans1 = 0.5d0*(ut1l(i)+ut1r(i))
          utrans2 = 0.5d0*(ut2l(i)+ut2r(i))
       endif
       wosq = (.5d0*(gamma-1.d0)*(pstar(i)+po)+pstar(i))*ro
       co = sqrt(gamma * po / ro)
       wo = sqrt(wosq)
       dpjmp = pstar(i)-po
       rstar = ro/(1.d0-ro*dpjmp/wosq)
       cstar = sqrt(gamma * pstar(i) / rstar)
       spout = co-sgnm*uo
       spin = cstar - sgnm*uo
       if(pstar(i) .ge. po)then
          spin = wo/ro-sgnm*uo
          spout = spin
       endif
       frac = 0.5d0*(1.d0+(spin+spout)/max(spout-spin,spin+spout,small*(cleft(i)+cright(i))))

       if(spout .lt. 0.d0)then
          rgdnv = ro 
          ugdnv = uo
          pgdnv = po
       elseif(spin .ge. 0.d0)then
          rgdnv = rstar
          ugdnv = ustar(i)
          pgdnv = pstar(i)
       else
          rgdnv = frac*rstar    + (1.d0 - frac)* ro
          ugdnv = frac*ustar(i) + (1.d0 - frac)* uo
          pgdnv = frac*pstar(i) + (1.d0 - frac)* po
       endif
    
       flux(i,j,k,1) = rgdnv*ugdnv
       flux(i,j,k,iu) = rgdnv*ugdnv*ugdnv+pgdnv
       flux(i,j,k,iut1) = rgdnv*ugdnv*utrans1
       flux(i,j,k,iut2) = rgdnv*ugdnv*utrans2
       flux(i,j,k,5) = ugdnv*(0.5d0*rgdnv*(ugdnv*ugdnv+utrans1*utrans1+utrans2*utrans2) &
            &        + pgdnv/(gamma -1.d0) + pgdnv)
    end do

  end subroutine analriem

end module riemann_module
