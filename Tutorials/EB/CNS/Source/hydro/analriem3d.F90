module riemann_module

  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private

  public :: analriem, hllriem

contains

  pure real(rt) function ff (a,b,c,d,e,f)
    real(rt), intent(in) :: a,b,c,d,e,f
    ff = (a*f + c*d - c*f*(e-b))/(c+f)
  end function ff

  subroutine analriem(gamma,statel,stater,smallp,smallr,flux,  debug)
    logical, intent(in), optional :: debug
    real(rt), intent(in) :: statel(5),stater(5), smallp, smallr, gamma
    real(rt), intent(out) :: flux(5)
    
    real(rt) :: pl,rl,ul,pr,rr,ur
    real(rt) :: ut1l,ut2l,ut1r,ut2r
    real(rt) :: clsql,clsqr,vbarl,vbarr
    real(rt) :: pstar,pstnm1,ustarp,ustarm,ustar
    real(rt) :: wl,wr,wlsq,wrsq
    real(rt) :: cleft,cright
    real(rt) :: dpditer,zp,zm,denom
    real(rt) :: ustnm1,ustnp1
    real(rt) :: uo,ro,rstar
    real(rt) :: po,pgdnv,rgdnv,ugdnv,dpjmp,cstar
    real(rt) :: spin,spout,sgnm,co,tauo,ushock,wo, wosq
    real(rt) :: frac
    real(rt) :: utrans1,utrans2
    integer :: itno,iter
    real(rt), parameter :: weakwv = 1.d-03
    real(rt), parameter :: small = 1.d-6

    rl = statel(1)
    ul = statel(2)
    pl = statel(3)
    ut1l = statel(4)
    ut2l = statel(5)
    
    rr = stater(1)
    ur = stater(2)
    pr = stater(3)
    ut1r = stater(4)
    ut2r = stater(5)
    
    pl = max(pl,smallp)
    pr = max(pr,smallp)
    rl = max(rl,smallr)
    rr = max(rr,smallr)
    clsql = gamma*pl*rl
    clsqr = gamma*pr*rr
    vbarl = 1.d0/rl
    vbarr = 1.d0/rr
    wl = sqrt(clsql)
    wr = sqrt(clsqr)
    cleft = wl/rl
    cright = wr/rr
    
    
    pstar=(wl*pr+wr*pl-wr*wl*(ur-ul))/(wl+wr)
    pstar=max(pstar,smallp)
    pstnm1 = pstar
    
    wlsq = (.5d0*(gamma-1.d0)*(pstar+pl)+pstar)/vbarl
    wrsq = (.5d0*(gamma-1.d0)*(pstar+pr)+pstar)/vbarr
    
    wl = sqrt(wlsq)
    wr = sqrt(wrsq)
    ustarp=ul-(pstar-pl)/wl
    ustarm=ur+(pstar-pr)/wr
    
    pstar = ff(pl,ul,wl,pr,ur,wr)
    pstar=max(pstar,smallp)
    
    itno = 3

    do iter = 1,itno

       wlsq = (.5d0*(gamma-1.d0)*(pstar+pl)+pstar)/vbarl
       wrsq = (.5d0*(gamma-1.d0)*(pstar+pr)+pstar)/vbarr
       
       wl = 1.d0/sqrt(wlsq)
       wr = 1.d0/sqrt(wrsq)
       
       ustnm1=ustarm
       ustnp1=ustarp
       
       ustarm=ur-(pr-pstar)*wr
       ustarp=ul+(pl-pstar)*wl
       
       dpditer=abs(pstnm1-pstar)
       zp=abs(ustarp-ustnp1)
       if (zp-weakwv*cleft .lt. 0.0 ) then
          zp = dpditer*wl
       endif
       zm=abs(ustarm-ustnm1)
       if (zm-weakwv*cright .lt. 0.0 ) then
          zm = dpditer*wr
       endif
       
       denom=dpditer/max(zp+zm,small*(cleft+cright))
       pstnm1 = pstar
       pstar=pstar-denom*(ustarm-ustarp)
       pstar=max(pstar,smallp)
       ustar = 0.5d0*(ustarm+ustarp)

    end do

    if (present(debug)) then
       if(debug)then
          write(6,*)"ustart,pstar", ustar,pstar
          write(6,*)" here"
          write(6,*)"lef", pl,rl,ul
          write(6,*)"rig", pr,rr,ur
       end if
    endif

    if(ustar .gt. 0d0)then
       ro = rl
       uo = ul
       po = pl
       tauo = vbarl
       sgnm = 1.d0
       utrans1 = ut1l
       utrans2 = ut2l
    elseif (ustar.le.0.d0)then
       ro = rr
       uo = ur
       po = pr
       tauo = vbarr
       sgnm =- 1.d0
       utrans1 = ut1r
       utrans2 = ut2r
    else
       uo = 0.5d0*(ur+ul)
       po = 0.5d0*(pr+pl)
       tauo = 0.5d0*(vbarl+vbarr)
       ro = 1.d0/tauo
       sgnm = 1.d0
       utrans1 = 0.5d0*(ut1l+ut1r)
       utrans2 = 0.5d0*(ut2l+ut2r)
    endif
    wosq = (.5d0*(gamma-1.d0)*(pstar+po)+pstar)/tauo
    co = sqrt(gamma * po / ro)
    wo = sqrt(wosq)
    dpjmp = pstar-po
    rstar = ro/(1.d0-ro*dpjmp/wosq)
    cstar = sqrt(gamma * pstar / rstar)
    spout = co-sgnm*uo
    spin = cstar - sgnm*uo
    ushock = wo/ro-sgnm*uo
    if(pstar .ge. po)then
       spin = ushock
       spout = ushock
    endif
    frac = 0.5d0*(1.d0+(spin+spout)/max(spout-spin,spin+spout,small*(cleft+cright)))

    if (present(debug)) then
       if(debug)then
          write(6,*)"cstar,sgnm,uo,co",cstar,sgnm,uo,co
          write(6,*)"dpjmp,spin,spout,frac", dpjmp,spin,spout,frac
       endif
    end if
      
    if(spout .lt. 0.d0)then
       rgdnv = ro 
       ugdnv = uo
       pgdnv = po
       
    elseif(spin .ge. 0.d0)then
       rgdnv = rstar 
       ugdnv = ustar
       pgdnv = pstar
       
    else
       rgdnv = frac*rstar + (1.d0 - frac)* ro
       ugdnv = frac*ustar + (1.d0 - frac)* uo
       pgdnv = frac*pstar + (1.d0 - frac)* po
       
    endif
    
    flux(1) = rgdnv*ugdnv
    flux(2) = rgdnv*ugdnv*ugdnv+pgdnv
    flux(3) = ugdnv*(0.5d0*rgdnv*(ugdnv*ugdnv+utrans1*utrans1+utrans2*utrans2) &
         + pgdnv/(gamma -1.d0) + pgdnv)
    flux(4) = rgdnv*ugdnv*utrans1
    flux(5) = rgdnv*ugdnv*utrans2
  end subroutine analriem

  subroutine hllriem(gamma,statel,stater,smallp,smallr,flux,  debug)
    logical, intent(in), optional :: debug
    real(rt), intent(in) :: statel(5),stater(5), smallp, smallr, gamma
    real(rt), intent(out) :: flux(5)
    
    real(rt) :: rl, pl, vnl, vtl, vttl, el, csl
    real(rt) :: rr, pr, vnr, vtr, vttr, er, csr
    real(rt) :: fl(5), fr(5), ul(5), ur(5)
    real(rt) :: alpha_plus, alpha_mins, alpha_pm, aainv

    rl = max(smallr,statel(1))
    vnl = statel(2)
    pl = max(smallp,statel(3))
    vtl = statel(4)
    vttl = statel(5)

    rr = max(smallr,stater(1))
    vnr = stater(2)
    pr = max(smallp,stater(3))
    vtr = stater(4)
    vttr = stater(5)

    el = pl/((gamma-1.d0)*rl)
    er = pr/((gamma-1.d0)*rr)
    csl = sqrt(gamma*pl/rl)
    csr = sqrt(gamma*pr/rr)

    ul(1) = rl
    ul(2) = rl*vnl
    ul(3) = rl*(el + 0.5d0*(vnl**2+vtl**2+vttl**2))
    ul(4) = rl*vtl
    ul(5) = rl*vttl

    ur(1) = rr
    ur(2) = rr*vnr
    ur(3) = rr*(er + 0.5d0*(vnr**2+vtr**2+vttr**2))
    ur(4) = rr*vtr
    ur(5) = rr*vttr

    alpha_plus = max(0.d0, csl+vnl, csr+vnr)
    alpha_mins = max(0.d0, csl-vnl, csr-vnr)

    ! LLF
    ! alpha_plus = max(alpha_plus,alpha_mins)
    ! alpha_mins = alpha_plus
    
    aainv = 1.d0/(alpha_plus+alpha_mins)
    alpha_pm = alpha_plus * alpha_mins * aainv
    alpha_plus = alpha_plus * aainv
    alpha_mins = 1.d0 - alpha_plus

    fl(1) = vnl*ul(1)
    fl(2) = vnl*ul(2)+pl
    fl(3) = vnl*(ul(3) + pl)
    fl(4) = vnl*ul(4)
    fl(5) = vnl*ul(5)

    fr(1) = vnr*ur(1)
    fr(2) = vnr*ur(2)+pr
    fr(3) = vnr*(ur(3) + pr)
    fr(4) = vnr*ur(4)
    fr(5) = vnr*ur(5)

    flux = alpha_plus * fl + alpha_mins * fr - alpha_pm*(ur-ul)

  end subroutine hllriem

end module riemann_module
