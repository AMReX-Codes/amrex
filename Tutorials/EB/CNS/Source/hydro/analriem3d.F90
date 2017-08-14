
      subroutine analriem(gamma,pl,rl,ul,ut1l,ut2l,pr,rr,ur,ut1r,ut2r,
      subroutine analriem(gamma,statel,stater,smallp,smallr,flux,  debug)

      implicit none

      logical debug
      double precision statel(5),stater(5)
      double precision pl,rl,ul,pr,rr,ur,smallp,smallr, gamma
      double precision ut1l,ut2l,ut1r,ut2r

      double precision clsql,clsqr,vbarl,vbarr
      double precision pstar,pstnm1,ustarp,ustarm,ustar
      double precision wl,wr,wlsq,wrsq
      double precision weakwv,cleft,cright
      double precision a1,b1,c1,d1,e1,f1,ff
      double precision dpditer,zp,zm,denom,small
      double precision ustnm1,ustnp1
      double precision function ff
      double precision uo,ro,rstar
      double precision flux(5)
      double precision po,pgdnv,rgdnv,ugdnv,dpjmp,cstar
      double precision spin,spout,sgnm,co,tauo,ushock,wo, wosq
      double precision frac
      double precision utrans1,utrans2

      integer itno,iter

      data weakwv/1.e-03/
      data small/1.e-06/

      ff(a1,b1,c1,d1,e1,f1)=(a1*f1 + c1*d1 - c1*f1*(e1-b1))/(c1+f1)

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
      wl = dsqrt(clsql)
      wr = dsqrt(clsqr)
      cleft = wl/rl
      cright = wr/rr


      pstar=(wl*pr+wr*pl-wr*wl
     1      *(ur-ul))/(wl+wr)
      pstar=max(pstar,smallp)
       pstnm1 = pstar

      wlsq = (.5d0*(gamma-1.d0)*(pstar+pl)+pstar)/vbarl
      wrsq = (.5d0*(gamma-1.d0)*(pstar+pr)+pstar)/vbarr

      wl = dsqrt(wlsq)
      wr = dsqrt(wrsq)
      ustarp=ul-(pstar-pl)/wl
      ustarm=ur+(pstar-pr)/wr

      pstar = ff(pl,ul,wl,pr,ur,wr)
      pstar=max(pstar,smallp)


      itno = 3

      do 300 iter = 1,itno

      wlsq = (.5d0*(gamma-1.d0)*(pstar+pl)+pstar)/vbarl
      wrsq = (.5d0*(gamma-1.d0)*(pstar+pr)+pstar)/vbarr

      wl = 1./sqrt(wlsq)
      wr = 1./sqrt(wrsq)

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

300   continue

      if(debug)then
           write(6,*)"ustart,pstar", ustar,pstar
           write(6,*)" here"
           write(6,*)"lef", pl,rl,ul
           write(6,*)"rig", pr,rr,ur
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
      frac = 0.5d0*(1.d0+(spin+spout)/max(spout-spin,spin+spout,
     1           small*(cleft+cright)))

      if(debug)then
c          write(6,*)"cstar,sgnm,uo,co",cstar,sgnm,uo,co
           write(6,*)"dpjmp,spin,spout,frac", dpjmp,spin,spout,frac
      endif
      

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
      flux(3) = ugdnv*(0.5d0*rgdnv*ugdnv*ugdnv + pgdnv/(gamma -1.d0)
     1            + pgdnv)
      flux(4) = rgdnv*ugdnv*utran1
      flux(5) = rgdnv*ugdnv*utran2

      return
      end
