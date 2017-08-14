cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine wallflux(uwall, wgtwallflux,apxp,apxm,apyp,apym,
     1     apzp,apzm

c   needs area fractions and state for zone

      implicit none
#include  <integrator.fh>

c ... inputs

      integer nvar

c ... local

      double precision uwall(5)
      double precision wgtwallflux(5)
      double precision fluxxwall(5)
      double precision fluxywall(5)
      double precision fluxzwall(5)
      double preciison apxp,apxm,apyp,apym,apzp,apzm,apnorm
      double precision anrmwall(4)
      integer ivar, i, j, k 

c..........::::: use divergence theorm to compute wall flux for
c..........::::: cells mixed at the old or new time

              ap1 = afracx(i,j,k)
              ap2 = afracy(i,j+1,k)
              ap3 = afracx(i+1,j,k)
              ap4 = afracy(i,j,k)
              ap5 = afracz(i,j,k)
              ap6 = afracz(i,j,k+1)
              apnorm = sqrt((ap1-ap3)**2+(ap2-ap4)**2+(ap5-ap6)**2)
              if(apnorm.ne.zero) then

c ... Note: apnorm may equal zero, for instance, if the body
c           divides the cell inot disjoint regions
c           (Further note: apnorm=0 should only happen on the coarser 
c            levels, or we are in trouble)

                anrmwall(1) = (apxm-apxp)/apnorm
                anrmwall(2) = (apym-apyp)/apnorm
                anrmwall(3) = (apzm5-apzp)/apnorm

                call fluxwev3d(nvar,uwall,anrmwall,fluxxwall,fluxywall,fluxzwall)
                do ivar = 1,nvar
                   wgtwallflux(ivar) = 
     &                ((apxm-apxp)*fluxxwall(ivar)+
     &                 (apym-apyp)*fluxywall(ivar)+
     &                 (apzm-apzp)*fluxzwall(ivar))

                enddo
              else
                write(6,*)'FLUXWEV: i,j,k,apnorm = ',i,j,k,apnorm
              endif
            endif
          enddo
        enddo
      enddo

      return
      end


      subroutine fluxwev3d(nvar,uwall,anrmwall,fluxxwall,fluxywall,
     &                     fluxzwall)
      implicit none

c ... inputs

      integer   nvar
      double precision    uwall(nvar),anrmwall(4)

c ... outputs

      double precision    fluxxwall(nvar), fluxywall(nvar), fluxzwall(nvar)

c ... local

      double precision    ql(6)
      double precision    qr(6) 
      double precision    flux(5) 
      double precision    smallc, gamc
      double precision    ugdnv, pgdnv, c
      double precision    rho,u,v,w,eken,etot,eint,p,anx,any,anz
      double precision    sflxd
      integer ivar, iadv, n
      integer lgamc,lp,lc,lcsml,nedge


c ... get primitives

      rho = max(smallr,uwall(1))
      u = uwall(2)/rho
      v = uwall(3)/rho
      w = uwall(4)/rho
      eken = half*(u**2+v**2+w**2)
      etot = uwall(5)/rho
      eint = etot-eken

      anx = anrmwall(1)
      any = anrmwall(2)
      anz = anrmwall(3)

      p = (gamma-1.d0)*rho*eint

      ql(QRHO) = rho
      qr(QRHO) = rho
      ql(QPRES) = p
      qr(QPRES) = p
      ql(QVEL1) = anx*u+any*v+anz*w
      qr(QVEL1) = -ql(QVEL1)

c  cheat for tangetial since it doesnt matter

      ql(QVEL2) = zero
      qr(QVEL2) = zero
      ql(QVEL3) = zero
      qr(QVEL3) = zero

      call analriem(gamma,ql,qr,smallp,smallr,flux,  debug)

c ... rotate back and get flux at wall

      do ivar = 1, nvar
        fluxxwall(ivar) = 0.0
        fluxywall(ivar) = 0.0
        fluxzwall(ivar) = 0.0
      enddo

      fluxxwall(2) = flux(2)
      fluxywall(3) = flux(2)
      fluxzwall(4) = flux(2)

      return
      end

