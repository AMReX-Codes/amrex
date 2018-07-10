
#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_ArrayLim.H"
#include "Derived.H"

#ifdef BL_DERIVE_IAMR
#define TRACER_POSITION 5
#else
#define TRACER_POSITION 6
#endif

      subroutine FORT_SETGAMMA (gam)
      REAL_T gam
#include "xxmeth.fh"	

      gamma = gam
      end

      subroutine FORT_SETALPHA (alp)
      REAL_T alp
#include "xxmeth.fh"	

      alpha = alp
      end

      subroutine FORT_SETQPOS (qpos)
      integer qpos
#include "xxmeth.fh"	

      qposition = qpos
      end
c
c     derive log density:  
c     
c     Inputs/Outputs
c     u            =>  (const)  state array
c     ulo,uhi      =>  (const)  index limits u
c     nvar         =>  (const)  number of variables in state array
c     xlo          =>  (const)  physical location of lo end of u
c     xhi          =>  (const)  physical location of hi end of u
c     dat         <=>  (modify) array holding derived data
c     dlo,dhi      =>  (const)  index limits of dat array
c     ovlo,ovhi    =>  (const)  subregion where derivation is done
c
      subroutine FORT_DERLGDN (u,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo, xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3, ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T u(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)

      integer i,j,k
      REAL_T smallr

      data smallr /1.0e-6/

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
                  dat(i,j,k) = log10(max(u(i,j,k,1),smallr))
               end do
            end do
         end do
      end
c
c     derive density:  
c     
c     Inputs/Outputs
c     u            =>  (const)  state array
c     ulo,uhi      =>  (const)  index limits of u
c     nvar         =>  (const)  number of variables in state array
c     xlo          =>  (const)  physical location of lo end of u
c     xhi          =>  (const)  physical location of hi end of u
c     dat         <=>  (modify) array holding derived data
c     dlo,dhi      =>  (const)  index limits of dat array
c     ovlo,ovhi    =>  (const)  subregion where derivation is done
c
      subroutine FORT_DERDEN (u,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo, xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3, ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T u(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)

      integer i,j,k

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
                  dat(i,j,k) = u(i,j,k,1)
               end do
            end do
         end do
      end
c
c     derive x velocity:  
c     
c     Inputs/Outputs
c     u            =>  (const)  state array
c     ulo,uhi      =>  (const)  index limits of u
c     nvar         =>  (const)  number of variables in state array
c     xlo          =>  (const)  physical location of lo end of u
c     xhi          =>  (const)  physical location of hi end of u
c     dat         <=>  (modify) array holding derived data
c     dlo,dhi      =>  (const)  index limits of dat array
c     ovlo,ovhi    =>  (const)  subregion where derivation is done
c
      subroutine FORT_DERXVEL (u,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo, xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3, ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T u(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3 :dhii3)

      integer i,j,k

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
                  dat(i,j,k) = u(i,j,k,2)/u(i,j,k,1)
               end do
            end do
         end do
      end
c
c     derive y velocity:  
c     
c     Inputs/Outputs
c     u            =>  (const)  state array
c     ulo,uhi      =>  (const)  index limits of u
c     nvar         =>  (const)  number of variables in state array
c     xlo          =>  (const)  physical location of lo end of u
c     xhi          =>  (const)  physical location of hi end of u
c     dat         <=>  (modify) array holding derived data
c     dlo,dhi      =>  (const)  index limits of dat array
c     ovlo,ovhi    =>  (const)  subregion where derivation is done
c
      subroutine FORT_DERYVEL (u,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo, xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3, ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T u(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)

      integer i,j,k

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
                  dat(i,j,k) = u(i,j,k,3)/u(i,j,k,1)
               end do
            end do
         end do
      end
c
c     derive z velocity:  
c     
c     Inputs/Outputs
c     u            =>  (const)  state array
c     ulo,uhi      =>  (const)  index limits of u
c     nvar         =>  (const)  number of variables in state array
c     xlo          =>  (const)  physical location of lo end of u
c     xhi          =>  (const)  physical location of hi end of u
c     dat         <=>  (modify) array holding derived data
c     dlo,dhi      =>  (const)  index limits of dat array
c     ovlo,ovhi    =>  (const)  subregion where derivation is done
c
      subroutine FORT_DERZVEL (u,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo, xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3, ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T u(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)

      integer i,j,k

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
                  dat(i,j,k) = u(i,j,k,4)/u(i,j,k,1)
               end do
            end do
         end do
      end
c
c     derive total energy  = (rho*E)/rho:  
c     
c     Inputs/Outputs
c     u            =>  (const)  state array
c     ulo,uhi      =>  (const)  index limits of u
c     nvar         =>  (const)  number of variables in state array
c     xlo          =>  (const)  physical location of lo end of u
c     xhi          =>  (const)  physical location of hi end of u
c     dat         <=>  (modify) array holding derived data
c     dlo,dhi      =>  (const)  index limits of dat array
c     ovlo,ovhi    =>  (const)  subregion where derivation is done
c
      subroutine FORT_DERTENG (u,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo,xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3,ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T u(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)

      integer i,j,k

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
                  dat(i,j,k) = u(i,j,k,5)/u(i,j,k,1)
               end do
            end do
         end do
      end
c
c     derive total energy  = (rho*(E-(u^2+v^2)/2)/rho:  
c     
c     Inputs/Outputs
c     u            =>  (const)  state array
c     ulo,uhi      =>  (const)  index limits of u
c     nvar         =>  (const)  number of variables in state array
c     xlo          =>  (const)  physical location of lo end of u
c     xhi          =>  (const)  physical location of hi end of u
c     dat         <=>  (modify) array holding derived data
c     dlo,dhi      =>  (const)  index limits of dat array
c     ovlo,ovhi    =>  (const)  subregion where derivation is done
c
      subroutine FORT_DERIENG (u,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo, xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3, ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T u(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)

      integer i,j,k
      REAL_T vx,vy,vz,eng

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
                  vx = u(i,j,k,2)/u(i,j,k,1)
                  vy = u(i,j,k,3)/u(i,j,k,1)
                  vz = u(i,j,k,4)/u(i,j,k,1)
                  eng = u(i,j,k,5)/u(i,j,k,1)
                  dat(i,j,k) = eng - (vx**2+vy**2+vz**2)*half
               end do
            end do
         end do
      end
c
c     derive kenetic energy  = rho(u**2 + v**2 + w**2)/2
c     
c     Inputs/Outputs
c     u            =>  (const)  state array
c     ulo,uhi      =>  (const)  index limits of u
c     nvar         =>  (const)  number of variables in state array
c     xlo          =>  (const)  physical location of lo end of u
c     xhi          =>  (const)  physical location of hi end of u
c     dat         <=>  (modify) array holding derived data
c     dlo,dhi      =>  (const)  index limits of dat array
c     ovlo,ovhi    =>  (const)  subregion where derivation is done
c
      subroutine FORT_DERKENG (u,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo,xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3,ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T u(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)

      integer i,j,k
      REAL_T eng

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
                  eng = u(i,j,k,2)**2 + u(i,j,k,3)**2 + u(i,j,k,4)**2
                  dat(i,j,k) = half*eng/u(i,j,k,1)
               end do
            end do
         end do
      end



c
c     derive vorticity
c     
c     Inputs/Outputs
c     u            =>  (const)  state array
c     ulo,uhi      =>  (const)  index limits of u
c     nvar         =>  (const)  number of variables in state array
c     xlo          =>  (const)  physical location of lo end of u
c     xhi          =>  (const)  physical location of hi end of u
c     dat         <=>  (modify) array holding derived data
c     dlo,dhi      =>  (const)  index limits of dat array
c     ovlo,ovhi    =>  (const)  subregion where derivation is done
c
c     Uses one-sided derivatives at edges.
c
      subroutine FORT_DERVORT (f,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo,xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3,ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T f(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)

      REAL_T uym,uyp,uy,uzm,uzp,uz
      REAL_T vxm,vxp,vx,vzm,vzp,vz
      REAL_T wxm,wxp,wx,wym,wyp,wy
      REAL_T uc,vc,wc
      REAL_T uy2m,uy2p,uz2m,uz2p
      REAL_T vx2m,vx2p,vz2m,vz2p
      REAL_T wx2m,wx2p,wy2m,wy2p
      REAL_T vorx,vory,vorz
      REAL_T ddx,ddy,ddz
      integer i,j,k

      ddx = half*float(uhii1-uloi1+1)/(xhi(1)-xlo(1))
      ddy = half*float(uhii2-uloi2+1)/(xhi(2)-xlo(2))
      ddz = half*float(uhii3-uloi3+1)/(xhi(3)-xlo(3))

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
		  uc = f(i,j,k,2)/f(i,j,k,1)

		  if (j.eq.ovlo(2).and.j+2.le.ovhi(2)) then
		     uy2p = f(i,j+2,k,2)/f(i,j+2,k,1)
		     uyp  = f(i,j+1,k,2)/f(i,j+1,k,1)
		     uy   = ddy*(-uy2p+4*uyp-3*uc)
		  else if (j.eq.ovlo(2)) then
		     uyp = f(i,j+1,k,2)/f(i,j+1,k,1)
		     uy  = 2*ddy*(uyp-uc)
		  else if (j.eq.ovhi(2).and.j-2.ge.ovlo(2)) then
		     uy2m = f(i,j-2,k,2)/f(i,j-2,k,1)
		     uym  = f(i,j-1,k,2)/f(i,j-1,k,1)
		     uy   = ddy*(+uy2m-4*uym+3*uc)
		  else if (j.eq.ovhi(2)) then
		     uym = f(i,j-1,k,2)/f(i,j-1,k,1)
		     uy  = 2*ddy*(-uym+uc)
		  else if (j.gt.ovlo(2).and.j.lt.ovhi(2) ) then
		     uyp = f(i,j+1,k,2)/f(i,j+1,k,1)
		     uym = f(i,j-1,k,2)/f(i,j-1,k,1)
		     uy  = ddy*(uyp-uym)
		  else
		     write(6,*) 'uy'
		     stop
		  end if

		  if (k.eq.ovlo(3).and.k+2.le.ovhi(3)) then
		     uz2p = f(i,j,k+2,2)/f(i,j,k+2,1)
		     uzp  = f(i,j,k+1,2)/f(i,j,k+1,1)
		     uz   = ddz*(-uz2p+4*uzp-3*uc)
		  else if (k.eq.ovlo(3)) then
		     uzp = f(i,j,k+1,2)/f(i,j,k+1,1)
		     uz  = 2*ddz*(uzp-uc)
		  else if (k.eq.ovhi(3).and.k-2.ge.ovlo(3)) then
		     uz2m = f(i,j,k-2,2)/f(i,j,k-2,1)
		     uzm  = f(i,j,k-1,2)/f(i,j,k-1,1)
		     uz   = ddz*(+uz2m-4*uzm+3*uc)
		  else if (k.eq.ovhi(3)) then
		     uzm = f(i,j,k-1,2)/f(i,j,k-1,1)
		     uz  = 2*ddz*(-uzm+uc)
		  else if (k.gt.ovlo(3) .and. k.lt.ovhi(3)) then
		     uzp = f(i,j,k+1,2)/f(i,j,k+1,1)
		     uzm = f(i,j,k-1,2)/f(i,j,k-1,1)
		     uz  = ddz*(uzp-uzm)
		  else
		     write(6,*) 'uz'
		     stop
		  end if

		  vc = f(i,j,k,3)/f(i,j,k,1)

		  if (i.eq.ovlo(1).and.i+2.le.ovhi(1)) then
		     vx2p = f(i+2,j,k,3)/f(i+2,j,k,1)
		     vxp = f(i+1,j,k,3)/f(i+1,j,k,1)
		     vx = ddx*(-vx2p+4*vxp-3*vc)
		  else if (i.eq.ovlo(1)) then
		     vxp = f(i+1,j,k,3)/f(i+1,j,k,1)
		     vx = 2*ddx*(vxp-vc)
		  else if (i.eq.ovhi(1).and.i-2.ge.ovlo(1)) then
		     vx2m = f(i-2,j,k,3)/f(i-2,j,k,1)
		     vxm = f(i-1,j,k,3)/f(i-1,j,k,1)
		     vx = ddx*(+vx2m-4*vxm+3*vc)
		  else if (i.eq.ovhi(1)) then
		     vxm = f(i-1,j,k,3)/f(i-1,j,k,1)
		     vx = 2*ddx*(-vxm+vc)
		  else if (i.gt.ovlo(1).and.i.lt.ovhi(1)) then
		     vxp = f(i+1,j,k,3)/f(i+1,j,k,1)
		     vxm = f(i-1,j,k,3)/f(i-1,j,k,1)
		     vx = ddx*(vxp-vxm)
		  else
		     write(6,*) 'vx'
		     stop
		  end if

		  if (k.eq.ovlo(3).and.k+2.le.ovhi(3)) then
		     vz2p = f(i,j,k+2,3)/f(i,j,k+2,1)
		     vzp  = f(i,j,k+1,3)/f(i,j,k+1,1)
		     vz   = ddz*(-vz2p+4*vzp-3*vc)
		  else if (k.eq.ovlo(3)) then
		     vzp = f(i,j,k+1,3)/f(i,j,k+1,1)
		     vz  = 2*ddz*(vzp-vc)
		  else if (k.eq.ovhi(3).and.k-2.ge.ovlo(3)) then
		     vz2m = f(i,j,k-2,3)/f(i,j,k-2,1)
		     vzm  = f(i,j,k-1,3)/f(i,j,k-1,1)
		     vz   = ddz*(+vz2m-4*vzm+3*vc)
		  else if (k.eq.ovhi(3)) then
		     vzm = f(i,j,k-1,3)/f(i,j,k-1,1)
		     vz  = 2*ddz*(-vzm+vc)
		  else if (k.gt.ovlo(3).and.k.lt.ovhi(3)) then
		     vzp = f(i,j,k+1,3)/f(i,j,k+1,1)
		     vzm = f(i,j,k-1,3)/f(i,j,k-1,1)
		     vz  = ddz*(vzp-vzm)
		  else
		     write(6,*) 'vz'
		     stop
		  end if

                  wc = f(i,j,k,4)/f(i,j,k,1)

		  if (i.eq.ovlo(1).and.i+2.le.ovhi(1)) then
		     wx2p = f(i+2,j,k,4)/f(i+2,j,k,1)
		     wxp  = f(i+1,j,k,4)/f(i+1,j,k,1)
		     wx   = ddx*(-wx2p+4*wxp-3*wc)
		  else if (i.eq.ovlo(1)) then
		     wxp = f(i+1,j,k,4)/f(i+1,j,k,1)
		     wx  = 2*ddx*(wxp-wc)
		  else if (i.eq.ovhi(1).and.i-2.ge.ovlo(1)) then
		     wx2m = f(i-2,j,k,4)/f(i-2,j,k,1)
		     wxm  = f(i-1,j,k,4)/f(i-1,j,k,1)
		     wx   = ddx*(+wx2m-4*wxm+3*wc)
		  else if (i.eq.ovhi(1)) then
		     wxm = f(i-1,j,k,4)/f(i-1,j,k,1)
		     wx  = 2*ddx*(-wxm+wc)
		  else if (i.gt.ovlo(1).and.i.lt.ovhi(1)) then
		     wxp = f(i+1,j,k,4)/f(i+1,j,k,1)
		     wxm = f(i-1,j,k,4)/f(i-1,j,k,1)
		     wx  = ddx*(wxp-wxm)
		  else
		     write(6,*) 'wx'
		     stop
		  end if

		  if (j.eq.ovlo(2).and.j+2.le.ovhi(2)) then
		     wy2p = f(i,j+2,k,4)/f(i,j+2,k,1)
		     wyp  = f(i,j+1,k,4)/f(i,j+1,k,1)
		     wy   = ddy*(-wy2p+4*wyp-3*wc)
		  else if (j.eq.ovlo(2)) then
		     wyp = f(i,j+1,k,4)/f(i,j+1,k,1)
		     wy  = 2*ddy*(wyp-wc)
		  else if (j.eq.ovhi(2).and.j-2.ge.ovlo(2)) then
		     wy2m = f(i,j-2,k,4)/f(i,j-2,k,1)
		     wym  = f(i,j-1,k,4)/f(i,j-1,k,1)
		     wy   = ddy*(+wy2m-4*wym+3*wc)
		  else if (j.eq.ovhi(2)) then
		     wym = f(i,j-1,k,4)/f(i,j-1,k,1)
		     wy  = 2*ddy*(-wym+wc)
		  else if (j.gt.ovlo(2).and.j.lt.ovhi(2)) then
		     wyp = f(i,j+1,k,4)/f(i,j+1,k,1)
		     wym = f(i,j-1,k,4)/f(i,j-1,k,1)
		     wy  = ddy*(wyp-wym)
		  else
		     write(6,*) 'wy'
		     stop
		  end if

                  vorx = wy-vz
                  vory = uz-wx
                  vorz = vx-uy

                  dat(i,j,k) = sqrt(vorx**2+vory**2+vorz**2)
               end do
            end do
         end do
      end

c
c     derive vorticity components
c     
c     Inputs/Outputs
c     u            =>  (const)  state array
c     ulo,uhi      =>  (const)  index limits of u
c     nvar         =>  (const)  number of variables in state array
c     xlo          =>  (const)  physical location of lo end of u
c     xhi          =>  (const)  physical location of hi end of u
c     dat         <=>  (modify) array holding derived data
c     dlo,dhi      =>  (const)  index limits of dat array
c     ovlo,ovhi    =>  (const)  subregion where derivation is done
c
c     Uses one-sided derivatives at edges.
c
      subroutine FORT_DERVORTCOMPS (f,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo,xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3,ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T f(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3, 3)

      REAL_T uym,uyp,uy,uzm,uzp,uz
      REAL_T vxm,vxp,vx,vzm,vzp,vz
      REAL_T wxm,wxp,wx,wym,wyp,wy
      REAL_T uc,vc,wc
      REAL_T uy2m,uy2p,uz2m,uz2p
      REAL_T vx2m,vx2p,vz2m,vz2p
      REAL_T wx2m,wx2p,wy2m,wy2p
      REAL_T vorx,vory,vorz
      REAL_T ddx,ddy,ddz
      integer i,j,k

      ddx = half*float(uhii1-uloi1+1)/(xhi(1)-xlo(1))
      ddy = half*float(uhii2-uloi2+1)/(xhi(2)-xlo(2))
      ddz = half*float(uhii3-uloi3+1)/(xhi(3)-xlo(3))

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
		  uc = f(i,j,k,2)/f(i,j,k,1)

		  if (j.eq.ovlo(2).and.j+2.le.ovhi(2)) then
		     uy2p = f(i,j+2,k,2)/f(i,j+2,k,1)
		     uyp  = f(i,j+1,k,2)/f(i,j+1,k,1)
		     uy   = ddy*(-uy2p+4*uyp-3*uc)
		  else if (j.eq.ovlo(2)) then
		     uyp = f(i,j+1,k,2)/f(i,j+1,k,1)
		     uy  = 2*ddy*(uyp-uc)
		  else if (j.eq.ovhi(2).and.j-2.ge.ovlo(2)) then
		     uy2m = f(i,j-2,k,2)/f(i,j-2,k,1)
		     uym  = f(i,j-1,k,2)/f(i,j-1,k,1)
		     uy   = ddy*(+uy2m-4*uym+3*uc)
		  else if (j.eq.ovhi(2)) then
		     uym = f(i,j-1,k,2)/f(i,j-1,k,1)
		     uy  = 2*ddy*(-uym+uc)
		  else if (j.gt.ovlo(2).and.j.lt.ovhi(2) ) then
		     uyp = f(i,j+1,k,2)/f(i,j+1,k,1)
		     uym = f(i,j-1,k,2)/f(i,j-1,k,1)
		     uy  = ddy*(uyp-uym)
		  else
		     write(6,*) 'uy'
		     stop
		  end if

		  if (k.eq.ovlo(3).and.k+2.le.ovhi(3)) then
		     uz2p = f(i,j,k+2,2)/f(i,j,k+2,1)
		     uzp  = f(i,j,k+1,2)/f(i,j,k+1,1)
		     uz   = ddz*(-uz2p+4*uzp-3*uc)
		  else if (k.eq.ovlo(3)) then
		     uzp = f(i,j,k+1,2)/f(i,j,k+1,1)
		     uz  = 2*ddz*(uzp-uc)
		  else if (k.eq.ovhi(3).and.k-2.ge.ovlo(3)) then
		     uz2m = f(i,j,k-2,2)/f(i,j,k-2,1)
		     uzm  = f(i,j,k-1,2)/f(i,j,k-1,1)
		     uz   = ddz*(+uz2m-4*uzm+3*uc)
		  else if (k.eq.ovhi(3)) then
		     uzm = f(i,j,k-1,2)/f(i,j,k-1,1)
		     uz  = 2*ddz*(-uzm+uc)
		  else if (k.gt.ovlo(3) .and. k.lt.ovhi(3)) then
		     uzp = f(i,j,k+1,2)/f(i,j,k+1,1)
		     uzm = f(i,j,k-1,2)/f(i,j,k-1,1)
		     uz  = ddz*(uzp-uzm)
		  else
		     write(6,*) 'uz'
		     stop
		  end if

		  vc = f(i,j,k,3)/f(i,j,k,1)

		  if (i.eq.ovlo(1).and.i+2.le.ovhi(1)) then
		     vx2p = f(i+2,j,k,3)/f(i+2,j,k,1)
		     vxp = f(i+1,j,k,3)/f(i+1,j,k,1)
		     vx = ddx*(-vx2p+4*vxp-3*vc)
		  else if (i.eq.ovlo(1)) then
		     vxp = f(i+1,j,k,3)/f(i+1,j,k,1)
		     vx = 2*ddx*(vxp-vc)
		  else if (i.eq.ovhi(1).and.i-2.ge.ovlo(1)) then
		     vx2m = f(i-2,j,k,3)/f(i-2,j,k,1)
		     vxm = f(i-1,j,k,3)/f(i-1,j,k,1)
		     vx = ddx*(+vx2m-4*vxm+3*vc)
		  else if (i.eq.ovhi(1)) then
		     vxm = f(i-1,j,k,3)/f(i-1,j,k,1)
		     vx = 2*ddx*(-vxm+vc)
		  else if (i.gt.ovlo(1).and.i.lt.ovhi(1)) then
		     vxp = f(i+1,j,k,3)/f(i+1,j,k,1)
		     vxm = f(i-1,j,k,3)/f(i-1,j,k,1)
		     vx = ddx*(vxp-vxm)
		  else
		     write(6,*) 'vx'
		     stop
		  end if

		  if (k.eq.ovlo(3).and.k+2.le.ovhi(3)) then
		     vz2p = f(i,j,k+2,3)/f(i,j,k+2,1)
		     vzp  = f(i,j,k+1,3)/f(i,j,k+1,1)
		     vz   = ddz*(-vz2p+4*vzp-3*vc)
		  else if (k.eq.ovlo(3)) then
		     vzp = f(i,j,k+1,3)/f(i,j,k+1,1)
		     vz  = 2*ddz*(vzp-vc)
		  else if (k.eq.ovhi(3).and.k-2.ge.ovlo(3)) then
		     vz2m = f(i,j,k-2,3)/f(i,j,k-2,1)
		     vzm  = f(i,j,k-1,3)/f(i,j,k-1,1)
		     vz   = ddz*(+vz2m-4*vzm+3*vc)
		  else if (k.eq.ovhi(3)) then
		     vzm = f(i,j,k-1,3)/f(i,j,k-1,1)
		     vz  = 2*ddz*(-vzm+vc)
		  else if (k.gt.ovlo(3).and.k.lt.ovhi(3)) then
		     vzp = f(i,j,k+1,3)/f(i,j,k+1,1)
		     vzm = f(i,j,k-1,3)/f(i,j,k-1,1)
		     vz  = ddz*(vzp-vzm)
		  else
		     write(6,*) 'vz'
		     stop
		  end if

                  wc = f(i,j,k,4)/f(i,j,k,1)

		  if (i.eq.ovlo(1).and.i+2.le.ovhi(1)) then
		     wx2p = f(i+2,j,k,4)/f(i+2,j,k,1)
		     wxp  = f(i+1,j,k,4)/f(i+1,j,k,1)
		     wx   = ddx*(-wx2p+4*wxp-3*wc)
		  else if (i.eq.ovlo(1)) then
		     wxp = f(i+1,j,k,4)/f(i+1,j,k,1)
		     wx  = 2*ddx*(wxp-wc)
		  else if (i.eq.ovhi(1).and.i-2.ge.ovlo(1)) then
		     wx2m = f(i-2,j,k,4)/f(i-2,j,k,1)
		     wxm  = f(i-1,j,k,4)/f(i-1,j,k,1)
		     wx   = ddx*(+wx2m-4*wxm+3*wc)
		  else if (i.eq.ovhi(1)) then
		     wxm = f(i-1,j,k,4)/f(i-1,j,k,1)
		     wx  = 2*ddx*(-wxm+wc)
		  else if (i.gt.ovlo(1).and.i.lt.ovhi(1)) then
		     wxp = f(i+1,j,k,4)/f(i+1,j,k,1)
		     wxm = f(i-1,j,k,4)/f(i-1,j,k,1)
		     wx  = ddx*(wxp-wxm)
		  else
		     write(6,*) 'wx'
		     stop
		  end if

		  if (j.eq.ovlo(2).and.j+2.le.ovhi(2)) then
		     wy2p = f(i,j+2,k,4)/f(i,j+2,k,1)
		     wyp  = f(i,j+1,k,4)/f(i,j+1,k,1)
		     wy   = ddy*(-wy2p+4*wyp-3*wc)
		  else if (j.eq.ovlo(2)) then
		     wyp = f(i,j+1,k,4)/f(i,j+1,k,1)
		     wy  = 2*ddy*(wyp-wc)
		  else if (j.eq.ovhi(2).and.j-2.ge.ovlo(2)) then
		     wy2m = f(i,j-2,k,4)/f(i,j-2,k,1)
		     wym  = f(i,j-1,k,4)/f(i,j-1,k,1)
		     wy   = ddy*(+wy2m-4*wym+3*wc)
		  else if (j.eq.ovhi(2)) then
		     wym = f(i,j-1,k,4)/f(i,j-1,k,1)
		     wy  = 2*ddy*(-wym+wc)
		  else if (j.gt.ovlo(2).and.j.lt.ovhi(2)) then
		     wyp = f(i,j+1,k,4)/f(i,j+1,k,1)
		     wym = f(i,j-1,k,4)/f(i,j-1,k,1)
		     wy  = ddy*(wyp-wym)
		  else
		     write(6,*) 'wy'
		     stop
		  end if

                  vorx = wy-vz
                  vory = uz-wx
                  vorz = vx-uy

CCCC     dat(i,j,k) = sqrt(vorx**2+vory**2+vorz**2)

CCCC              dat(i,j,k,1) = vorx
CCCC              dat(i,j,k,2) = vory
CCCC              dat(i,j,k,3) = vorz

                  dat(i,j,k,1) = sqrt(vorx**2)
                  dat(i,j,k,2) = sqrt(vory**2)
                  dat(i,j,k,3) = sqrt(vorz**2)
CCC      dat(i,j,k,1) = sqrt(vorx**2+vory**2+vorz**2)
CCC      dat(i,j,k,2) = sqrt(vorx**2+vory**2+vorz**2)
CCC      dat(i,j,k,3) = sqrt(vorx**2+vory**2+vorz**2)

               end do
            end do
         end do
      end



c
c     derive concentration variable
c
c     Inputs/Outputs
c     u            =>  (const)  state array
c     ulo,uhi      =>  (const)  index limits of u
c     nvar         =>  (const)  number of variables in state array
c     xlo          =>  (const)  physical location of lo end of u
c     xhi          =>  (const)  physical location of hi end of u
c     dat         <=>  (modify) array holding derived data
c     dlo,dhi      =>  (const)  index limits of dat array
c     ovlo,ovhi    =>  (const)  subregion where derivation is done
c
      subroutine FORT_DERCON (u,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo,xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3,ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T u(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)

      integer i,j,k

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
c                 dat(i,j,k) = u(i,j,k,nvar)/u(i,j,k,1)
c JBB
                  dat(i,j,k) = u(i,j,k,nvar)
               end do
            end do
         end do
      end
c
c     INTEGRATE (3D version)
c   
c     integral = sum { field(i,j,k) * dx * dy * dz }
c   
c     Inputs / Outputs:
c   
c     field       => field of values to integrate
c     DIMS(field) => index limits of field
c     DIMS(valid) => index limits of valid data in field
c     delta       => cell size
c     integral    <= integral
c
      subroutine FORT_INTEGRATE
     +   (field, DIMS(field), DIMS(valid), delta, integral)
      integer DIMDEC(field)
      integer DIMDEC(valid)
      REAL_T  delta(3), integral
      REAL_T  field(DIMV(field))

      integer i, j, k
      REAL_T  sum
      
      sum = zero
      
      do k = ARG_L3(valid), ARG_H3(valid)
         do j = ARG_L2(valid), ARG_H2(valid)
            do i = ARG_L1(valid), ARG_H1(valid)
               sum = sum + field(i,j,k)
            end do
         end do
      end do

      integral = delta(1) * delta(2) * delta(3) * sum

      end
c
c     JPDF (3D version)
c   
c     Inputs / Outputs:
c   
c     field       => field of values to bin (1 = horizontal, 2 = vertical)
c     visible     => flags indicating visibility of cells (0.0 no, 1.0 yes)
c     DIMS(field) => index limits of field
c     DIMS(valid) => index limits of visible, and of valid data in field
c     hBins       => number of horizontal bins
c     hMax        => upper limit for horizontal bins
c     hMin        => lower limit for horizontal bins
c     vBins       => number of vertical bins
c     vMax        => upper limit for vertical bins
c     vMin        => lower limit for vertical bins
c     excessCount <= count of cells outside the bin limits
c     jointCount  <= count of cells inside the individual, two-way bins
c
      subroutine FORT_JPDF
     +     (field, visible, DIMS(field), DIMS(valid), 
     +     hBins, hMax, hMin, vBins, vMax, vMin,
     +     excessCount, jointCount)
      double precision  foobar
      integer DIMDEC(field)
      integer DIMDEC(valid)
      integer hBins, vBins
      integer excessCount
      integer jointCount
      REAL_T  hMax, hMin, vMax, vMin
      REAL_T  field(DIMV(field), 2)
c     the visible box is assumed sized to the valid region
      REAL_T  visible(DIMV(valid))
      dimension jointCount (0 : hBins * vBins - 1)

      integer bin, hBin, i, j, k, vBin
      REAL_T h1, h2, hVal, v1, v2, vVal

      h1 = (hBins - half) / (hMax - hMin)
      h2 = (- half) / (hMax - hMin)
      v1 = (vBins - half) / (vMax - vMin)
      v2 = (- half) / (vMax - vMin)

      do k = ARG_L3(valid), ARG_H3(valid)
         do j = ARG_L2(valid), ARG_H2(valid)
            do i = ARG_L1(valid), ARG_H1(valid)
               if (visible(i, j, k) .eq. 1) then
                  hVal = field(i, j, k, 1)
                  if (hVal .eq. hMin) then
                     hBin = 0
                  else if (hVal .eq. hMax) then
                     hBin = hBins - 1
                  else
c     if either of the Min and Max limits is an extrema, then
c     rounding may choose an out of range bin when none is 
c     intended.  hence the need to deal with this case, above. 
                     hBin = nint ((hVal - hMin) * h1 + (hMax - hVal) * h2)
                  end if
                  
                  vVal = field(i, j, k, 2)
                  if (vVal .eq. vMin) then
                     vBin = 0
                  else if (vVal .eq. vMax) then
                     vBin = vBins - 1
                  else
c     if either of the Min and Max limits is an extrema, then
c     rounding may choose an out of range bin when none is 
c     intended.  hence the need to deal with this case, above. 
                     vBin = nint ((vVal - vMin) * v1 + (vMax - vVal) * v2)
                  end if
                  
                  if (0 .le. hBin .and. hBin .lt. hBins .and.
     +                 0 .le. vBin .and. vBin .lt. vBins) then
                     bin = hBin + hBins * vBin
                     jointCount(bin) = jointCount(bin) + 1
                  else
                     excessCount = excessCount + 1
                  end if
               end if
            end do
         end do
      end do
      
      return
      end
c
c     LIMITS (3D version)
c   
c     maximum = max { field(i,j,k) }
c     minimum = min { field(i,j,k) }
c   
c     Inputs / Outputs:
c   
c     field       => field of values to integrate
c     visible     => flags indicating visibility of cells (0.0 no, 1.0 yes)
c     DIMS(field) => index limits of field
c     DIMS(valid) => index limits of visible, and of valid data in field
c     maximum     <= max of starting value and field cells that are visible
c     minimum     <= min of starting value and field cells that are visible
c
      subroutine FORT_LIMITS
     +   (field, visible, DIMS(field), DIMS(valid), maximum, minimum)
      integer DIMDEC(field)
      integer DIMDEC(valid)
      integer ifound
      REAL_T  maximum, minimum
      REAL_T  field(DIMV(field))
c     the visible box is assumed sized to the valid region
      REAL_T  visible(DIMV(valid))

      integer i, j, k

      do k = ARG_L3(valid), ARG_H3(valid)
         do j = ARG_L2(valid), ARG_H2(valid)
            do i = ARG_L1(valid), ARG_H1(valid)
               if (visible(i, j, k) .eq. 1) then
                  maximum = max (maximum, field(i, j, k))
                  minimum = min (minimum, field(i, j, k))
               end if
            end do
         end do
      end do

      end
c
c     SUMMASS
c   
c     MASS = sum{ vol(i,j)*rho(i,j) }
c   
c     Inputs / Outputs:
c   
c     rho        => density field
c     rlo,rhi    => index limits of rho aray
c     lo,hi      => index limits of grid interior
c     delta      => cell size
c     sum        <=  total mass
c     tmp        => temp column array
c
      subroutine FORT_SUMMASS (rho,DIMS(rho),DIMS(grid),delta,sum,tmp)
      integer DIMDEC(rho)
      integer DIMDEC(grid)
      REAL_T  sum, delta(3)
      REAL_T  rho(DIMV(rho))
      REAL_T  tmp(DIM2(grid))
      
      integer i, j, k
      REAL_T  dx, dy, dz, vol
      
      dx = delta(1)
      dy = delta(2)
      dz = delta(3)
      vol = dx*dy*dz
      
      do j = ARG_L2(grid), ARG_H2(grid)
         tmp(j) = zero
      end do
      
      do k = ARG_L3(grid), ARG_H3(grid)
         do i = ARG_L1(grid), ARG_H1(grid)
            do j = ARG_L2(grid), ARG_H2(grid)
               tmp(j) = tmp(j) + vol*rho(i,j,k)
            end do
         end do
      end do
      
      sum = zero
      do j = ARG_L2(grid), ARG_H2(grid)
         sum = sum + tmp(j)
      end do
      end
c     
c     SUMCONC
c     
c     CONC = sum{ vol(i,j,k)*conc(i,j,k)*(1-conc(i,j,k) }
c     
c     Inputs / Outputs:
c     
c     conc       => concentration
c     rlo,rhi    => index limits of rho aray
c     lo,hi      => index limits of grid interior
c     delta      => cell size
c     sum        <= concentration sum
c     tmp        => temp column array
c     
      subroutine FORT_SUMCONC (conc,DIMS(conc),DIMS(grid),delta,sum,tmp)
      integer DIMDEC(conc)
      integer DIMDEC(grid)
      REAL_T  sum, delta(3)
      REAL_T  conc(DIMV(conc))
      REAL_T  tmp(DIM2(grid))
      
      integer i, j, k
      REAL_T  dx, dy, dz, vol
      
      dx = delta(1)
      dy = delta(2)
      dz = delta(3)
      vol = dx*dy*dz
      
      do j = ARG_L2(grid), ARG_H2(grid)
         tmp(j) = zero
      end do
      
      do k = ARG_L3(grid), ARG_H3(grid)
         do i = ARG_L1(grid), ARG_H1(grid)
            do j = ARG_L2(grid), ARG_H2(grid)
c     tmp(j) = tmp(j) + vol * conc(i,j,k) * (1-conc(i,j,k))
c     JBB
               tmp(j) = tmp(j) + vol * conc(i,j,k) 
            end do
         end do
      end do
      
      sum = zero
      do j = ARG_L2(grid), ARG_H2(grid)
         sum = sum + tmp(j)
      end do
      end
c     
c     SUMVORT
c     
c     VORT = sum{ vol(i,j,k)*vort(i,j,k)**2 }
c     
c     Inputs / Outputs:
c     
c     vort       => vorticity
c     vlo,vhi    => index limits of vorticity array
c     lo,hi      => index limits of grid interior
c     delta      => cell size
c     sum        <= vorticity sum
c     tmp        => temp column array
c     
      subroutine FORT_SUMVORT (vort,DIMS(vort),DIMS(grid),delta,sum,tmp)
      integer DIMDEC(vort)
      integer DIMDEC(grid)
      REAL_T  sum, delta(3)
      REAL_T  vort(DIMV(vort))
      REAL_T  tmp(DIM2(grid))
      
      integer i, j, k
      REAL_T  dx, dy, dz, vol
      
      dx = delta(1)
      dy = delta(2)
      dz = delta(3)
      vol = dx*dy*dz
      
      do j = ARG_L2(grid), ARG_H2(grid)
         tmp(j) = zero
      end do
      
      do k = ARG_L3(grid), ARG_H3(grid)
         do i = ARG_L1(grid), ARG_H1(grid)
            do j = ARG_L2(grid), ARG_H2(grid)
               tmp(j) = tmp(j) + vol*vort(i,j,k) * vort(i,j,k)
            end do
         end do
      end do
      
      sum = zero
      do j = ARG_L2(grid), ARG_H2(grid)
         sum = sum + tmp(j)
      end do
      end
      
c     
c     derive dilitation
c     
c     Inputs/Outputs
c     u            =>  (const)  state array
c     ulo,uhi      =>  (const)  index limits of u
c     nvar         =>  (const)  number of variables in state array
c     xlo          =>  (const)  physical location of lo end of u
c     xhi          =>  (const)  physical location of hi end of u
c     dat         <=>  (modify) array holding derived data
c     dlo,dhi      =>  (const)  index limits of dat array
c     ovlo,ovhi    =>  (const)  subregion where derivation is done
c
c     Uses one-sided derivatives at edges.
c
      subroutine FORT_DERDILI (f,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo,xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3,ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T f(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)

      REAL_T ux2m,ux2p,uxm,uxp,ux
      REAL_T vy2m,vy2p,vym,vyp,vy
      REAL_T wz2m,wz2p,wzm,wzp,wz
      REAL_T uc,vc,wc
      REAL_T ddx,ddy,ddz
      integer i,j,k

      ddx = half*float(uhii1-uloi1+1)/(xhi(1)-xlo(1))
      ddy = half*float(uhii2-uloi2+1)/(xhi(2)-xlo(2))
      ddz = half*float(uhii3-uloi3+1)/(xhi(3)-xlo(3))

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
		  uc = f(i,j,k,2)/f(i,j,k,1)

		  if (i.eq.ovlo(1).and.i+2.le.ovhi(1)) then
		     ux2p = f(i+2,j,k,2)/f(i+2,j,k,1)
		     uxp = f(i+1,j,k,2)/f(i+1,j,k,1)
		     ux = ddx*(-ux2p+4*uxp-3*uc)
		  else if (i.eq.ovlo(1)) then
		     uxp = f(i+1,j,k,2)/f(i+1,j,k,1)
		     ux = 2*ddx*(uxp-uc)
		  else if (i.eq.ovhi(1).and.i-2.ge.ovlo(1)) then
		     ux2m = f(i-2,j,k,2)/f(i-2,j,k,1)
		     uxm = f(i-1,j,k,2)/f(i-1,j,k,1)
		     ux = ddx*(+ux2m-4*uxm+3*uc)
		  else if (i.eq.ovhi(1)) then
		     uxm = f(i-1,j,k,2)/f(i-1,j,k,1)
		     ux = 2*ddx*(-uxm+uc)
		  else if (i.gt.ovlo(1).and.i.lt.ovhi(1)) then
		     uxp = f(i+1,j,k,2)/f(i+1,j,k,1)
		     uxm = f(i-1,j,k,2)/f(i-1,j,k,1)
		     ux = ddx*(uxp-uxm)
		  else
		     write(6,*) 'ux'
		     stop
		  end if

		  vc = f(i,j,k,3)/f(i,j,k,1)

		  if (j.eq.ovlo(2).and.j+2.le.ovhi(2)) then
		     vy2p = f(i,j+2,k,3)/f(i,j+2,k,1)
		     vyp  = f(i,j+1,k,3)/f(i,j+1,k,1)
		     vy   = ddy*(-vy2p+4*vyp-3*vc)
		  else if (j.eq.ovlo(2)) then
		     vyp = f(i,j+1,k,3)/f(i,j+1,k,1)
		     vy  = 2*ddy*(vyp-vc)
		  else if (j.eq.ovhi(2).and.j-2.ge.ovlo(2)) then
		     vy2m = f(i,j-2,k,3)/f(i,j-2,k,1)
		     vym  = f(i,j-1,k,3)/f(i,j-1,k,1)
		     vy   = ddy*(+vy2m-4*vym+3*vc)
		  else if (j.eq.ovhi(2)) then
		     vym = f(i,j-1,k,3)/f(i,j-1,k,1)
		     vy  = 2*ddy*(-vym+vc)
		  else if (j.gt.ovlo(2).and.j.lt.ovhi(2) ) then
		     vyp = f(i,j+1,k,3)/f(i,j+1,k,1)
		     vym = f(i,j-1,k,3)/f(i,j-1,k,1)
		     vy  = ddy*(vyp-vym)
		  else
		     write(6,*) 'vy'
		     stop
		  end if

		  wc = f(i,j,k,4)/f(i,j,k,1)

		  if (k.eq.ovlo(3).and.k+2.le.ovhi(3)) then
		     wz2p = f(i,j,k+2,4)/f(i,j,k+2,1)
		     wzp  = f(i,j,k+1,4)/f(i,j,k+1,1)
		     wz   = ddz*(-wz2p+4*wzp-3*wc)
		  else if (k.eq.ovlo(3)) then
		     wzp = f(i,j,k+1,4)/f(i,j,k+1,1)
		     wz  = 2*ddz*(wzp-wc)
		  else if (k.eq.ovhi(3).and.k-2.ge.ovlo(3)) then
		     wz2m = f(i,j,k-2,4)/f(i,j,k-2,1)
		     wzm  = f(i,j,k-1,4)/f(i,j,k-1,1)
		     wz   = ddz*(+wz2m-4*wzm+3*wc)
		  else if (k.eq.ovhi(3)) then
		     wzm = f(i,j,k-1,4)/f(i,j,k-1,1)
		     wz  = 2*ddz*(-wzm+wc)
		  else if (k.gt.ovlo(3) .and. k.lt.ovhi(3)) then
		     wzp = f(i,j,k+1,4)/f(i,j,k+1,1)
		     wzm = f(i,j,k-1,4)/f(i,j,k-1,1)
		     wz  = ddz*(wzp-wzm)
		  else
		     write(6,*) 'wz'
		     stop
		  end if

                  dat(i,j,k) = ux + vy + wz
               end do
            end do
         end do
      end

c
c      The following routine does NOT have the signature of a derived function.
c
c      BINTRAC
c   
c      Bin the tracer for concentration in the following intervals:
c
c           (0.25, 0.75), (0.1, 0.9) & (0.05, 0.95)
c   
c      Inputs / Outputs:
c   
c      trac       => tracer (rho(i,j,k)*conc(i,j,k))
c      rho        => density
c      DIMS(c)    => index limits of rho and trac arrays
c      lo,hi      => index limits of grid interior
c      delta      => cell size
c      sum        <= concentration sum
c
       subroutine FORT_BINTRAC (trac,DIMS(c),rho,lo,hi,delta,sum)

       integer DIMDEC(c)
       integer lo(3), hi(3)
       REAL_T  sum(3), delta(3)
       REAL_T  trac(DIMV(c))
       REAL_T  rho(DIMV(c))

       integer i, j, k
       REAL_T  vol,conc

       vol = delta(1)*delta(2)*delta(3)

       do i = 1, 3
          sum(i) = zero
       end do

       do k = lo(3), hi(3)
          do i = lo(1), hi(1)
             do j = lo(2), hi(2)

                if (rho(i,j,k) .gt. 0) then

                   conc = trac(i,j,k) / rho(i,j,k)
                
                   if ((conc .gt. .25) .and. (conc .lt. .75)) then
                      sum(1) = sum(1) + vol * trac(i,j,k)
                   end if

                   if ((conc .gt. .1) .and. (conc .lt. .9)) then
                      sum(2) = sum(2) + vol * trac(i,j,k)
                   end if

                   if ((conc .gt. .05) .and. (conc .lt. .95)) then
                      sum(3) = sum(3) + vol * trac(i,j,k)
                   end if

                end if

             end do
          end do
       end do

       end
c
c     derive tracer:  
c     
c     Inputs/Outputs
c     u            =>  (const)  state array
c     ulo,uhi      =>  (const)  index limits of u
c     nvar         =>  (const)  number of variables in state array
c     xlo          =>  (const)  physical location of lo end of u
c     xhi          =>  (const)  physical location of hi end of u
c     dat         <=>  (modify) array holding derived data
c     dlo,dhi      =>  (const)  index limits of dat array
c     ovlo,ovhi    =>  (const)  subregion where derivation is done
c
      subroutine FORT_DERTRAC (u,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo, xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3, ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T u(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)

      integer i,j,k

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
                  dat(i,j,k) = u(i,j,k,TRACER_POSITION)
               end do
            end do
         end do
      end



c
      subroutine FORT_DERTMASK (u,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo, xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3, ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T u(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)
      REAL_T ryco, ryco2

      integer i,j,k

#include "xxmeth.fh"	

c     write(6,*) 'FORT_DERTMASK:  using alpha = ', alpha

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
                  ryco2 = u(i,j,k,9)
                  ryco  = u(i,j,k,11)
		  if (abs(ryco2).gt.alpha*abs(ryco)) then
                    dat(i,j,k) = u(i,j,k,6)
		  else
                    dat(i,j,k) = 0.0
		  endif
               end do
            end do
         end do
      end


c
      subroutine FORT_DERSTATE (u,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo, xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3, ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T u(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)

      integer i,j,k

#include "xxmeth.fh"	

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
                  dat(i,j,k) = u(i,j,k,qposition)
               end do
            end do
         end do
      end




c
c                          sqrt(u**2 + v**2 + w**2)
c     derive mach number = ------------------------
c                           sqrt(gamma * pres / rho)
c     
c     Inputs/Outputs
c     u            =>  (const)  state array
c     ulo,uhi      =>  (const)  index limits of u
c     nvar         =>  (const)  number of variables in state array
c     xlo          =>  (const)  physical location of lo end of u
c     xhi          =>  (const)  physical location of hi end of u
c     mach        <=>  (modify) array holding derived data
c     dlo,dhi      =>  (const)  index limits of mach array
c     ovlo,ovhi    =>  (const)  subregion where derivation is done
c
      subroutine FORT_DERMACH (u,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo,xhi,mach,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3,ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T vx, vy, vz
      REAL_T xlo(3), xhi(3)
      REAL_T u(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T mach(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)

      integer i,j,k
      REAL_T vvtmp, rho, c, eng, pres
#include "xxmeth.fh"	

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
		  rho = u(i,j,k,1)
                  vx = u(i,j,k,2)/rho
                  vy = u(i,j,k,3)/rho
                  vz = u(i,j,k,4)/rho

                  vvtmp = vx**2 + vy**2 + vz**2
		  eng = u(i,j,k,5) - half * rho * vvtmp
		  pres = eng * (gamma - one)
		  c = sqrt(gamma * pres / rho)
                  mach(i,j,k) = sqrt(vvtmp) / c
               end do
            end do
         end do
      end

       subroutine FORT_SELECTPOS (prob,DIMS(c),rn,dx,dy,dz,x,y,z)

       implicit none

       integer DIMDEC(c)
       REAL_T  prob(DIMV(c))
       REAL_T  dx,dy,dz,x,y,z,rn

       integer ix, iy, iz

       do iz = c_l3, c_h3-1
           if (rn.le.prob(c_h1,c_h2,iz)) go to 100
       enddo

100    continue

       do iy = c_l2, c_h2-1
           if (rn.le.prob(c_h1,iy,iz)) go to 200
       enddo

200    continue

       do ix = c_l1, c_h1-1
          if (rn.le.prob(ix,iy,iz)) go to 300
       enddo

300    continue

       x = (ix+half)*dx
       y = (iy+half)*dy
       z = (iz+half)*dz

       end

      subroutine FORT_PROBFAB (rf,DIMS(c),rr,prob,DIMS(p),
     &                         nspec,nreac,ispec,
     &                         edges,edgeslen,pedges,isrz,dx,dy,dz)

       implicit none

       integer DIMDEC(c), DIMDEC(p), nspec, nreac, ispec, isrz
       integer edgeslen, edges(0:edgeslen-1), pedges(0:nspec-1)
       REAL_T  rf(DIMV(c),0:nreac-1)
       REAL_T  rr(DIMV(c),0:nreac-1)
       REAL_T  prob(DIMV(p)), dx, dy, dz

       integer ix, iy, iz, ie, nedges, rxnid, factor, nu
       REAL_T  lambda, lmax

       if (ispec.lt.zero.or.ispec.gt.nspec-1) stop 'invalid ispec'

       nedges = edges(pedges(ispec))

       do iz = p_l3, p_h3
          do iy = p_l2, p_h2
             do ix = p_l1, p_h1

                lambda = 0

                do ie = 0, nedges-1

                   rxnid  = edges(pedges(ispec)+ie*4+1)
                   factor = edges(pedges(ispec)+ie*4+2)
                   nu     = edges(pedges(ispec)+ie*4+4)

                   if (factor.gt.zero) then
                      lambda = nu*lambda + rf(ix,iy,iz,rxnid)/factor
                   else
                      lambda = nu*lambda - rr(ix,iy,iz,rxnid)/factor
                   end if

                end do

                prob(ix,iy,iz) = lambda*dx*dy*dz

             end do
          end do
       enddo

       end

       subroutine FORT_ACCUMPROB (prob, DIMS(p), totreact, cumprob)

       implicit none

       integer DIMDEC(p)
       REAL_T  prob(DIMV(p)), totreact, cumprob

       integer ix, iy, iz

       do iz = p_l3, p_h3
          do iy = p_l2, p_h2
             do ix = p_l1, p_h1

                prob(ix,iy,iz) = cumprob+prob(ix,iy,iz)/totreact
                cumprob        = prob(ix,iy,iz)

             enddo
          enddo
       enddo

       end





      subroutine FORT_TEMPQMODEL (u,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo, xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3, ovlo,ovhi)
c  14
c  density
c  xmom
c  ymom
c  zmom
c  eden
c  tracer
c  fuel
c  oxidizer
c  products
c  pressure
c  xvel
c  yvel
c  zvel
c  vfrac

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T u(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)

      integer i,j,k
      REAL_T vx,vy,vz,eng

      REAL_T am,bm,cm,Rm,yf,ya,yp,Tm,Tmin
      REAL_T rgdnv, pgdnv, fuel_R, air_R, products_R

      call bl_abort('FORT_TEMPQMODEL not implemented.')
#if 0
c     do k = ovlo(3),ovhi(3) 
c         do j = ovlo(2),ovhi(2) 
c             do i = ovlo(1),ovhi(1) 
c                 vx = u(i,j,k,2)/u(i,j,k,1)
c                 vy = u(i,j,k,3)/u(i,j,k,1)
c                 vz = u(i,j,k,4)/u(i,j,k,1)
c                 eng = u(i,j,k,5)/u(i,j,k,1)
c                 dat(i,j,k) = eng - (vx**2+vy**2+vz**2)*half
c              end do
c           end do
c        end do

      data  Tmin/200.0/
      data  fuel_R/3.39e+06/
      data  air_R/2.87e+06/
      data  products_R/2.77e+06/


      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
		rgdnv = u(i,j,k,1)
		pgdnv = u(i,j,k,10)
                yf = u(i,j,k,7)/u(i,j,k,1)
                ya = u(i,j,k,8)/u(i,j,k,1)
                yp = u(i,j,k,9)/u(i,j,k,1)
CCCCCCCCCC      //ia = 1
CCCCCCCCCC      //am = yf*fuel_a + ya*air_a(ia) + yp*products_a
CCCCCCCCCC      //bm = yf*fuel_b + ya*air_b(ia) + yp*products_b
CCCCCCCCCC      //cm = yf*fuel_c + ya*air_c(ia) + yp*products_c
                Rm = yf*fuel_R + ya*air_R + yp*products_R
                Tm = pgdnv / (rgdnv*Rm)
                if(Tm .ge. air_T1)then
CCCCCCCCCC        //ia = 2
CCCCCCCCCC        //am = yf*fuel_a + ya*air_a(ia) + yp*products_a
CCCCCCCCCC        //bm = yf*fuel_b + ya*air_b(ia) + yp*products_b
CCCCCCCCCC        //cm = yf*fuel_c + ya*air_c(ia) + yp*products_c
                  Rm = yf*fuel_R + ya*air_R + yp*products_R
                  Tm = pgdnv / (rgdnv * Rm)
                endif
                Tm = max(Tm,Tmin)
                dat(i,j,k) = Tm
CCCCCCCCCC      egdnv(i,j,k) = am*Tm**2 + bm*Tm + cm
CCCCCCCCCC      egdnv(i,j,k) = rgdnv(i,j,k)*egdnv(i,j,k)

               end do
            end do
         end do
#endif

       end


      subroutine FORT_FUELFRAC (u,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo, xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3, ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T u(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)

      integer i,j,k

      REAL_T yf,ya,yp

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
                yf = u(i,j,k,7)/u(i,j,k,1)
CCCCC           ya = u(i,j,k,8)/u(i,j,k,1)
CCCCC           yp = u(i,j,k,9)/u(i,j,k,1)
                dat(i,j,k) = yf

               end do
            end do
         end do

      end


      subroutine FORT_OXFRAC (u,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo, xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3, ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T u(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)

      integer i,j,k

      REAL_T yf,ya,yp

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
CCCCC           yf = u(i,j,k,7)/u(i,j,k,1)
                ya = u(i,j,k,8)/u(i,j,k,1)
CCCCC           yp = u(i,j,k,9)/u(i,j,k,1)
                dat(i,j,k) = ya

               end do
            end do
         end do

      end


      subroutine FORT_PRODFRAC (u,uloi1,uloi2,uloi3,uhii1,uhii2,uhii3,
     $ nvar,xlo, xhi,dat,dloi1,dloi2,dloi3,dhii1,dhii2,dhii3, ovlo,ovhi)

      integer uloi1,uloi2,uloi3
      integer uhii1,uhii2,uhii3
      integer dloi1,dloi2,dloi3
      integer dhii1,dhii2,dhii3

      integer ovlo(3), ovhi(3)
      integer nvar
      REAL_T xlo(3), xhi(3)
      REAL_T u(uloi1:uhii1, uloi2:uhii2, uloi3:uhii3, nvar)
      REAL_T dat(dloi1:dhii1, dloi2:dhii2, dloi3:dhii3)

      integer i,j,k

      REAL_T yf,ya,yp

      do k = ovlo(3),ovhi(3) 
          do j = ovlo(2),ovhi(2) 
              do i = ovlo(1),ovhi(1) 
CCCCC           yf = u(i,j,k,7)/u(i,j,k,1)
CCCCC           ya = u(i,j,k,8)/u(i,j,k,1)
                yp = u(i,j,k,9)/u(i,j,k,1)
                dat(i,j,k) = yp

               end do
            end do
         end do

      end




c
c      KENG
c   
c      MASS = sum{ vol(i,j)*rho(i,j) }
c   
c      Inputs / Outputs:
c   
c      rho        => density field
c      rlo,rhi    => index limits of rho aray
c      lo,hi      => index limits of grid interior
c      delta      => cell size
c      sum        <=  total mass
c      tmp        => temp column array
c
       subroutine FORT_KENG (vels,DIMS(vels),DIMS(grid),delta,tmp)
       integer DIMDEC(vels)
       integer DIMDEC(grid)
       REAL_T  delta(3)
       REAL_T  vels(DIMV(vels), 3)
       REAL_T  tmp(DIM2(grid))

       integer i, j, k
       REAL_T  dx, dy, dz, vol

       dx = delta(1)
       dy = delta(2)
       dz = delta(3)
       vol = dx*dy*dz

       do k = ARG_L3(grid), ARG_H3(grid)
          tmp(k) = float(k)
       end do

       do k = ARG_L3(grid), ARG_H3(grid)
          do i = ARG_L1(grid), ARG_H1(grid)
             do j = ARG_L2(grid), ARG_H2(grid)
CCC             tmp(j) = tmp(j) + vol*rho(i,j,k)
             end do
          end do
       end do

      end

