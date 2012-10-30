#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "DarcySNES_F.H"
#include "ArrayLim.H"

#define SDIM 2

      subroutine FORT_RS_PDOTRES(Res,DIMS(Res),RSn,DIMS(RSn),
     &     RSnp1,DIMS(RSnp1),dt,phi,lo,hi,nc)
      implicit none      
      integer DIMDEC(Res)
      integer DIMDEC(RSn)
      integer DIMDEC(RSnp1)
      integer lo(SDIM),hi(SDIM)
      integer nc
      REAL_T  Res(DIMV(Res),nc)
      REAL_T  RSn(DIMV(RSn),nc)
      REAL_T  RSnp1(DIMV(RSnp1),nc)
      REAL_T  dt, phi
      
      integer i, j, n

c     On input, Res contains Div(rho.DarcyVel)
      do n=1,nc
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               Res(i,j,n) = (RSnp1(i,j,n) - RSn(i,j,n))*phi/dt + Res(i,j,n)
            enddo
         enddo
      enddo
      end


      subroutine FORT_DARCYFLUX(p,DIMS(p),kr,DIMS(kr),idx,DIMS(idx),fluxx, fluxy, 
     &     fluxx_lo1, fluxx_lo2, fluxy_lo1, fluxy_lo2, 
     &     fluxx_hi1, fluxx_hi2, fluxy_hi1, fluxy_hi2,
     &     lo, hi, kappax, kappaz, rho, mu, g, dx)

      implicit none
      integer DIMDEC(p)
      integer DIMDEC(kr)
      integer DIMDEC(idx)
      integer lo(SDIM),hi(SDIM)
      integer fluxx_lo1, fluxx_lo2, fluxy_lo1, fluxy_lo2, fluxx_hi1, fluxx_hi2, fluxy_hi1, fluxy_hi2
      REAL_T  p(DIMV(p))
      REAL_T  kr(DIMV(kr))
      integer idx(DIMV(idx))
      REAL_T  fluxx(fluxx_lo1:fluxx_hi1,fluxx_lo2:fluxx_hi2)
      REAL_T  fluxy(fluxy_lo1:fluxy_hi1,fluxy_lo2:fluxy_hi2)
      REAL_T  kappax, kappaz, rho, mu, g(2), dx(2)
      REAL_T  F, dxInv, dyInv, kr_upwind, lamx, lamy

      integer i, j
c     NOTE: If the cc on either side of a face is covered (idx<0), we do not compute the
c           Darcy flux there.  If this is a crse edge at a crse-fine boundary, the flux
c           will be filled by averaging down the fine values.  Otherwise, the flux at 
c           such a location will never be used.

c     X-edges
      dxInv = 1.d0 / dx(1)
      lamx = kappax / mu
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            if (idx(i-1,j).ge.0  .and.  idx(i,j).ge.0) then
               F = - ( (p(i,j)-p(i-1,j))*dxInv + rho*g(1) )
               if (F.gt.0.d0) then
                  kr_upwind = kr(i-1,j)
               else
                  kr_upwind = kr(i,j)
               endif
               fluxx(i,j) = F * kr_upwind * lamx * rho
            endif
         enddo
      enddo

c     Y-edges
      dyInv = 1.d0 / dx(2)
      lamy = kappaz / mu
      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            if (idx(i,j-1).ge.0  .and.  idx(i,j).ge.0) then
               F = - ( (p(i,j)-p(i,j-1))*dyInv + rho*g(2) )
               if (F.gt.0.d0) then
                  kr_upwind = kr(i,j-1)
               else
                  kr_upwind = kr(i,j)
               endif
               fluxy(i,j) = F * kr_upwind * lamy * rho
            endif
         enddo
      enddo
      end

