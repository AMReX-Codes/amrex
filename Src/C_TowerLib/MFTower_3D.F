#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "MFTower_F.H"
#include "ArrayLim.H"

#define SDIM 3

      subroutine FORT_CC_TO_EC_GRAD(grad,DIMS(grad),
     &     c,DIMS(c),lo,hi,dx,mult,dir,nc)
      implicit none      
      integer DIMDEC(grad)
      integer DIMDEC(c)
      integer lo(SDIM),hi(SDIM)
      REAL_T  dx(0:SDIM-1)
      REAL_T  grad(DIMV(grad),0:nc-1)
      REAL_T  c(DIMV(c),0:nc-1)
      REAL_T  mult
      integer dir, nc
      
      integer i, j, k, n
      REAL_T tmult
      tmult = mult/dx(dir)
      if (dir.eq.0) then
         do n=0,nc-1
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)+1
                     grad(i,j,k,n) = tmult*(c(i,j,k,n) - c(i-1,j,k,n))
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=0,nc-1
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)+1
                  do i = lo(1),hi(1)
                     grad(i,j,k,n) = tmult*(c(i,j,k,n) - c(i,j-1,k,n))
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.2) then
         do n=0,nc-1
            do k = lo(3),hi(3)+1
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     grad(i,j,k,n) = tmult*(c(i,j,k,n) - c(i,j,k-1,n))
                  enddo
               enddo
            enddo
         enddo
      else
         print *, 'FORT_CC_TO_EC_GRAD::bad dir: ',dir
         call bl_pd_abort()
      endif
      end     

      subroutine FORT_EC_TO_CC_DIV(div,DIMS(div),
     &     e,DIMS(e),a,DIMS(a),v,DIMS(v),lo,hi,mult,dir,nc)
      implicit none      
      integer DIMDEC(div)
      integer DIMDEC(e)
      integer DIMDEC(a)
      integer DIMDEC(v)
      integer lo(SDIM),hi(SDIM)
      REAL_T  div(DIMV(div),0:nc-1)
      REAL_T  e(DIMV(e),0:nc-1)
      REAL_T  a(DIMV(a))
      REAL_T  v(DIMV(v))
      REAL_T  mult
      integer dir, nc
      
      integer i, j, k, n
      if (dir.eq.0) then
         do n=0,nc-1
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     div(i,j,k,n) = div(i,j,kn) + mult*(e(i+1,j,k,n)*a(i+1,j,k) - e(i,j,k,n)*a(i,j,k))/v(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=0,nc-1
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     div(i,j,k,n) = div(i,j,k,n) + mult*(e(i,j+1,k,n)*a(i,j+1,k) - e(i,j,k,n)*a(i,j,k))/v(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.2) then
         do n=0,nc-1
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     div(i,j,k,n) = div(i,j,k,n) + mult*(e(i,j,k+1,n)*a(i,j,k+1) - e(i,j,k,n)*a(i,j,k))/v(i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else
         print *, 'FORT_EC_TO_CC_DIV::bad dir: ',dir
         call bl_pd_abort()
      endif
      end     

      subroutine FORT_COARSEN_EC(ce,DIMS(ce),
     &     fe,DIMS(fe),celo,cehi,ratio,dir,nc)
      implicit none      
      integer ratio(SDIM), dir, nc
      integer DIMDEC(ce)
      integer DIMDEC(fe)
      integer celo(SDIM),cehi(SDIM)
      REAL_T  ce(DIMV(ce),nc)
      REAL_T  fe(DIMV(fe),nc)
      
      integer i, j, k, ii, jj, kk, n, rim1, rjm1, rkm1, ir, jr, kr
      REAL_T areaInv

      rim1 = ratio(1) - 1
      rjm1 = ratio(2) - 1
      rkm1 = ratio(3) - 1
      
      if (dir.eq.0) then
         areaInv = 1.d0/(ratio(2)*ratio(3))
      else if (dir.eq.1) then
         areaInv = 1.d0/(ratio(1)*ratio(3))
      else if (dir.eq.2) then
         areaInv = 1.d0/(ratio(1)*ratio(2))
      else
         print *, 'FORT_COARSEN_EC::bad dir: ',dir
         call bl_pd_abort()
      endif

      if (dir.eq.0) then
         do n=1,nc
            do jr=0,rjm1
               do kr=0,rkm1
                  do k = celo(3),cehi(3)
                     kk = ratio(3)*k + kr
                     do j = celo(2),cehi(2)
                        jj = ratio(2)*j + jr
                        do i = celo(1),cehi(1)
                           ii = ratio(1)*i
                           ce(i,j,k,n) = ce(i,j,k,n) + fe(ii,jj,kk,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=1,nc
            do ir=0,rim1
               do kr=0,rkm1
                  do k = celo(3),cehi(3)
                     kk = ratio(3)*k + kr
                     do i = celo(1),cehi(1)
                        ii = ratio(1)*i + ir
                        do j = celo(2),cehi(2)
                           jj = ratio(2)*j
                           ce(i,j,k,n) = ce(i,j,k,n) + fe(ii,jj,kk,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.2) then
         do n=1,nc
            do ir=0,rim1
               do jr=0,rjm1
                  do i = celo(1),cehi(1)
                     ii = ratio(1)*i + ir
                     do j = celo(2),cehi(2)
                        jj = ratio(2)*j + jr
                        do k = celo(3),cehi(3)
                           kk = ratio(3)*k
                           ce(i,j,k,n) = ce(i,j,k,n) + fe(ii,jj,kk,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
      
      do n=1,nc
         do k = celo(3),cehi(3)
            do j = celo(2),cehi(2)
               do i = celo(1),cehi(1)
                  ce(i,j,k,n) = ce(i,j,k,n) * areaInv
               enddo
            enddo
         enddo
      enddo

      end     

