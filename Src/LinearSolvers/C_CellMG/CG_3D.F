c     $Id: CG_3D.F,v 1.3 1998-07-14 20:23:42 car Exp $
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include "CG_F.H"
#include "ArrayLim.H"

c-----------------------------------------------------------------------
c      
c     CGUPDATE: Modify the input arrays as follows:
c     
c     phi = phi + a*pp
c     rr  = rr  - a*ww    through range lo:hi
c     
c     phi <=>
c     rr  <=>
c     a   <=
c     ww  <=
c     pp  <=
c     
c-----------------------------------------------------------------------
      subroutine FORT_CGUPDATE(
     $     phi, DIMS(phi),
     $     rr,  DIMS(rr),
     $     a,
     $     ww,  DIMS(ww),
     $     pp,  DIMS(pp),
     $     lo, hi, nc
     $     )
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(rr)
      REAL_T rr(DIMV(rr),nc)
      integer DIMDEC(ww)
      REAL_T ww(DIMV(ww),nc)
      integer DIMDEC(pp)
      REAL_T pp(DIMV(pp),nc)
      integer DIMDEC(phi)
      REAL_T phi(DIMV(phi),nc)
      REAL_T a
c ::: wyc change -- why was a a vector?
c     REAL_T a(nc)
c
      integer i,j,k,n
c
      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  phi(i,j,k,n) = phi(i,j,k,n) + a * pp(i,j,k,n)
                  rr(i,j,k,n) = rr(i,j,k,n) - a * ww(i,j,k,n)
c ::: wyc change -- why was a a vector?
c                 phi(i,j,k,n) = phi(i,j,k,n) + a(n) * pp(i,j,k,n)
c                 rr(i,j,k,n) = rr(i,j,k,n) - a(n) * ww(i,j,k,n)
               end do
            end do
         end do
      end do
c
      end
      
c-----------------------------------------------------------------------
c      
c     CGADVCP: Modify the input arrays as follows:
c     
c     pp = rr + b*pp    through range lo:hi
c     
c     pp  <=>
c     rr  <=>
c     b   <=
c     
c-----------------------------------------------------------------------
      subroutine FORT_CGADVCP(
     $     pp, DIMS(pp),
     $     rr, DIMS(rr),
     $     b,
     $     lo, hi, nc
     $     )
c
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(rr)
      REAL_T rr(DIMV(rr),nc)
      integer DIMDEC(pp)
      REAL_T pp(DIMV(pp),nc)
      REAL_T b
c ::: wyc change
c     REAL_T b(nc)
c
c
      integer i,j,k,n
c
      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  pp(i,j,k,n) = rr(i,j,k,n) + b*pp(i,j,k,n)
c                 pp(i,j,k,n) = rr(i,j,k,n) + b(n)*pp(i,j,k,n)
               end do
            end do
         end do
      end do
c
      end
      
c-----------------------------------------------------------------------
c      
c     CGAXP: Modify the input as follows:
c     
c     pw = Transpose(pp) . ww    through range lo:hi
c     
c     pw   =>
c     pp  <=
c     ww  <=
c     
c-----------------------------------------------------------------------
      subroutine FORT_CGXDOTY(
     $     pw,
     $     pp, DIMS(pp),
     $     ww, DIMS(ww),
     $     lo, hi, nc
     $     )
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(ww)
      REAL_T ww(DIMV(ww),nc)
      integer DIMDEC(pp)
      REAL_T pp(DIMV(pp),nc)
      REAL_T pw(nc)
c
      integer i, j, k, n
c
      do n = 1, nc
         pw(n) = 0.0
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  pw(n) = pw(n) + pp(i,j,k,n)*ww(i,j,k,n)
               end do
            end do
         end do
      end do
c
      end
