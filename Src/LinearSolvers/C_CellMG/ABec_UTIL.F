c
c $Id: ABec_UTIL.F,v 1.2 1998-07-14 20:23:41 car Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>

c-----------------------------------------------------------------------
c
c     Tridiagonal solve
c
      subroutine tridiag(a,b,c,r,u,n)
 
      integer n
      integer nmax
  
      REAL_T a(n)
      REAL_T b(n)
      REAL_T c(n)
      REAL_T r(n)
      REAL_T u(n)

      parameter (nmax = 4098)

      integer j
      REAL_T bet
      REAL_T gam(nmax)
      if (b(1) .eq. 0) print *,'CANT HAVE B(1) = ZERO'

      bet = b(1)
      u(1) = r(1)/bet

      do j = 2,n
        gam(j) = c(j-1)/bet
        bet = b(j) - a(j)*gam(j)
        if (bet .eq. 0) then
          print *,'TRIDIAG FAILED '
          stop
        end if
        u(j) = (r(j)-a(j)*u(j-1))/bet
      end do

      do j = n-1,1,-1
        u(j) = u(j) - gam(j+1)*u(j+1)
      end do

      return
      end

      
