c
c $Id: LO_UTIL.F,v 1.3 2002-10-19 18:33:57 car Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif
      
#undef TESTING_POLYNOMINTERPCOEFFS
#ifdef TESTING_POLYNOMINTERPCOEFFS
#define NORDER 3
#define REAL_T real*8
#define zero 0.d0
#define one 1.d0
#else
#include <CONSTANTS.H>
#include <REAL.H>
#endif
    
      

c     polyInterpCoeff:
c  
c     This routine returns the Lagrange interpolating coefficients for a
c     polynomial through N points, evaluated at xInt (see Numerical Recipes,
c     v2, p102, e.g.):
c
c            (x-x2)(x-x3)...(x-xN)              (x-x1)(x-x2)...(x-x(N-1))
c    P(x) = ----------------------- y1  + ... + ------------------------  yN
c           (x1-x2)(x1-x3)...(x1-xN)            (x1-x2)(x1-x3)...(x1-xN)
c
c     P(xInt) = sum_(i=1)^(N) y[i]*c[i]
c
      subroutine polyInterpCoeff(xInt, x, N, c)
      implicit none
      integer N, i, j
      REAL_T xInt, x(N), c(N), num, den
      do j=1,N
         num = one
         den = one
         do i = 1,j-1
            num = num*(xInt - x(i))
            den = den*(x(j) - x(i))
         end do
         do i = j+1,N
            num = num*(xInt - x(i))
            den = den*(x(j) - x(i))
         end do
#ifndef NDEBUG
         if (den .eq. zero) STOP 'polyInterpCoeff::invalid data'
#endif         
         c(j) = num/den
      end do
      return
      end

      

#ifdef TESTING_POLYNOMINTERPCOEFFS

c     
c     This is a test driver for the routine polyInterpCoeff.  Sample data
c     is created from the statement function, and the location of the 
c     boundary node and internal nodes are set, as apporpriate.  The
c     number of points created is equal to the test NORDER set at the
c     top of this file through a define.  The coefficients are computed,
c     and then the ghost cell value is constructed from the resulting
c     coefficients and written out.
c

      program polyInterpCoeffTest
      integer i, j
      REAL_T c(NORDER), y(NORDER), yInt
      REAL_T x(NORDER), const, yfunc, xx, xInt
c
c     Make data, locations
      yfunc(xx) = (2*xx-one)**(NORDER-1) + const
      const = -one
      x(1) = zero
      y(1) = yfunc(x(1))
      do j=2,NORDER
         x(j) = j - 1.5D0*one
         y(j) = yfunc(x(j))
      end do
c      
c     Set interpolation point to ghost cell location
      xInt = (x(2) - x(3))*0.5D0
c      
c     Call routine
      call polyInterpCoeff(xInt,x,NORDER,c)
c
c     Evaluate result
      yInt = zero
      do j=1,NORDER
         yInt = yInt + c(j)*y(j)
      end do
c
c     Dump output
      write(6,*) 'x = ',x
      write(6,*) 'y = ',y
      write(6,*) 'c = ',c
      write(6,*) 'Interpolated y = ',yInt
      end

#endif
