c
c $Id: LO_2D.F,v 1.2 1999-02-03 00:16:54 lijewski Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include "LO_BCTYPES.H"
#include "LO_F.H"
#include "ArrayLim.H"

c-----------------------------------------------------------------------
      subroutine FORT_RESIDL (
     $     res, DIMS(res), 
     $     rhs, DIMS(rhs),
     $     phi, DIMS(phi),
     $     lo, hi, nc
     $     )
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(phi)
      REAL_T phi(DIMV(phi),nc)
      integer DIMDEC(rhs)
      REAL_T rhs(DIMV(rhs),nc)
      integer DIMDEC(res)
      REAL_T res(DIMV(res),nc)
c
      integer i
      integer j
      integer n
c
      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               res(i,j,n) = rhs(i,j,n) - phi(i,j,n)
            end do
         end do
      end do
c
      end
c-----------------------------------------------------------------------
      subroutine FORT_HARMONIC_AVERAGEEC (
     $     c, DIMS(c),
     $     f, DIMS(f),
     $     lo, hi, nc,
     $     cdir
     $     )
c
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cdir
      integer DIMDEC(f)
      REAL_T f(DIMV(f),nc)
      integer DIMDEC(c)
      REAL_T c(DIMV(c),nc)
c
      REAL_T factor, den
      parameter(factor=2.00)
      integer n
      integer i
      integer j
c
      if ( cdir .eq. 0 ) then
         do n = 1, nc
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)+1
                  den = f(2*i,2*j,n) + f(2*i,2*j+1,n)
                  if (den .ne. 0.0) then
                    c(i,j,n) =  factor*f(2*i,2*j,n)*f(2*i,2*j+1,n)/den
                  else
                    c(i,j,n) =  0.0
                  end if
               end do
            end do
         end do
      else if (cdir .eq. 1 ) then
         do n = 1, nc
            do j = lo(2), hi(2)+1
               do i = lo(1), hi(1)
                  den = f(2*i,2*j,n) + f(2*i+1,2*j,n)
                  if (den .ne. 0.0) then
                    c(i,j,n) =  factor*f(2*i,2*j,n)*f(2*i+1,2*j,n)/den
                  else
                    c(i,j,n) =  0.0
                  end if
               end do
            end do
         end do
      end if
c
      end
c-----------------------------------------------------------------------
      subroutine FORT_AVERAGEEC (
     $     c, DIMS(c),
     $     f, DIMS(f),
     $     lo, hi, nc,
     $     cdir
     $     )
c
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cdir
      integer DIMDEC(f)
      REAL_T f(DIMV(f),nc)
      integer DIMDEC(c)
      REAL_T c(DIMV(c),nc)
c
      integer n
      integer i
      integer j
      REAL_T denom
      parameter(denom=half)
c
      if (cdir .eq. 0 ) then
         do n = 1, nc
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)+1
                  c(i,j,n) = (f(2*i,2*j,n) + f(2*i,2*j+1,n))*denom
               end do
            end do
         end do
      else if (cdir .eq. 1) then
         do n = 1, nc
            do j = lo(2), hi(2)+1
               do i = lo(1), hi(1)
                  c(i,j,n) = (f(2*i,2*j,n) + f(2*i+1,2*j,n))*denom
               end do
            end do
         end do
      end if
c
      end
c-----------------------------------------------------------------------
      subroutine FORT_AVERAGECC (
     $     c, DIMS(c),
     $     f, DIMS(f),
     $     lo, hi, nc
     $     )
c
      integer nc
      integer DIMDEC(f)
      integer DIMDEC(c)
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)
c
      integer i
      integer j
      integer n
      REAL_T denom
      parameter(denom=fourth)
c
      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               c(i,j,n) =  (
     $              f(2*i+1,2*j+1,n) + f(2*i  ,2*j+1,n)
     $              + f(2*i+1,2*j  ,n) + f(2*i  ,2*j  ,n))*denom
            end do
         end do
      end do
c
      end
c-----------------------------------------------------------------------
      subroutine FORT_APPLYBC (
     $     flagden, flagbc, maxorder,
     $     phi,   DIMS(phi),
     $     cdir, bct, bcl,
     $     bcval, DIMS(bcval),
     $     mask,  DIMS(mask),
     $     den,   DIMS(den),
     $     lo, hi, nc,
     $     h
     $     )
c
c     If the boundary is of Neumann type, set the ghost cell value to
c     that of the outermost point in the valid data (2nd order accurate)
c     and then fill the "den" array with the value "1"
c     
c     
c     If flagbc==1:
c     
c     If the boundary is of Dirichlet type, construct a polynomial
c     interpolation through the boundary location and internal points
c     (at locations x(-1:len-2) that generates the ghost cell value (at
c     location xInt).  Then fill the ghost cell with the interpolated value.
c     If flagden==1, load the "den" array with the interpolation
c     coefficient corresponding to outermost point in the valid region
c     ( the coef(0) corresponding to the location x(0) )
c      
c     Note: 
c     The bc type = LO_REFLECT_ODD is a special type of boundary condition.
c      
      integer maxorder
      integer nc, cdir, flagden, flagbc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(phi)
      REAL_T phi(DIMV(phi),nc)
      integer DIMDEC(den)
      REAL_T den(DIMV(den))
      integer DIMDEC(bcval)
      REAL_T bcval(DIMV(bcval),nc)
      integer DIMDEC(mask)
      integer mask(DIMV(mask))
      integer bct
      REAL_T bcl
      REAL_T h(BL_SPACEDIM)
c
      integer i
      integer j
      integer n
      logical is_dirichlet
      logical is_neumann
      logical is_odd
      REAL_T xb
c
      integer lenx
      integer leny
      integer m
c
      integer Lmaxorder
      integer maxmaxorder
      parameter(maxmaxorder=4)
      REAL_T x(-1:maxmaxorder-2)
      REAL_T coef(-1:maxmaxorder-2)
      REAL_T xInt
c
      is_dirichlet(i) = ( i .eq. LO_DIRICHLET )
      is_neumann(i)   = ( i .eq. LO_NEUMANN )
c
      if ( maxorder .eq. -1 ) then
         Lmaxorder = maxmaxorder
      else
         Lmaxorder = MIN(maxorder,maxmaxorder)
      end if
      lenx = MIN(hi(1)-lo(1), Lmaxorder-2)
      leny = MIN(hi(2)-lo(2), Lmaxorder-2)
c
c     TODO:
c     In order for this to work with growing multigrid, must
c     sort xa[] because it is possible for the xb value to lay
c     within this range.
c     
c     The Left face of the grid
c
      if(cdir .eq. 0) then
         if (is_neumann(bct)) then
            do n = 1, nc
               do j = lo(2), hi(2)
                  phi(lo(1)-1,j,n) = cvmgt(
     $                 phi(lo(1),j,n),
     $                 phi(lo(1)-1,j,n),
     $                 mask(lo(1)-1,j) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1) then
               do j = lo(2), hi(2)
                  den(lo(1),j) = 1.0
               end do
            end if
         else if (is_dirichlet(bct)) then
            do m=0,lenx
               x(m) = m + 0.5
            end do
            x(-1) = - bcl/h(1)
            xInt = - 0.5
            call polyInterpCoeff(xInt, x, lenx+2, coef)
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do j = lo(2), hi(2)
                     phi(lo(1)-1, j, n) = cvmgt(
     $                    bcval(lo(1)-1,j,n)*coef(-1),
     $                    phi(lo(1)-1, j, n),
     $                    mask(lo(1)-1,j) .gt. 0)
                  end do
               else
                  do j = lo(2), hi(2)
                     phi(lo(1)-1, j, n) = cvmgt(
     $                    0.0,
     $                    phi(lo(1)-1, j, n),
     $                    mask(lo(1)-1,j) .gt. 0)
                  end do
               end if
               do m = 0, lenx
                  do j = lo(2), hi(2)
                     phi(lo(1)-1,j,n) = cvmgt(
     $                    phi(lo(1)-1,j,n)
     $                    + phi(lo(1)+m, j, n)*coef(m),
     $                    phi(lo(1)-1,j,n),
     $                    mask(lo(1)-1,j) .gt. 0)
                  end do
               end do
            end do
            if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  den(lo(1),j) = cvmgt(coef(0), 0.0,
     $                 mask(lo(1)-1,j) .gt. 0)
               end do
            end if

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do j = lo(2), hi(2)
                  phi(lo(1)-1, j, n) = cvmgt(
     $                -phi(lo(1),j,n),
     $                 phi(lo(1)-1, j, n),
     $                 mask(lo(1)-1,j) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  den(lo(1),j) = cvmgt(-1.0, 0.0,
     $                 mask(lo(1)-1,j) .gt. 0)
               end do
            end if

         else 
            print *,'UNKNOWN BC ON LEFT FACE IN APPLYBC'
            stop
         end if
      end if
c     
c     The Right face of the grid
c
      if(cdir .eq. 2) then
         if(is_neumann(bct)) then
            do n = 1, nc
               do j = lo(2), hi(2)
                  phi(hi(1)+1,j,n) = cvmgt(
     $                 phi(hi(1), j, n),
     $                 phi(hi(1)+1, j, n),
     $                 mask(hi(1)+1,j) .gt. 0)
               end do
            end do
	    if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  den(hi(1),j) = 1.0
               end do
	    end if
         else if (is_dirichlet(bct)) then
            do m=0,lenx
               x(m) = m + 0.5
            end do
            x(-1) = - bcl/h(1)
            xInt = - 0.5
            call polyInterpCoeff(xInt, x, lenx+2, coef)
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do j = lo(2), hi(2)
                     phi(hi(1)+1,j,n) = cvmgt(
     $                    bcval(hi(1)+1,j,n)*coef(-1),
     $                    phi(hi(1)+1,j,n),
     $                    mask(hi(1)+1,j) .gt. 0)
                  end do
               else
                  do j = lo(2), hi(2)
                     phi(hi(1)+1,j,n) = cvmgt(
     $                    0.0,
     $                    phi(hi(1)+1,j,n),
     $                    mask(hi(1)+1,j) .gt. 0)
                  end do
               end if
               do m = 0, lenx
                  do j = lo(2), hi(2)
                     phi(hi(1)+1,j,n) = cvmgt(
     $                    phi(hi(1)+1,j,n)
     $                    + phi(hi(1)-m,j,n)*coef(m),
     $                    phi(hi(1)+1,j,n),
     $                    mask(hi(1)+1,j) .gt. 0)
                  end do
               end do
            end do
            if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  den(hi(1),j)   = cvmgt(coef(0), 0.0,
     $                 mask(hi(1)+1,j) .gt. 0)
               end do
            end if

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do j = lo(2), hi(2)
                  phi(hi(1)+1, j, n) = cvmgt(
     $                -phi(hi(1),j,n),
     $                 phi(hi(1)+1, j, n),
     $                 mask(hi(1)+1,j) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               do j = lo(2), hi(2)
                  den(hi(1),j) = cvmgt(-1.0, 0.0,
     $                 mask(hi(1)+1,j) .gt. 0)
               end do
            end if

         else
            print *,'UNKNOWN BC ON RIGHT FACE IN APPLYBC'
            stop
         end if
      end if
c
c     The Bottom of the Grid
c
      if(cdir .eq. 1) then
         if(is_neumann(bct)) then
            do n = 1, nc
               do i = lo(1),hi(1)
                  phi(i,lo(2)-1,n) = cvmgt(
     $                 phi(i,lo(2),n),
     $                 phi(i,lo(2)-1,n),
     $                 mask(i,lo(2)-1) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1),hi(1)
                  den(i,lo(2))   = 1.0
               end do
            end if
         else if (is_dirichlet(bct)) then
            do m=0,leny
               x(m) = m + 0.5
            end do
            x(-1) = - bcl/h(2)
            xInt = - 0.5
            call polyInterpCoeff(xInt, x, leny+2, coef)
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do i = lo(1), hi(1)
                     phi(i,lo(2)-1,n) = cvmgt(
     $                    bcval(i,lo(2)-1,n)*coef(-1),
     $                    phi(i,lo(2)-1,n),
     $                    mask(i,lo(2)-1) .gt. 0)
                  end do
               else
                  do i = lo(1), hi(1)
                     phi(i,lo(2)-1,n) = cvmgt(
     $                    0.0,
     $                    phi(i,lo(2)-1,n),
     $                    mask(i,lo(2)-1) .gt. 0)
                  end do
               end if
               do m = 0, leny
                  do i = lo(1), hi(1)
                     phi(i, lo(2)-1, n) = cvmgt(
     $                    phi(i, lo(2)-1,n)
     $                    + phi(i, lo(2)+m,n)*coef(m),
     $                    phi(i, lo(2)-1, n),
     $                    mask(i, lo(2)-1) .gt. 0)
                  end do
               end do
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1), hi(1)
                  den(i, lo(2))   = cvmgt(coef(0), 0.0,
     $                 mask(i, lo(2)-1) .gt. 0)
               end do
            end if

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do i = lo(1), hi(1)
                  phi(i,lo(2)-1,n) = cvmgt(
     $                -phi(i,lo(2),n),
     $                 phi(i,lo(2)-1,n),
     $                 mask(i,lo(2)-1) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1), hi(1)
                  den(i,lo(2)) = cvmgt(-1.0, 0.0,
     $                 mask(i,lo(2)-1) .gt. 0)
               end do
            end if

         else
            print *,'UNKNOWN BC ON BOTTOM FACE IN APPLYBC'
            stop
         end if
      end if
c
c     The top of the grid
c
      if (cdir .eq. 3) then
         if(is_neumann(bct)) then
            do n = 1, nc
               do i = lo(1), hi(1)
                  phi(i,hi(2)+1,n) = cvmgt(
     $                 phi(i,hi(2),n),
     $                 phi(i,hi(2)+1,n),
     $                 mask(i,hi(2)+1) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1), hi(1)
                  den(i,hi(2))   = 1.0
               end do
            end if
         else if (is_dirichlet(bct)) then
            if ( bct .eq. LO_REFLECT_ODD ) leny = 0
            do m=0,leny
               x(m) = m + 0.5
            end do
            x(-1) = - bcl/h(2)
            xInt = - 0.5
            call polyInterpCoeff(xInt, x, leny+2, coef)
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  do i = lo(1), hi(1)
                     phi(i,hi(2)+1,n) = cvmgt(
     $                    bcval(i,hi(2)+1,n)*coef(-1),
     $                    phi(i,hi(2)+1,n),
     $                    mask(i,hi(2)+1) .gt. 0)
                  end do
               else
                  do i = lo(1), hi(1)
                     phi(i,hi(2)+1,n) = cvmgt(
     $                    0.0,
     $                    phi(i,hi(2)+1,n),
     $                    mask(i,hi(2)+1) .gt. 0)
                  end do
               end if
               do m = 0, leny
                  do i = lo(1), hi(1)
                     phi(i, hi(2)+1,n) = cvmgt(
     $                    phi(i,hi(2)+1,n)
     $                    + phi(i, hi(2)-m,n)*coef(m),
     $                    phi(i,hi(2)+1,n),
     $                    mask(i,hi(2)+1) .gt. 0)
                  end do
               end do
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1), hi(1)
                  den(i,hi(2))   = cvmgt(coef(0), 0.0,
     $                 mask(i,hi(2)+1) .gt. 0)
               end do
            end if

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               do i = lo(1), hi(1)
                  phi(i,hi(2)+1,n) = cvmgt(
     $                -phi(i,hi(2),n),
     $                 phi(i,hi(2)+1,n),
     $                 mask(i,hi(2)+1) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               do i = lo(1), hi(1)
                  den(i,hi(2)) = cvmgt(-1.0, 0.0,
     $                 mask(i,hi(2)+1) .gt. 0)
               end do
            end if

         else
            print *,'UNKNOWN BC ON TOP FACE IN APPLYBC'
            stop
         end if
      end if
c
      end
