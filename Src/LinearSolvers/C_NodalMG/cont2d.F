c General graphics routine for plotting contours given a matrix of values.
c
c by Louis H. Howell, 11/12/86
c
c subroutine cont(a,mdim,m,n,x1,x2,y1,y2,idir,value,ifail)
c
c   This routine takes a matrix a as an argument and plots the curves where
c the function defined by this matrix is equal to the input variable value.
c The matrix a contains values of a function of two variables on the nodes of
c a grid.  Linear interpolation between the values at these nodes is used to
c determine where contour lines cross the lines between nodes.
c
c   a:  This is a double precision matrix dimensioned a(0:mdim,0:*).  The
c portion actually plotted is the subset (0:m,0:n).
c
c   mdim,m,n,idir:  all integer input parameters
c
c   x1,x2,y1,y2:  all double precision input parameters
c   These parameters determine where on the plotting surface the contours
c will be drawn.  Unlike the surface graphics routines this subroutine does
c not set the dimensions or clipping limits on the plotting surface.  The
c dimensions must be set by the user before this subroutine is called; the
c contours are then plotted in the box defined by x1,x2,y1,y2.  Clipping is
c unnecessary since no contour lines will be drawn beyond these limits.
c
c   idir:  This parameter should be either 0 or 1, depending on which subscript
c of the matrix a should correspond to which axis of the plotting surface.
c Let the first subscript be i and the second be j, so an element of a is
c a(i,j).  If idir=0 i increases downward and j increases to the right, as in
c a matrix.  If idir=1 i increases to the right and j increases upward, as in
c the xy plane.  These directions are correct provided that x1<x2 and y1<y2.
c It is quite legal, however, to reverse one or both of these inequalities in
c order to produce an apparent reflection of the grid.
c
c   value:  This double precision input parameter is the function value for
c which the contour is to be plotted.  Normally the subroutine would be called
c several times with different values for this variable.  Contours may be any
c shape down to the resolution of the grid, and do not need to be connected.
c The plot will show a saddle point if and only if a contour crosses all four
c sides of a grid square.
c
c   ifail:  This integer output parameter is returned as 0 if the program
c runs normally.  It is returned as 1 if the requested contour value is exactly
c equal to the values of the function at two orthogonally adjacent grid points.
c In this case the contour is not well defined.  The program will not crash in
c this case, but the plot may show a split or multiple contour line at the
c offending points.
c
c   For demonstration purposes a demo program for this routine is included
c in /usr/local/tools/demo.  The sample makefile in that directory links
c to all the necessary libraries required for the program to compile.  Note
c that several other "tools" routines must be called to set up the plotting
c surface before this subroutine may be called, see /usr/local/tools/README
c for more information.  A single precision version of this subroutine, scont,
c is also available.
      subroutine cont2d(a,mlo,mhi,nlo,nhi,ilo,ihi,jlo,jhi,
     @                  x1,x2,y1,y2,idir,value,ifail)
      integer mlo, mhi, nlo, nhi, ilo, ihi, jlo, jhi, idir, ifail
      real*8 x1, x2, y1, y2, value
      real*8 a(mlo:mhi,nlo:nhi)
      logical left, right, top, bottom, btwn
      real*8 xstep, ystep, lt, lb, rt, rb
      real*8 xl, xr, xt, xb, yl, yr, yt, yb
      integer i, j
      ifail=0
      if (idir .eq. 0)then
         xstep=(x2-x1)/(jhi-jlo)
         ystep=(y2-y1)/(ihi-ilo)
      else
         xstep=(x2-x1)/(ihi-ilo)
         ystep=(y2-y1)/(jhi-jlo)
      endif
      do 100 i=ilo+1,ihi
         do 100 j=jlo+1,jhi
            if (idir .eq. 0)then
               lt=a(i-1,j-1)
               lb=a(i,j-1)
               rt=a(i-1,j)
               rb=a(i,j)
               xl=x1+(j-1)*xstep
               yb=y1+(ihi-i)*ystep
            else
               lt=a(i-1,j)
               lb=a(i-1,j-1)
               rt=a(i,j)
               rb=a(i,j-1)
               xl=x1+(i-ilo-1)*xstep
               yb=y1+(j-jlo-1)*ystep
            endif
            left=btwn(lb,value,lt)
            if (left)then
               if (lb .ne. lt)then
                  yl=yb+(value-lb)/(lt-lb)*ystep
               else
                  yl=yb
                  ifail=1
               endif
            endif
            right=btwn(rb,value,rt)
            if (right)then
               xr=xl+xstep
               if (rb .ne. rt)then
                  yr=yb+(value-rb)/(rt-rb)*ystep
               else
                  yr=yb+ystep
                  ifail=1
               endif
            endif
            top=btwn(lt,value,rt)
            if (top)then
               yt=yb+ystep
               if (lt .ne. rt)then
                  xt=xl+(value-lt)/(rt-lt)*xstep
               else
                  xt=xl
                  ifail=1
               endif
            endif
            bottom=btwn(lb,value,rb)
            if (bottom)then
               if (lb .ne. rb)then
                  xb=xl+(value-lb)/(rb-lb)*xstep
               else
                  xb=xl+xstep
                  ifail=1
               endif
            endif
            if (left .and. right .and. top .and. bottom)then
               call pm(xl,yl)
               call pc(xr,yr)
               call pm(xt,yt)
               call pc(xb,yb)
            else if (top .and. bottom)then
               call pm(xt,yt)
               call pc(xb,yb)
            else if (left)then
               call pm(xl,yl)
               if (right)then
                  call pc(xr,yr)
               else if (top)then
                  call pc(xt,yt)
               else
                  call pc(xb,yb)
               endif
            else if (right)then
               call pm(xr,yr)
               if (top)then
                  call pc(xt,yt)
               else
                  call pc(xb,yb)
               endif
            endif
 100        continue
      return
      end
