c
c 3D graphics routines for plotting contours given a matrix of values.
c
c Louis H. Howell, last change (slices) 9/8/94
c
c Comments?  Try the 2D version, it may help a little.
c
c Slice stuff uses 2D version:

#include "cont2d.F"

c Routines to cut slice:

      subroutine getsl0(a,mlo,mhi,nlo,nhi,olo,ohi,
     @  ilo,ihi,jlo,jhi,klo,khi,b,islice)
      implicit real*8 (a-h,l,p-z), integer(i-k,m-o)
      real*8 a(mlo:mhi,nlo:nhi,olo:ohi)
      real*8 b(jlo:jhi,klo:khi)
      do 10 k = klo, khi
         do 10 j = jlo, jhi
 10         b(j,k) = a(islice,j,k)
      end

      subroutine getsl1(a,mlo,mhi,nlo,nhi,olo,ohi,
     @  ilo,ihi,jlo,jhi,klo,khi,b,jslice)
      implicit real*8 (a-h,l,p-z), integer(i-k,m-o)
      real*8 a(mlo:mhi,nlo:nhi,olo:ohi)
      real*8 b(ilo:ihi,klo:khi)
      do 10 k = klo, khi
         do 10 i = ilo, ihi
 10         b(i,k) = a(i,jslice,k)
      end

      subroutine getsl2(a,mlo,mhi,nlo,nhi,olo,ohi,
     @  ilo,ihi,jlo,jhi,klo,khi,b,kslice)
      implicit real*8 (a-h,l,p-z), integer(i-k,m-o)
      real*8 a(mlo:mhi,nlo:nhi,olo:ohi)
      real*8 b(ilo:ihi,jlo:jhi)
      do 10 j = jlo, jhi
         do 10 i = ilo, ihi
 10         b(i,j) = a(i,j,kslice)
      end

c Begin full 3D version:

      subroutine cont3d(a,mlo,mhi,nlo,nhi,olo,ohi,
     @  ilo,ihi,jlo,jhi,klo,khi,rx,ry,rz,
     @  xo,yo,zo,theta,phi,value,ifail)
      implicit real*8 (a-h,l,p-z), integer(i-k,m-o)
      real*8 a(mlo:mhi,nlo:nhi,olo:ohi)
      integer rx, ry, rz
      parameter(pi=3.14159265358979323846d0)
      ifail=0
      ct = cos(theta*pi/180)
      st = sin(theta*pi/180)
      cp = cos(phi*pi/180)
      sp = sin(phi*pi/180)
      ctsp = ct * sp
      stsp = st * sp
      hx = 1.d0 / rx
      hy = 1.d0 / ry
      hz = 1.d0 / rz
      do 10 k=klo,khi
         zp = k * hz - zo
         do 10 j=jlo,jhi-1
            yp = j * hy - yo
            do 10 i=ilo,ihi-1
               xp = i * hx - xo
               xlb = xp * ct - yp * st
               xlt = xlb - st * hy
               xrb = xlb + ct * hx
               xrt = xrb - st * hy
               ylb = xp * stsp + yp * ctsp + zp * cp
               ylt = ylb + ctsp * hy
               yrb = ylb + stsp * hx
               yrt = yrb + ctsp * hy
               lb=a(i,j,k)
               lt=a(i,j+1,k)
               rb=a(i+1,j,k)
               rt=a(i+1,j+1,k)
               call face(xlb,xlt,xrb,xrt,ylb,ylt,yrb,yrt,
     @                   lb,lt,rb,rt,value,ifail)
 10            continue
      do 20 j=jlo,jhi
         yp = j * hy - yo
         do 20 k=klo,khi-1
            zp = k * hz - zo
            do 20 i=ilo,ihi-1
               xp = i * hx - xo
               xlb = xp * ct - yp * st
               xlt = xlb
               xrb = xlb + ct * hx
               xrt = xrb
               ylb = xp * stsp + yp * ctsp + zp * cp
               ylt = ylb + cp * hz
               yrb = ylb + stsp * hx
               yrt = yrb + cp * hz
               lb=a(i,j,k)
               lt=a(i,j,k+1)
               rb=a(i+1,j,k)
               rt=a(i+1,j,k+1)
               call face(xlb,xlt,xrb,xrt,ylb,ylt,yrb,yrt,
     @                   lb,lt,rb,rt,value,ifail)
 20            continue
      do 30 i=ilo,ihi
         xp = i * hx - xo
         do 30 k=klo,khi-1
            zp = k * hz - zo
            do 30 j=jlo,jhi-1
               yp = j * hy - yo
               xlb = xp * ct - yp * st
               xlt = xlb
               xrb = xlb - st * hy
               xrt = xrb
               ylb = xp * stsp + yp * ctsp + zp * cp
               ylt = ylb + cp * hz
               yrb = ylb + ctsp * hy
               yrt = yrb + cp * hz
               lb=a(i,j,k)
               lt=a(i,j,k+1)
               rb=a(i,j+1,k)
               rt=a(i,j+1,k+1)
               call face(xlb,xlt,xrb,xrt,ylb,ylt,yrb,yrt,
     @                   lb,lt,rb,rt,value,ifail)
 30            continue
      return
      end

      subroutine face(xlb,xlt,xrb,xrt,ylb,ylt,yrb,yrt,
     @                lb,lt,rb,rt,value,ifail)
      implicit real*8 (a-h,l,p-z), integer(i-k,m-o)
      logical left,right,top,bottom,btwn
      left=btwn(lb,value,lt)
      if (left)then
         if (lb .ne. lt)then
            xl=xlb + (xlt - xlb) * (value-lb) / (lt-lb)
            yl=ylb + (ylt - ylb) * (value-lb) / (lt-lb)
         else
            xl=xlb
            yl=ylb
            ifail=1
         endif
      endif
      right=btwn(rb,value,rt)
      if (right)then
         if (rb .ne. rt)then
            xr=xrb + (xrt - xrb) * (value-rb) / (rt-rb)
            yr=yrb + (yrt - yrb) * (value-rb) / (rt-rb)
         else
            xr = xrb
            yr = yrb
            ifail=1
         endif
      endif
      top=btwn(lt,value,rt)
      if (top)then
         if (lt .ne. rt)then
            xt=xlt + (xrt - xlt) * (value-lt) / (rt-lt)
            yt=ylt + (yrt - ylt) * (value-lt) / (rt-lt)
         else
            xt=xlt
            yt=ylt
            ifail=1
         endif
      endif
      bottom=btwn(lb,value,rb)
      if (bottom)then
         if (lb .ne. rb)then
            xb=xlb + (xrb - xlb) * (value-lb) / (rb-lb)
            yb=ylb + (yrb - ylb) * (value-lb) / (rb-lb)
         else
            xb=xlb
            yb=ylb
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
      return
      end
