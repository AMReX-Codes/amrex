
#include "hg_version.H"

#define dims(a) a l0, a h0, a l1, a h1
#define dimdec(a) a l0:a h0,a l1:a h1
c#define dimdecu(a) (a h0-a l0+1)*(a h1-a l1+1)
#define dimdecu(a) *

      subroutine hggrad(gpx, gpy, dims(gp),
     @ dest, dims(dest),
     @ dims(freg),
     @ hx, hy)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(gp)
      integer dims(dest)
      integer dims(freg)
      real*8 gpx(dimdec(gp))
      real*8 gpy(dimdec(gp))
      real*8 dest(dimdec(dest))
      hxm1h = 0.5d0 / hx
      hym1h = 0.5d0 / hy
      do 10 j = fregl1, fregh1
         do 10 i = fregl0, fregh0
            gpx(i,j) = hxm1h * (dest(i+1,j) + dest(i+1,j+1) -
     @                          dest(i,j) - dest(i,j+1))
 10         gpy(i,j) = hym1h * (dest(i,j+1) + dest(i+1,j+1) -
     @                          dest(i,j) - dest(i+1,j))
      end

      subroutine hgdiv(src, dims(src),
     @ uf, vf, dims(f),
     @ dims(freg),
     @ hx, hy)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(src)
      integer dims(f)
      integer dims(freg)
      real*8 src(dimdec(src))
      real*8 uf(dimdec(f))
      real*8 vf(dimdec(f))
      real*8 hx, hy
      hxm1 = 1.d0 / hx
      hym1 = 1.d0 / hy
      fac = 0.5d0
      do 10 j = fregl1, fregh1
         do 10 i = fregl0, fregh0
 10         src(i,j) = fac *
     @        (hxm1 * (uf(i,j-1) - uf(i-1,j-1) +
     @                 uf(i,j)   - uf(i-1,j)) +
     @         hym1 * (vf(i-1,j) - vf(i-1,j-1) +
     @                 vf(i,j)   - vf(i,j-1)))
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgfdiv(src, dims(src),
     @ uc, vc, dims(c),
     @ uf, vf, dims(f),
     @ dims(creg),
     @ hx, hy, irat, idim, idir)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 uc(dimdec(c))
      real*8 vc(dimdec(c))
      real*8 uf(dimdec(f))
      real*8 vf(dimdec(f))
      real*8 hx, hy
      integer irat, idim, idir
      hxm1 = 1.d0 / hx
      hym1 = 1.d0 / hy
      fac0 = 1.d0 / (irat+1)
      if (idim .eq. 0) then
         i = cregl0
         if (idir .eq. 1) then
            iuc = i - 1
            iuf = i * irat
         else
            iuc = i
            iuf = i * irat - 1
         endif
         do 10 j = cregl1, cregh1
 10         src(i*irat,j*irat) = fac0 *
     @        (hym1 * (vc(iuc,j) - vc(iuc,j-1)) -
     @         hxm1 * idir * (uc(iuc,j) + uc(iuc,j-1)))
         fac0 = fac0 / (irat*irat)
         i = i * irat
         do 20 n = 0, irat-1
            fac1 = (irat-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do 20 j = irat*cregl1, irat*cregh1, irat
 20            src(i,j) = src(i,j) + fac1 *
     @           (hxm1 * idir * (uf(iuf,j-n) + uf(iuf,j-n-1) +
     @                           uf(iuf,j+n) + uf(iuf,j+n-1)) +
     @            hym1 * (vf(iuf,j-n) - vf(iuf,j-n-1) +
     @                    vf(iuf,j+n) - vf(iuf,j+n-1)))
      else
         j = cregl1
         if (idir .eq. 1) then
            juc = j - 1
            juf = j * irat
         else
            juc = j
            juf = j * irat - 1
         endif
         do 30 i = cregl0, cregh0
 30         src(i*irat,j*irat) = fac0 *
     @        (hxm1 * (uc(i,juc) - uc(i-1,juc)) -
     @         hym1 * idir * (vc(i,juc) + vc(i-1,juc)))
         fac0 = fac0 / (irat*irat)
         j = j * irat
         do 40 m = 0, irat-1
            fac1 = (irat-m) * fac0
            if (m .eq. 0) fac1 = 0.5d0 * fac1
            do 40 i = irat*cregl0, irat*cregh0, irat
 40            src(i,j) = src(i,j) + fac1 *
     @           (hxm1 * (uf(i-m,juf) - uf(i-m-1,juf) +
     @                    uf(i+m,juf) - uf(i+m-1,juf)) +
     @            hym1 * idir * (vf(i-m,juf) + vf(i-m-1,juf) +
     @                           vf(i+m,juf) + vf(i+m-1,juf)))
      endif
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgodiv(src, dims(src),
     @ uc, vc, dims(c),
     @ uf, vf, dims(f),
     @ dims(creg),
     @ hx, hy, irat, idir, jdir)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 uc(dimdec(c))
      real*8 vc(dimdec(c))
      real*8 uf(dimdec(f))
      real*8 vf(dimdec(f))
      real*8 hx, hy
      integer irat, idir, jdir
      hxm1 = 1.d0 / hx
      hym1 = 1.d0 / hy
      fac0 = dble(irat) / ((1.5d0*irat+1.d0)*irat-0.5d0)
      i = cregl0
      j = cregl1
      if (idir .eq. 1) then
         iuc = i - 1
         iuf = i * irat
      else
         iuc = i
         iuf = i * irat - 1
      endif
      if (jdir .eq. 1) then
         juc = j - 1
         juf = j * irat
      else
         juc = j
         juf = j * irat - 1
      endif
      src(i*irat,j*irat) = fac0 *
     @  (hxm1 * idir * (uc(iuc+idir,juc) - uc(iuc,juc) -
     @                  uc(iuc,juc+jdir)) +
     @   hym1 * jdir * (vc(iuc,juc+jdir) - vc(iuc,juc) -
     @                  vc(iuc+idir,juc)))
      fac0 = fac0 / irat
      i = i * irat
      j = j * irat
      src(i,j) = src(i,j) + fac0 *
     @  (hxm1 * idir * uf(iuf,juf) + hym1 * jdir * vf(iuf,juf))
      fac0 = fac0 / irat
      do 10 n = jdir, jdir*(irat-1), jdir
         fac1 = (irat-abs(n)) * fac0
         m = idir*abs(n)
 10      src(i,j) = src(i,j) + fac1 *
     @     (hxm1 * idir * (uf(iuf,juf+n) + uf(iuf,juf+n-jdir) +
     @                     uf(iuf+m,juf) - uf(iuf+m-idir,juf)) +
     @      hym1 * jdir * (vf(iuf+m,juf) + vf(iuf+m-idir,juf) +
     @                     vf(iuf,juf+n) - vf(iuf,juf+n-jdir)))
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgidiv(src, dims(src),
     @ uc, vc, dims(c),
     @ uf, vf, dims(f),
     @ dims(creg),
     @ hx, hy, irat, idir, jdir)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 uc(dimdec(c))
      real*8 vc(dimdec(c))
      real*8 uf(dimdec(f))
      real*8 vf(dimdec(f))
      real*8 hx, hy
      integer irat, idir, jdir
      hxm1 = 1.d0 / hx
      hym1 = 1.d0 / hy
      fac0 = dble(irat) / ((0.5d0*irat+1.d0)*irat+0.5d0)
      i = cregl0
      j = cregl1
      if (idir .eq. 1) then
         iuc = i - 1
         iuf = i * irat
      else
         iuc = i
         iuf = i * irat - 1
      endif
      if (jdir .eq. 1) then
         juc = j - 1
         juf = j * irat
      else
         juc = j
         juf = j * irat - 1
      endif
      src(i*irat,j*irat) = -fac0 *
     @  (hxm1 * idir * uc(iuc,juc) + hym1 * jdir * vc(iuc,juc))
      fac0 = fac0 / irat
      i = i * irat
      j = j * irat
      src(i,j) = src(i,j) + fac0 *
     @  (hxm1 * idir * (uf(iuf,juf) - uf(iuf-idir,juf) +
     @                  uf(iuf,juf-jdir)) +
     @   hym1 * jdir * (vf(iuf,juf) + vf(iuf-idir,juf) -
     @                  vf(iuf,juf-jdir)))
      fac0 = fac0 / irat
      do 10 n = jdir, jdir*(irat-1), jdir
         fac1 = (irat-abs(n)) * fac0
         m = idir*abs(n)
 10      src(i,j) = src(i,j) + fac1 *
     @     (hxm1 * idir * (uf(iuf,juf-n) + uf(iuf,juf-n-jdir) +
     @                     uf(iuf-m,juf) - uf(iuf-m-idir,juf)) +
     @      hym1 * jdir * (vf(iuf-m,juf) + vf(iuf-m-idir,juf) +
     @                     vf(iuf,juf-n) - vf(iuf,juf-n-jdir)))
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgddiv(src, dims(src),
     @ uc, vc, dims(c),
     @ uf, vf, dims(f),
     @ dims(creg),
     @ hx, hy, irat, jdir)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 uc(dimdec(c))
      real*8 vc(dimdec(c))
      real*8 uf(dimdec(f))
      real*8 vf(dimdec(f))
      real*8 hx, hy
      integer irat, jdir
      hxm1 = 1.d0 / hx
      hym1 = 1.d0 / hy
      fac0 = dble(irat) / ((irat+2.d0)*irat-1.d0)
      i = cregl0
      j = cregl1
      if (jdir .eq. 1) then
         juc = j
         juf = j * irat
      else
         juc = j - 1
         juf = j * irat - 1
      endif
      src(i*irat,j*irat) = fac0 *
     @  (hxm1 * (uc(i,juc-jdir) - uc(i-1,juc)) +
     @   hym1 * jdir * (vc(i-1,juc) - vc(i,juc-jdir)))
      fac0 = fac0 / irat
      i = i * irat
      j = j * irat
      src(i,j) = src(i,j) + fac0 *
     @  (hxm1 * (uf(i,juf) - uf(i-1,juf-jdir)) +
     @   hym1 * jdir * (vf(i,juf) - vf(i-1,juf-jdir)))
      fac0 = fac0 / irat
      do 10 n = jdir, jdir*(irat-1), jdir
         fac1 = (irat-abs(n)) * fac0
         m = abs(n)
 10      src(i,j) = src(i,j) + fac1 *
     @     (hxm1 * (uf(i,juf+n) + uf(i,juf+n-jdir) -
     @              uf(i-1,juf-n) - uf(i-1,juf-n-jdir) +
     @              uf(i+m,juf) - uf(i+m-1,juf) +
     @              uf(i-m,juf-jdir) - uf(i-m-1,juf-jdir)) +
     @      hym1 * jdir * (vf(i,juf+n) - vf(i,juf+n-jdir) +
     @                     vf(i-1,juf-n) - vf(i-1,juf-n-jdir) +
     @                     vf(i+m,juf) + vf(i+m-1,juf) -
     @                     vf(i-m,juf-jdir) - vf(i-m-1,juf-jdir)))
      end
