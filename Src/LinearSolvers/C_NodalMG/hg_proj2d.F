
#include "hg_version.H"

#define dims(a) a l0, a h0, a l1, a h1
#define dimdec(a) a l0:a h0,a l1:a h1
c#define dimdecu(a) (a h0-a l0+1)*(a h1-a l1+1)
#define dimdecu(a) *

      subroutine hggrad(gpx, gpy, dims(gp),
     @ dest, dims(dest),
     @ dims(freg),
     @ hx, hy, irz)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(gp)
      integer dims(dest)
      integer dims(freg)
      real*8 gpx(dimdec(gp))
      real*8 gpy(dimdec(gp))
      real*8 dest(dimdec(dest))
      real*8 hx, hy
      integer irz
      hxm1h = 0.5d0 / hx
      hym1h = 0.5d0 / hy
      do 10 j = fregl1, fregh1
         do 10 i = fregl0, fregh0
            gpx(i,j) = hxm1h * (dest(i+1,j) + dest(i+1,j+1) -
     @                          dest(i,j) - dest(i,j+1))
 10         gpy(i,j) = hym1h * (dest(i,j+1) + dest(i+1,j+1) -
     @                          dest(i,j) - dest(i+1,j))
      if (irz .eq. 1) then
         fac = hx / (12.d0 * hy)
         do 20 j = fregl1, fregh1
            r = (fregl0 - 0.5d0) * hx
            do 20 i = fregl0, fregh0
               r = r + hx
 20            gpy(i,j) = gpy(i,j) + (fac / r) *
     @                               (dest(i+1,j+1) - dest(i+1,j) -
     @                                dest(i,j+1) + dest(i,j))
      endif
      end

      subroutine hgdiv(src, dims(src),
     @ uf, vf, dims(f),
     @ dims(freg),
     @ hx, hy, irz, imax)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(src)
      integer dims(f)
      integer dims(freg)
      real*8 src(dimdec(src))
      real*8 uf(dimdec(f))
      real*8 vf(dimdec(f))
      real*8 hx, hy
      integer irz, imax
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
      if (irz .eq. 1) then
         fac = hx / (12.d0 * hy)
         r1  = (fregl0 - 0.5d0) * hx
         r1m = 1.d0 / r1
         do 20 i = fregl0, fregh0
            r0m = r1m
            r1  = r1 + hx
            if (i .lt. imax) then
               r1m = 1.d0 / r1
            else
               r1m = -r0m
            endif
            do 20 j = fregl1, fregh1
 20            src(i,j) = src(i,j) + fac *
     @            (r0m * (vf(i-1,j) - vf(i-1,j-1)) +
     @             r1m * (vf(i,j-1) - vf(i,j)))
      endif
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgfdiv(src, dims(src),
     @ uc, vc, dims(c),
     @ uf, vf, dims(f),
     @ dims(creg),
     @ hx, hy, irat, idim, idir, irz, imax)
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
      integer irat, idim, idir, irz, imax
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
         if (irz .eq. 1) then
            r = (iuc + 0.5d0) * (hx * irat)
            rfac = hx / (6.d0 * hy * r)
            do 20 j = cregl1, cregh1
 20            src(i*irat,j*irat) = src(i*irat,j*irat) + irat * fac0 *
     @           (rfac * idir * (vc(iuc,j) - vc(iuc,j-1)))
            r = (iuf + 0.5d0) * hx
            rfac = hx / (6.d0 * hy * r)
         endif
         fac0 = fac0 / (irat*irat)
         i = i * irat
         do 50 n = 0, irat-1
            fac1 = (irat-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do 30 j = irat*cregl1, irat*cregh1, irat
 30            src(i,j) = src(i,j) + fac1 *
     @           (hxm1 * idir * (uf(iuf,j-n) + uf(iuf,j-n-1) +
     @                           uf(iuf,j+n) + uf(iuf,j+n-1)) +
     @            hym1 * (vf(iuf,j-n) - vf(iuf,j-n-1) +
     @                    vf(iuf,j+n) - vf(iuf,j+n-1)))
            if (irz .eq. 1) then
               do 40 j = irat*cregl1, irat*cregh1, irat
 40               src(i,j) = src(i,j) - fac1 *
     @              (rfac * idir * (vf(iuf,j-n) - vf(iuf,j-n-1) +
     @                              vf(iuf,j+n) - vf(iuf,j+n-1)))
            endif
 50         continue
      else
         j = cregl1
         if (idir .eq. 1) then
            juc = j - 1
            juf = j * irat
         else
            juc = j
            juf = j * irat - 1
         endif
         do 60 i = cregl0, cregh0
 60         src(i*irat,j*irat) = fac0 *
     @        (hxm1 * (uc(i,juc) - uc(i-1,juc)) -
     @         hym1 * idir * (vc(i,juc) + vc(i-1,juc)))
         if (irz .eq. 1 .and. cregh0 .lt. imax) then
            do 70 i = cregl0, cregh0
               r = (i + 0.5d0) * (hx * irat)
               rfac0 = hx / (6.d0 * hy * (r - hx * irat))
               rfac1 = hx / (6.d0 * hy * r)
 70            src(i*irat,j*irat) = src(i*irat,j*irat) + irat * fac0 *
     @           (idir * (rfac1 * vc(i,juc) - rfac0 * vc(i-1,juc)))
         else if (irz .eq. 1) then
c This should only occur with a corner at the outer boundary:
            i = cregh0
            r = (i - 0.5d0) * (hx * irat)
            rfac0 = hx / (6.d0 * hy * r)
            rfac1 = -rfac0
            src(i*irat,j*irat) = src(i*irat,j*irat) + irat * fac0 *
     @           (idir * (rfac1 * vc(i,juc) - rfac0 * vc(i-1,juc)))
         endif
         fac0 = fac0 / (irat*irat)
         j = j * irat
         do 100 m = 0, irat-1
            fac1 = (irat-m) * fac0
            if (m .eq. 0) fac1 = 0.5d0 * fac1
            do 80 i = irat*cregl0, irat*cregh0, irat
 80            src(i,j) = src(i,j) + fac1 *
     @           (hxm1 * (uf(i-m,juf) - uf(i-m-1,juf) +
     @                    uf(i+m,juf) - uf(i+m-1,juf)) +
     @            hym1 * idir * (vf(i-m,juf) + vf(i-m-1,juf) +
     @                           vf(i+m,juf) + vf(i+m-1,juf)))
            if (irz .eq. 1 .and. cregh0 .lt. imax) then
               do 90 i = irat*cregl0, irat*cregh0, irat
                  r = (i + 0.5d0) * hx
                  rfac0m = hx / (6.d0 * hy * (r - (m + 1) * hx))
                  rfac1m = hx / (6.d0 * hy * (r - m * hx))
                  rfac0p = hx / (6.d0 * hy * (r + (m - 1) * hx))
                  rfac1p = hx / (6.d0 * hy * (r + m * hx))
 90               src(i,j) = src(i,j) - fac1 *
     @              (idir * (rfac1m * vf(i-m,juf) - rfac0m * vf(i-m-1,juf) +
     @                       rfac1p * vf(i+m,juf) - rfac0p * vf(i+m-1,juf)))
            else if (irz .eq. 1) then
c This should only occur with a corner at the outer boundary:
               i = irat*cregh0
               r = (i + 0.5d0) * hx
               rfac0m = hx / (6.d0 * hy * (r - (m + 1) * hx))
               if (m .eq. 0) then
                  rfac1m = -rfac0m
               else
                  rfac1m = hx / (6.d0 * hy * (r - m * hx))
               endif
               rfac0p = -rfac1m
               rfac1p = -rfac0m
               src(i,j) = src(i,j) - fac1 *
     @              (idir * (rfac1m * vf(i-m,juf) - rfac0m * vf(i-m-1,juf) +
     @                       rfac1p * vf(i+m,juf) - rfac0p * vf(i+m-1,juf)))
            endif
 100        continue
      endif
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgodiv(src, dims(src),
     @ uc, vc, dims(c),
     @ uf, vf, dims(f),
     @ dims(creg),
     @ hx, hy, irat, idir, jdir, irz)
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
      integer irat, idir, jdir, irz
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
      if (irz .eq. 1) then
         r = (iuc + 0.5d0) * (hx * irat)
         rfac0 = hx / (6.d0 * hy * (r + idir * hx * irat))
         rfac1 = hx / (6.d0 * hy * r)
         src(i*irat,j*irat) = src(i*irat,j*irat) + irat * fac0 *
     @     (idir * jdir * (rfac1 * (vc(iuc,juc+jdir) - vc(iuc,juc)) +
     @                     rfac0 *  vc(iuc+idir,juc)))
      endif
      fac0 = fac0 / irat
      i = i * irat
      j = j * irat
      src(i,j) = src(i,j) + fac0 *
     @  (hxm1 * idir * uf(iuf,juf) + hym1 * jdir * vf(iuf,juf))
      if (irz .eq. 1) then
         r = (iuf + 0.5d0) * hx
         rfac = hx / (6.d0 * hy * r)
         src(i,j) = src(i,j) - fac0 *
     @     (idir * jdir * rfac * vf(iuf,juf))
      endif
      fac0 = fac0 / irat
      do 10 n = jdir, jdir*(irat-1), jdir
         fac1 = (irat-abs(n)) * fac0
         m = idir*abs(n)
         src(i,j) = src(i,j) + fac1 *
     @     (hxm1 * idir * (uf(iuf,juf+n) + uf(iuf,juf+n-jdir) +
     @                     uf(iuf+m,juf) - uf(iuf+m-idir,juf)) +
     @      hym1 * jdir * (vf(iuf+m,juf) + vf(iuf+m-idir,juf) +
     @                     vf(iuf,juf+n) - vf(iuf,juf+n-jdir)))
         if (irz .eq. 1) then
            r = (iuf + 0.5d0) * hx
            rfac1  = hx / (6.d0 * hy * r)
            rfac0p = hx / (6.d0 * hy * (r + (m - idir) * hx))
            rfac1p = hx / (6.d0 * hy * (r + m * hx))
            src(i,j) = src(i,j) - fac1 *
     @        (idir * jdir * (rfac1p * vf(iuf+m,juf) -
     @                        rfac0p * vf(iuf+m-idir,juf) +
     @                        rfac1 * (vf(iuf,juf+n) - vf(iuf,juf+n-jdir))))
         endif
 10      continue
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgidiv(src, dims(src),
     @ uc, vc, dims(c),
     @ uf, vf, dims(f),
     @ dims(creg),
     @ hx, hy, irat, idir, jdir, irz)
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
      integer irat, idir, jdir, irz
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
      if (irz .eq. 1) then
         r = (iuc + 0.5d0) * (hx * irat)
         rfac = hx / (6.d0 * hy * r)
         src(i*irat,j*irat) = src(i*irat,j*irat) - irat * fac0 *
     @     (idir * jdir * rfac * vc(iuc,juc))
      endif
      fac0 = fac0 / irat
      i = i * irat
      j = j * irat
      src(i,j) = src(i,j) + fac0 *
     @  (hxm1 * idir * (uf(iuf,juf) - uf(iuf-idir,juf) +
     @                  uf(iuf,juf-jdir)) +
     @   hym1 * jdir * (vf(iuf,juf) + vf(iuf-idir,juf) -
     @                  vf(iuf,juf-jdir)))
      if (irz .eq. 1) then
         r = (iuf + 0.5d0) * hx
         rfac0 = hx / (6.d0 * hy * (r - idir * hx))
         rfac1 = hx / (6.d0 * hy * r)
         src(i,j) = src(i,j) + fac0 *
     @     (idir * jdir * (rfac1 * (vf(iuf,juf-jdir) - vf(iuf,juf)) +
     @                     rfac0 *  vf(iuf-idir,juf)))
      endif
      fac0 = fac0 / irat
      do 10 n = jdir, jdir*(irat-1), jdir
         fac1 = (irat-abs(n)) * fac0
         m = idir*abs(n)
         src(i,j) = src(i,j) + fac1 *
     @     (hxm1 * idir * (uf(iuf,juf-n) + uf(iuf,juf-n-jdir) +
     @                     uf(iuf-m,juf) - uf(iuf-m-idir,juf)) +
     @      hym1 * jdir * (vf(iuf-m,juf) + vf(iuf-m-idir,juf) +
     @                     vf(iuf,juf-n) - vf(iuf,juf-n-jdir)))
         if (irz .eq. 1) then
            r = (iuf + 0.5d0) * hx
            rfac0m = hx / (6.d0 * hy * (r - (m + idir) * hx))
            rfac1m = hx / (6.d0 * hy * (r - m * hx))
            rfac1  = hx / (6.d0 * hy * r)
            src(i,j) = src(i,j) - fac1 *
     @        (idir * jdir * (rfac1m * vf(iuf-m,juf) -
     @                        rfac0m * vf(iuf-m-idir,juf) +
     @                        rfac1 * (vf(iuf,juf-n) - vf(iuf,juf-n-jdir))))
         endif
 10      continue
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgddiv(src, dims(src),
     @ uc, vc, dims(c),
     @ uf, vf, dims(f),
     @ dims(creg),
     @ hx, hy, irat, jdir, irz)
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
      integer irat, jdir, irz
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
      if (irz .eq. 1) then
         r = (i + 0.5d0) * (hx * irat)
         rfac0 = hx / (6.d0 * hy * (r - hx * irat))
         rfac1 = hx / (6.d0 * hy * r)
         src(i*irat,j*irat) = src(i*irat,j*irat) + irat * fac0 *
     @     (jdir * (rfac1 * vc(i,juc-jdir) +
     @              rfac0 * vc(i-1,juc)))
      endif
      fac0 = fac0 / irat
      i = i * irat
      j = j * irat
      src(i,j) = src(i,j) + fac0 *
     @  (hxm1 * (uf(i,juf) - uf(i-1,juf-jdir)) +
     @   hym1 * jdir * (vf(i,juf) - vf(i-1,juf-jdir)))
      if (irz .eq. 1) then
         r = (i + 0.5d0) * hx
         rfac0 = hx / (6.d0 * hy * (r - hx))
         rfac1 = hx / (6.d0 * hy * r)
         src(i,j) = src(i,j) - fac0 *
     @     (jdir * (rfac1 * vf(i,juf) +
     @              rfac0 * vf(i-1,juf-jdir)))
      endif
      fac0 = fac0 / irat
      do 10 n = jdir, jdir*(irat-1), jdir
         fac1 = (irat-abs(n)) * fac0
         m = abs(n)
         src(i,j) = src(i,j) + fac1 *
     @     (hxm1 * (uf(i,juf+n) + uf(i,juf+n-jdir) -
     @              uf(i-1,juf-n) - uf(i-1,juf-n-jdir) +
     @              uf(i+m,juf) - uf(i+m-1,juf) +
     @              uf(i-m,juf-jdir) - uf(i-m-1,juf-jdir)) +
     @      hym1 * jdir * (vf(i,juf+n) - vf(i,juf+n-jdir) +
     @                     vf(i-1,juf-n) - vf(i-1,juf-n-jdir) +
     @                     vf(i+m,juf) + vf(i+m-1,juf) -
     @                     vf(i-m,juf-jdir) - vf(i-m-1,juf-jdir)))
         if (irz .eq. 1) then
            r = (i + 0.5d0) * hx
            rfac0m = hx / (6.d0 * hy * (r - (m + 1) * hx))
            rfac1m = hx / (6.d0 * hy * (r - m * hx))
            rfac0  = hx / (6.d0 * hy * (r - hx))
            rfac1  = hx / (6.d0 * hy * r)
            rfac0p = hx / (6.d0 * hy * (r + (m - 1) * hx))
            rfac1p = hx / (6.d0 * hy * (r + m * hx))
            src(i,j) = src(i,j) - fac1 *
     @        (jdir * (rfac1p * vf(i+m,juf) -
     @                 rfac0p * vf(i+m-1,juf) +
     @                 rfac1 * (vf(i,juf+n) - vf(i,juf+n-jdir)) -
     @                 rfac0 * (vf(i-1,juf-n) - vf(i-1,juf-n-jdir)) -
     @                 rfac1m * vf(i-m,juf-jdir) +
     @                 rfac0m * vf(i-m-1,juf-jdir)))
         endif
 10      continue
      end
