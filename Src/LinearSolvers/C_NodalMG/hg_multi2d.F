
#include "hg_version.H"

#define dims(a) a l0, a h0, a l1, a h1
#define dimdec(a) a l0:a h0,a l1:a h1
c#define dimdecu(a) (a h0-a l0+1)*(a h1-a l1+1)
#define dimdecu(a) *

#ifdef CROSS_STENCIL

#ifndef CONSTANT

c variable density versions:

c Note---assumes fdst linearly interpolated from cdst along edge
      subroutine hgfres(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ sigmaf, dims(sf),
     @ sigmac, dims(sc),
     @ dims(reg),
     @ hx, hy, ir, jr, idim, idir)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(sf)
      integer dims(sc)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 sigmaf(dimdec(sf))
      real*8 sigmac(dimdec(sc))
      real*8 hx, hy
      integer ir, jr, idim, idir
      if (idim .eq. 0) then
         i = regl0
         if (idir .eq. 1) then
            is = i - 1
         else
            is = i
         endif
         fac0 = ir / (ir + 1.d0)
         hxm2 = 1.d0 / (ir * ir * hx * hx)
         hym2 = 1.d0 / (jr * jr * hy * hy)
         do 10 j = regl1, regh1
 10         res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     @        (hxm2 *
     @          ((sigmac(is,j-1) + sigmac(is,j)) *
     @            (cdst(i-idir,j) - cdst(i,j))) +
     @         hym2 *
     @          (sigmac(is,j-1) *
     @            (cdst(i,j-1) - cdst(i,j)) +
     @           sigmac(is,j) *
     @            (cdst(i,j+1) - cdst(i,j))))
         fac0 = fac0 / (ir * jr * jr)
         hxm2 = ir * ir * hxm2
         hym2 = jr * jr * hym2
         i = i * ir
         if (idir .eq. 1) then
            is = i
         else
            is = i - 1
         endif
         do 20 n = 0, jr-1
            fac1 = (jr-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do 20 j = jr*regl1, jr*regh1, jr
               tmp = hxm2 *
     @           ((sigmaf(is,j-n-1) + sigmaf(is,j-n)) *
     @             (fdst(i+idir,j-n) - fdst(i,j-n)) +
     @		  (sigmaf(is,j+n-1) + sigmaf(is,j+n)) *
     @             (fdst(i+idir,j+n) - fdst(i,j+n)))
 20            res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     @		 (sigmaf(is,j-n-1) *
     @             (fdst(i,j-n-1) - fdst(i,j-n)) +
     @		  sigmaf(is,j-n) *
     @             (fdst(i,j-n+1) - fdst(i,j-n)) +
     @		  sigmaf(is,j+n-1) *
     @             (fdst(i,j+n-1) - fdst(i,j+n)) +
     @		  sigmaf(is,j+n) *
     @             (fdst(i,j+n+1) - fdst(i,j+n))))
      else
         j = regl1
         if (idir .eq. 1) then
            js = j - 1
         else
            js = j
         endif
         fac0 = jr / (jr + 1.d0)
         hxm2 = 1.d0 / (ir * ir * hx * hx)
         hym2 = 1.d0 / (jr * jr * hy * hy)
         do 30 i = regl0, regh0
 30         res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     @        (hxm2 *
     @          (sigmac(i-1,js) *
     @            (cdst(i-1,j) - cdst(i,j)) +
     @           sigmac(i,js) *
     @            (cdst(i+1,j) - cdst(i,j))) +
     @         hym2 *
     @          ((sigmac(i-1,js) + sigmac(i,js)) *
     @            (cdst(i,j-idir) - cdst(i,j))))
         fac0 = fac0 / (ir * ir * jr)
         hxm2 = ir * ir * hxm2
         hym2 = jr * jr * hym2
         j = j * jr
         if (idir .eq. 1) then
            js = j
         else
            js = j - 1
         endif
         do 40 m = 0, ir-1
            fac1 = (ir-m) * fac0
            if (m .eq. 0) fac1 = 0.5d0 * fac1
            do 40 i = ir*regl0, ir*regh0, ir
               tmp = hxm2 *
     @		 (sigmaf(i-m-1,js) *
     @             (fdst(i-m-1,j) - fdst(i-m,j)) +
     @		  sigmaf(i-m,js) *
     @             (fdst(i-m+1,j) - fdst(i-m,j)) +
     @		  sigmaf(i+m-1,js) *
     @             (fdst(i+m-1,j) - fdst(i+m,j)) +
     @		  sigmaf(i+m,js) *
     @             (fdst(i+m+1,j) - fdst(i+m,j)))
 40            res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     @           ((sigmaf(i-m-1,js) + sigmaf(i-m,js)) *
     @             (fdst(i-m,j+idir) - fdst(i-m,j)) +
     @		  (sigmaf(i+m-1,js) + sigmaf(i+m,js)) *
     @             (fdst(i+m,j+idir) - fdst(i+m,j))))
      endif
      end

c Note---assumes fdst linearly interpolated from cdst along edges
      subroutine hgcres(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ sigmaf, dims(sf),
     @ sigmac, dims(sc),
     @ dims(reg),
     @ hx, hy, ir, jr, ga)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(sf)
      integer dims(sc)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 sigmaf(dimdec(sf))
      real*8 sigmac(dimdec(sc))
      real*8 hx, hy
      integer ir, jr, ga(0:1,0:1)
      hxm2c = 1.d0 / (ir * ir * hx * hx)
      hym2c = 1.d0 / (jr * jr * hy * hy)
      hxm2 = ir * ir * hxm2c
      hym2 = jr * jr * hym2c
      ic = regl0
      jc = regl1
      if = ic * ir
      jf = jc * jr
      sum = 0.d0
      center = 0.d0
c quadrants
      ffac = 0.5d0
      cfac = 0.5d0 * ir * jr
      do 10 ji = 0, 1
         jdir = 2 * ji - 1
         do 10 ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,ji) .eq. 1) then
               center = center + ffac
               sum = sum + sigmaf(if+ii-1,jf+ji-1) *
     @           (hxm2 * (fdst(if+idir,jf) - fdst(if,jf)) +
     @            hym2 * (fdst(if,jf+jdir) - fdst(if,jf)))
            else
               center = center + cfac
               sum = sum + ir * jr * sigmac(ic+ii-1,jc+ji-1) *
     @           (hxm2c * (cdst(ic+idir,jc) - cdst(ic,jc)) +
     @            hym2c * (cdst(ic,jc+jdir) - cdst(ic,jc)))
            endif
 10         continue
c edges
      do 40 ji = 0, 1
         jdir = 2 * ji - 1
         do 40 ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,ji) - ga(ii,1-ji) .eq. 1) then
               fac1 = 1.d0 / ir
               ffac = 0.5d0 * (ir-1)
               center = center + ffac
               do 20 m = idir, idir*(ir-1), idir
                  fac = (ir-abs(m)) * fac1
 20               sum = sum + fac *
     @              (hxm2 * (sigmaf(if+m-1,jf+ji-1) *
     @                        (fdst(if+m-1,jf) - fdst(if+m,jf)) +
     @                       sigmaf(if+m,jf+ji-1) *
     @                        (fdst(if+m+1,jf) - fdst(if+m,jf))) +
     @               hym2 *
     @                 (sigmaf(if+m-1,jf+ji-1) + sigmaf(if+m,jf+ji-1)) *
     @                 (fdst(if+m,jf+jdir) - fdst(if+m,jf)))
            endif
            if (ga(ii,ji) - ga(1-ii,ji) .eq. 1) then
               fac1 = 1.d0 / jr
               ffac = 0.5d0 * (jr-1)
               center = center + ffac
               do 30 n = jdir, jdir*(jr-1), jdir
                  fac = (jr-abs(n)) * fac1
 30               sum = sum + fac *
     @              (hxm2 *
     @                 (sigmaf(if+ii-1,jf+n-1) + sigmaf(if+ii-1,jf+n)) *
     @                 (fdst(if+idir,jf+n) - fdst(if,jf+n)) +
     @               hym2 * (sigmaf(if+ii-1,jf+n-1) *
     @                        (fdst(if,jf+n-1) - fdst(if,jf+n)) +
     @                       sigmaf(if+ii-1,jf+n) *
     @                        (fdst(if,jf+n+1) - fdst(if,jf+n))))
            endif
 40         continue
c weighting
      res(if,jf) = src(if,jf) - sum / center
      end

#else

c constant density versions:

c Note---assumes fdst linearly interpolated from cdst along edge
      subroutine hgfres(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ dims(reg),
     @ hx, irat, idim, idir)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 hx
      integer irat, idim, idir
      ir = irat
      jr = irat
      hm2 = 1.d0 / (hx*hx)
      fac0 = 1.d0 / (irat*irat+irat)
      if (idim .eq. 0) then
         i = regl0
         do 10 j = regl1, regh1
 10         res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 * hm2 *
     @        (2.d0 * cdst(i-idir,j) +
     @         cdst(i,j-1) + cdst(i,j+1) -
     @         4.d0 * cdst(i,j))
         fac0 = fac0 / irat
         i = i * irat
         do 20 n = 0, irat-1
            fac1 = (irat-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do 20 j = irat*regl1, irat*regh1, irat
 20            res(i,j) = res(i,j) - fac1 * hm2 *
     @           (2.d0 *
     @            (fdst(i+idir,j-n) + fdst(i+idir,j+n)) +
     @            (fdst(i,j-n-1) + fdst(i,j-n+1) +
     @             fdst(i,j+n-1) + fdst(i,j+n+1)) -
     @            4.d0 *
     @            (fdst(i,j-n) + fdst(i,j+n)))
      else
         j = regl1
         do 30 i = regl0, regh0
 30         res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 * hm2 *
     @        (2.d0 * cdst(i,j-idir) +
     @         cdst(i-1,j) + cdst(i+1,j) -
     @         4.d0 * cdst(i,j))
         fac0 = fac0 / irat
         j = j * irat
         do 40 m = 0, irat-1
            fac1 = (irat-m) * fac0
            if (m .eq. 0) fac1 = 0.5d0 * fac1
            do 40 i = irat*regl0, irat*regh0, irat
 40            res(i,j) = res(i,j) - fac1 * hm2 *
     @           ((fdst(i-m-1,j) + fdst(i-m+1,j) +
     @             fdst(i+m-1,j) + fdst(i+m+1,j)) +
     @            2.d0 *
     @            (fdst(i-m,j+idir) + fdst(i+m,j+idir)) -
     @            4.d0 *
     @            (fdst(i-m,j) + fdst(i+m,j)))
      endif
      end

c Note---assumes fdst linearly interpolated from cdst along edges
      subroutine hgcres(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ dims(reg),
     @ hx, irat, ga)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 hx, hy
      integer irat, ga(0:1,0:1)
      hm2 = 1.d0 / (hx*hx)
      ic = regl0
      jc = regl1
      if = ic * irat
      jf = jc * irat
      sum = 0.d0
      center = 0.d0
c quadrants
      ffac = 0.5d0
      cfac = 0.5d0 * irat * irat
      do 10 ji = 0, 1
         jdir = 2 * ji - 1
         do 10 ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,ji) .eq. 1) then
               center = center + ffac
               sum = sum + hm2 *
     @           (fdst(if+idir,jf) + fdst(if,jf+jdir) -
     @            2.d0 * fdst(if,jf))
            else
               center = center + cfac
               sum = sum + hm2 *
     @           (cdst(ic+idir,jc) + cdst(ic,jc+jdir) -
     @            2.d0 * cdst(ic,jc))
            endif
 10         continue
c edges
      fac1 = 1.d0 / irat
      ffac = 0.5d0 * (irat-1)
      do 40 ji = 0, 1
         jdir = 2 * ji - 1
         do 40 ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,ji) - ga(ii,1-ji) .eq. 1) then
               center = center + ffac
               do 20 m = idir, idir*(irat-1), idir
                  fac = (irat-abs(m)) * fac1
 20               sum = sum + fac * hm2 *
     @              ((fdst(if+m-1,jf) + fdst(if+m+1,jf)) +
     @               2.d0 * fdst(if+m,jf+jdir) -
     @               4.d0 * fdst(if+m,jf))
            endif
            if (ga(ii,ji) - ga(1-ii,ji) .eq. 1) then
               center = center + ffac
               do 30 n = jdir, jdir*(irat-1), jdir
                  fac = (irat-abs(n)) * fac1
 30               sum = sum + fac * hm2 *
     @              (2.d0 * fdst(if+idir,jf+n) +
     @               (fdst(if,jf+n-1) + fdst(if,jf+n+1)) -
     @               4.d0 * fdst(if,jf+n))
            endif
 40         continue
c weighting
      res(if,jf) = src(if,jf) - sum / center
      end

#endif

#else

c 9-point stencil versions:

#ifndef CONSTANT

c#define TURNOFF

c variable density versions:

c Note---assumes fdst linearly interpolated from cdst along edge
      subroutine hgfres(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ sigmaf, dims(sf),
     @ sigmac, dims(sc),
     @ dims(reg),
#ifdef TURNOFF
     @ hx, hy, ir, jr, idim, idir, irz1, imax)
#else
     @ hx, hy, ir, jr, idim, idir, irz, imax)
#endif
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(sf)
      integer dims(sc)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 sigmaf(dimdec(sf))
      real*8 sigmac(dimdec(sc))
      real*8 hx, hy
      integer ir, jr, idim, idir, irz, imax
#ifdef TURNOFF
      irz = 0
#endif
      if (idim .eq. 0) then
         i = regl0
         if (idir .eq. 1) then
            is = i - 1
         else
            is = i
         endif
         fac0 = ir / (3.d0 * (ir + 1.d0))
         hxm2 = 1.d0 / (ir * ir * hx * hx)
         hym2 = 1.d0 / (jr * jr * hy * hy)
         do 10 j = regl1, regh1
 10         res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     @        (hxm2 *
     @          (sigmac(is,j-1) *
     @            (cdst(i-idir,j-1) - cdst(i,j-1) +
     @             2.d0 * (cdst(i-idir,j) - cdst(i,j))) +
     @           sigmac(is,j) *
     @            (cdst(i-idir,j+1) - cdst(i,j+1) +
     @             2.d0 * (cdst(i-idir,j) - cdst(i,j)))) +
     @         hym2 *
     @          (sigmac(is,j-1) *
     @            (cdst(i-idir,j-1) - cdst(i-idir,j) +
     @             2.d0 * (cdst(i,j-1) - cdst(i,j))) +
     @           sigmac(is,j) *
     @            (cdst(i-idir,j+1) - cdst(i-idir,j) +
     @             2.d0 * (cdst(i,j+1) - cdst(i,j)))))
         if (irz .eq. 1) then
            r = (is + 0.5d0) * (hx * ir)
            rfac = idir * ir * hx * hym2 / (2.d0 * r)
            do 20 j = regl1, regh1
 20            res(i*ir,j*jr) = res(i*ir,j*jr) - fac0 *
     @           (rfac * (sigmac(is,j-1) * (cdst(i,j-1) - cdst(i,j)) +
     @                    sigmac(is,j)   * (cdst(i,j+1) - cdst(i,j))))
         endif
         fac0 = fac0 / (ir * jr * jr)
         hxm2 = ir * ir * hxm2
         hym2 = jr * jr * hym2
         i = i * ir
         if (idir .eq. 1) then
            is = i
         else
            is = i - 1
         endif
         if (irz .eq. 1) then
            r = (is + 0.5d0) * hx
            rfac = idir * hx * hym2 / (2.d0 * r)
         endif
         do 50 n = 0, jr-1
            fac1 = (jr-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do 30 j = jr*regl1, jr*regh1, jr
               tmp = hxm2 *
     @           (sigmaf(is,j-n-1) *
     @             (fdst(i+idir,j-n-1) - fdst(i,j-n-1) +
     @              2.d0 * (fdst(i+idir,j-n) - fdst(i,j-n))) +
     @		  sigmaf(is,j-n) *
     @             (fdst(i+idir,j-n+1) - fdst(i,j-n+1) +
     @              2.d0 * (fdst(i+idir,j-n) - fdst(i,j-n))) +
     @		  sigmaf(is,j+n-1) *
     @             (fdst(i+idir,j+n-1) - fdst(i,j+n-1) +
     @              2.d0 * (fdst(i+idir,j+n) - fdst(i,j+n))) +
     @		  sigmaf(is,j+n) *
     @             (fdst(i+idir,j+n+1) - fdst(i,j+n+1) +
     @              2.d0 * (fdst(i+idir,j+n) - fdst(i,j+n))))
 30            res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     @		 (sigmaf(is,j-n-1) *
     @             (fdst(i+idir,j-n-1) - fdst(i+idir,j-n) +
     @              2.d0 * (fdst(i,j-n-1) - fdst(i,j-n))) +
     @		  sigmaf(is,j-n) *
     @             (fdst(i+idir,j-n+1) - fdst(i+idir,j-n) +
     @              2.d0 * (fdst(i,j-n+1) - fdst(i,j-n))) +
     @		  sigmaf(is,j+n-1) *
     @             (fdst(i+idir,j+n-1) - fdst(i+idir,j+n) +
     @              2.d0 * (fdst(i,j+n-1) - fdst(i,j+n))) +
     @		  sigmaf(is,j+n) *
     @             (fdst(i+idir,j+n+1) - fdst(i+idir,j+n) +
     @              2.d0 * (fdst(i,j+n+1) - fdst(i,j+n)))))
            if (irz .eq. 1) then
               do 40 j = jr*regl1, jr*regh1, jr
 40               res(i,j) = res(i,j) + fac1 *
     @           (rfac * (sigmaf(is,j-n-1) * (fdst(i,j-n-1) - fdst(i,j-n)) +
     @                    sigmaf(is,j-n)   * (fdst(i,j-n+1) - fdst(i,j-n)) +
     @                    sigmaf(is,j+n-1) * (fdst(i,j+n-1) - fdst(i,j+n)) +
     @                    sigmaf(is,j+n)   * (fdst(i,j+n+1) - fdst(i,j+n))))
            endif
 50         continue
      else
         j = regl1
         if (idir .eq. 1) then
            js = j - 1
         else
            js = j
         endif
         fac0 = jr / (3.d0 * (jr + 1.d0))
         hxm2 = 1.d0 / (ir * ir * hx * hx)
         hym2 = 1.d0 / (jr * jr * hy * hy)
         do 60 i = regl0, regh0
 60         res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     @        (hxm2 *
     @          (sigmac(i-1,js) *
     @            (cdst(i-1,j-idir) - cdst(i,j-idir) +
     @             2.d0 * (cdst(i-1,j) - cdst(i,j))) +
     @           sigmac(i,js) *
     @            (cdst(i+1,j-idir) - cdst(i,j-idir) +
     @             2.d0 * (cdst(i+1,j) - cdst(i,j)))) +
     @         hym2 *
     @          (sigmac(i-1,js) *
     @            (cdst(i-1,j-idir) - cdst(i-1,j) +
     @             2.d0 * (cdst(i,j-idir) - cdst(i,j))) +
     @           sigmac(i,js) *
     @            (cdst(i+1,j-idir) - cdst(i+1,j) +
     @             2.d0 * (cdst(i,j-idir) - cdst(i,j)))))
         if (irz .eq. 1 .and. regh0 .lt. imax) then
            do 70 i = regl0, regh0
               r = (i + 0.5d0) * (hx * ir)
               rfac0 = ir * hx * hym2 / (2.d0 * (r - hx * ir))
               rfac1 = ir * hx * hym2 / (2.d0 * r)
 70            res(i*ir,j*jr) = res(i*ir,j*jr) - fac0 *
     @           (rfac0 * sigmac(i-1,js) - rfac1 * sigmac(i,js)) *
     @           (cdst(i,j-idir) - cdst(i,j))
         else if (irz .eq. 1) then
c This should only occur with a corner at the outer boundary:
            i = regh0
            r = (i - 0.5d0) * (hx * ir)
            rfac0 = ir * hx * hym2 / (2.d0 * r)
            rfac1 = -rfac0
            res(i*ir,j*jr) = res(i*ir,j*jr) - fac0 *
     @           (rfac0 * sigmac(i-1,js) - rfac1 * sigmac(i,js)) *
     @           (cdst(i,j-idir) - cdst(i,j))
         endif
         fac0 = fac0 / (ir * ir * jr)
         hxm2 = ir * ir * hxm2
         hym2 = jr * jr * hym2
         j = j * jr
         if (idir .eq. 1) then
            js = j
         else
            js = j - 1
         endif
         do 100 m = 0, ir-1
            fac1 = (ir-m) * fac0
            if (m .eq. 0) fac1 = 0.5d0 * fac1
            do 80 i = ir*regl0, ir*regh0, ir
               tmp = hxm2 *
     @		 (sigmaf(i-m-1,js) *
     @             (fdst(i-m-1,j+idir) - fdst(i-m,j+idir) +
     @              2.d0 * (fdst(i-m-1,j) - fdst(i-m,j))) +
     @		  sigmaf(i-m,js) *
     @             (fdst(i-m+1,j+idir) - fdst(i-m,j+idir) +
     @              2.d0 * (fdst(i-m+1,j) - fdst(i-m,j))) +
     @		  sigmaf(i+m-1,js) *
     @             (fdst(i+m-1,j+idir) - fdst(i+m,j+idir) +
     @              2.d0 * (fdst(i+m-1,j) - fdst(i+m,j))) +
     @		  sigmaf(i+m,js) *
     @             (fdst(i+m+1,j+idir) - fdst(i+m,j+idir) +
     @              2.d0 * (fdst(i+m+1,j) - fdst(i+m,j))))
 80            res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     @           (sigmaf(i-m-1,js) *
     @             (fdst(i-m-1,j+idir) - fdst(i-m-1,j) +
     @              2.d0 * (fdst(i-m,j+idir) - fdst(i-m,j))) +
     @		  sigmaf(i-m,js) *
     @             (fdst(i-m+1,j+idir) - fdst(i-m+1,j) +
     @              2.d0 * (fdst(i-m,j+idir) - fdst(i-m,j))) +
     @		  sigmaf(i+m-1,js) *
     @             (fdst(i+m-1,j+idir) - fdst(i+m-1,j) +
     @              2.d0 * (fdst(i+m,j+idir) - fdst(i+m,j))) +
     @		  sigmaf(i+m,js) *
     @             (fdst(i+m+1,j+idir) - fdst(i+m+1,j) +
     @              2.d0 * (fdst(i+m,j+idir) - fdst(i+m,j)))))
         if (irz .eq. 1 .and. regh0 .lt. imax) then
            do 90 i = ir*regl0, ir*regh0, ir
               r = (i + 0.5d0) * hx
               rfac0m = hx * hym2 / (2.d0 * (r - (m + 1) * hx))
               rfac1m = hx * hym2 / (2.d0 * (r - m * hx))
               rfac0p = hx * hym2 / (2.d0 * (r + (m - 1) * hx))
               rfac1p = hx * hym2 / (2.d0 * (r + m * hx))
 90            res(i,j) = res(i,j) - fac1 *
     @           ((rfac0m * sigmaf(i-m-1,js) - rfac1m * sigmaf(i-m,js)) *
     @            (fdst(i-m,j+idir) - fdst(i-m,j)) +
     @            (rfac0p * sigmaf(i+m-1,js) - rfac1p * sigmaf(i+m,js)) *
     @            (fdst(i+m,j+idir) - fdst(i+m,j)))
         else if (irz .eq. 1) then
c This should only occur with a corner at the outer boundary:
            i = ir * regh0
            r = (i + 0.5d0) * hx
            rfac0m = hx * hym2 / (2.d0 * (r - (m + 1) * hx))
            if (m .eq. 0) then
               rfac1m = -rfac0m
            else
               rfac1m = hx * hym2 / (2.d0 * (r - m * hx))
            endif
            rfac0p = -rfac1m
            rfac1p = -rfac0m
            res(i,j) = res(i,j) - fac1 *
     @           ((rfac0m * sigmaf(i-m-1,js) - rfac1m * sigmaf(i-m,js)) *
     @            (fdst(i-m,j+idir) - fdst(i-m,j)) +
     @            (rfac0p * sigmaf(i+m-1,js) - rfac1p * sigmaf(i+m,js)) *
     @            (fdst(i+m,j+idir) - fdst(i+m,j)))
         endif
 100     continue
      endif
      end

c Note---assumes fdst linearly interpolated from cdst along edge
      subroutine hgores(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ sigmaf, dims(sf),
     @ sigmac, dims(sc),
     @ dims(reg),
#ifdef TURNOFF
     @ hx, hy, ir, jr, idir, jdir, irz1)
#else
     @ hx, hy, ir, jr, idir, jdir, irz)
#endif
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(sf)
      integer dims(sc)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 sigmaf(dimdec(sf))
      real*8 sigmac(dimdec(sc))
      real*8 hx, hy
      integer ir, jr, idir, jdir, irz
#ifdef TURNOFF
      irz = 0
#endif
      i = regl0
      j = regl1
      if (idir .eq. 1) then
         is = i - 1
      else
         is = i
      endif
      if (jdir .eq. 1) then
         js = j - 1
      else
         js = j
      endif
      hxm2 = 1.d0 / (ir * ir * hx * hx)
      hym2 = 1.d0 / (jr * jr * hy * hy)
      fac0 = (ir * jr) / (4.5d0*ir*jr + 1.5d0*ir + 1.5d0*jr - 1.5d0)
      res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     @  (hxm2 * (sigmac(is,js+jdir) *
     @            (cdst(i-idir,j+jdir) - cdst(i,j+jdir) +
     @             2.d0 * (cdst(i-idir,j) - cdst(i,j))) +
     @           sigmac(is,js) *
     @            (cdst(i-idir,j-jdir) - cdst(i,j-jdir) +
     @             2.d0 * (cdst(i-idir,j) - cdst(i,j))) +
     @           sigmac(is+idir,js) *
     @            (cdst(i+idir,j-jdir) - cdst(i,j-jdir) +
     @             2.d0 * (cdst(i+idir,j) - cdst(i,j)))) +
     @   hym2 * (sigmac(is,js+jdir) *
     @            (cdst(i-idir,j+jdir) - cdst(i-idir,j) +
     @             2.d0 * (cdst(i,j+jdir) - cdst(i,j))) +
     @           sigmac(is,js) *
     @            (cdst(i-idir,j-jdir) - cdst(i-idir,j) +
     @             2.d0 * (cdst(i,j-jdir) - cdst(i,j))) +
     @           sigmac(is+idir,js) *
     @            (cdst(i+idir,j-jdir) - cdst(i+idir,j) +
     @             2.d0 * (cdst(i,j-jdir) - cdst(i,j)))))
      if (irz .eq. 1) then
         r = (is + 0.5d0) * (hx * ir)
         rfac0 = ir * hx * hym2 / (2.d0 * (r + idir * hx * ir))
         rfac1 = ir * hx * hym2 / (2.d0 * r)
         res(i*ir,j*jr) = res(i*ir,j*jr) + idir * fac0 *
     @     ((rfac0 * sigmac(is+idir,js) - rfac1 * sigmac(is,js)) *
     @      (cdst(i,j-jdir) - cdst(i,j)) -
     @      rfac1 * sigmac(is,js+jdir) * (cdst(i,j+jdir) - cdst(i,j)))
      endif
      fac0 = fac0 / (ir * jr)
      hxm2 = ir * ir * hxm2
      hym2 = jr * jr * hym2
      i = i * ir
      j = j * jr
      if (idir .eq. 1) then
         is = i
      else
         is = i - 1
      endif
      if (jdir .eq. 1) then
         js = j
      else
         js = j - 1
      endif
      res(i,j) = res(i,j) - fac0 * sigmaf(is,js) *
     @  (hxm2 * (fdst(i+idir,j+jdir) - fdst(i,j+jdir) +
     @           2.d0 * (fdst(i+idir,j) - fdst(i,j))) +
     @   hym2 * (fdst(i+idir,j+jdir) - fdst(i+idir,j) +
     @           2.d0 * (fdst(i,j+jdir) - fdst(i,j))))
      if (irz .eq. 1) then
         r = (is + 0.5d0) * hx
         rfac = hx * hym2 / (2.d0 * r)
         res(i,j) = res(i,j) + idir * fac0 *
     @      rfac * sigmaf(is,js) * (fdst(i,j+jdir) - fdst(i,j))
      endif
      fac0 = fac0 / ir
      do 10 m = idir, idir*(ir-1), idir
         fac1 = (ir-abs(m)) * fac0
         tmp = hxm2 *
     @     (sigmaf(i+m-1,js) *
     @       (fdst(i+m-1,j+jdir) - fdst(i+m,j+jdir) +
     @        2.d0 * (fdst(i+m-1,j) - fdst(i+m,j))) +
     @      sigmaf(i+m,js) *
     @       (fdst(i+m+1,j+jdir) - fdst(i+m,j+jdir) +
     @        2.d0 * (fdst(i+m+1,j) - fdst(i+m,j))))
         res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     @     (sigmaf(i+m-1,js) *
     @       (fdst(i+m-1,j+jdir) - fdst(i+m-1,j) +
     @        2.d0 * (fdst(i+m,j+jdir) - fdst(i+m,j))) +
     @	    sigmaf(i+m,js) *
     @       (fdst(i+m+1,j+jdir) - fdst(i+m+1,j) +
     @        2.d0 * (fdst(i+m,j+jdir) - fdst(i+m,j)))))
         if (irz .eq. 1) then
            r = (i + m + 0.5d0) * hx
            rfac0p = hx * hym2 / (2.d0 * (r - hx))
            rfac1p = hx * hym2 / (2.d0 * r)
            res(i,j) = res(i,j) - fac1 *
     @         (rfac0p * sigmaf(i+m-1,js) - rfac1p * sigmaf(i+m,js)) *
     @          (fdst(i+m,j+jdir) - fdst(i+m,j))
         endif
 10      continue
      fac0 = ir * fac0 / jr
      do 20 n = jdir, jdir*(jr-1), jdir
         fac1 = (jr-abs(n)) * fac0
         tmp = hxm2 *
     @     (sigmaf(is,j+n-1) *
     @       (fdst(i+idir,j+n-1) - fdst(i,j+n-1) +
     @        2.d0 * (fdst(i+idir,j+n) - fdst(i,j+n))) +
     @      sigmaf(is,j+n) *
     @       (fdst(i+idir,j+n+1) - fdst(i,j+n+1) +
     @        2.d0 * (fdst(i+idir,j+n) - fdst(i,j+n))))
         res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     @     (sigmaf(is,j+n-1) *
     @       (fdst(i+idir,j+n-1) - fdst(i+idir,j+n) + 
     @        2.d0 * (fdst(i,j+n-1) - fdst(i,j+n))) +
     @	    sigmaf(is,j+n) *
     @       (fdst(i+idir,j+n+1) - fdst(i+idir,j+n) +
     @        2.d0 * (fdst(i,j+n+1) - fdst(i,j+n)))))
         if (irz .eq. 1) then
            r = (is + 0.5d0) * hx
            rfac1 = -idir * hx * hym2 / (2.d0 * r)
            res(i,j) = res(i,j) - fac1 *
     @        rfac1 * (sigmaf(is,j+n-1) * (fdst(i,j+n-1) - fdst(i,j+n)) +
     @                 sigmaf(is,j+n)   * (fdst(i,j+n+1) - fdst(i,j+n)))
         endif
 20      continue
      end

c Note---assumes fdst linearly interpolated from cdst along edge
      subroutine hgires(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ sigmaf, dims(sf),
     @ sigmac, dims(sc),
     @ dims(reg),
#ifdef TURNOFF
     @ hx, hy, ir, jr, idir, jdir, irz1)
#else
     @ hx, hy, ir, jr, idir, jdir, irz)
#endif
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(sf)
      integer dims(sc)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 sigmaf(dimdec(sf))
      real*8 sigmac(dimdec(sc))
      real*8 hx, hy
      integer ir, jr, idir, jdir, irz
#ifdef TURNOFF
      irz = 0
#endif
      i = regl0
      j = regl1
      if (idir .eq. 1) then
         is = i - 1
      else
         is = i
      endif
      if (jdir .eq. 1) then
         js = j - 1
      else
         js = j
      endif
      hxm2 = 1.d0 / (ir * ir * hx * hx)
      hym2 = 1.d0 / (jr * jr * hy * hy)
      fac0 = (ir * jr) / (1.5d0*ir*jr + 1.5d0*ir + 1.5d0*jr + 1.5d0)
      res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 * sigmac(is,js) *
     @  (hxm2 * (cdst(i-idir,j-jdir) - cdst(i,j-jdir) +
     @           2.d0 * (cdst(i-idir,j) - cdst(i,j))) +
     @   hym2 * (cdst(i-idir,j-jdir) - cdst(i-idir,j) +
     @           2.d0 * (cdst(i,j-jdir) - cdst(i,j))))
      if (irz .eq. 1) then
         r = (is + 0.5d0) * (hx * ir)
         rfac = ir * hx * hym2 / (2.d0 * r)
         res(i*ir,j*jr) = res(i*ir,j*jr) - idir * fac0 *
     @     rfac * sigmac(is,js) * (cdst(i,j-jdir) - cdst(i,j))
      endif
      fac0 = fac0 / (ir * jr)
      hxm2 = ir * ir * hxm2
      hym2 = jr * jr * hym2
      i = i * ir
      j = j * jr
      if (idir .eq. 1) then
         is = i
      else
         is = i - 1
      endif
      if (jdir .eq. 1) then
         js = j
      else
         js = j - 1
      endif
      res(i,j) = res(i,j) - fac0 *
     @  (hxm2 * (sigmaf(is,js-jdir) *
     @            (fdst(i+idir,j-jdir) - fdst(i,j-jdir) +
     @             2.d0 * (fdst(i+idir,j) - fdst(i,j))) +
     @           sigmaf(is,js) *
     @            (fdst(i+idir,j+jdir) - fdst(i,j+jdir) +
     @             2.d0 * (fdst(i+idir,j) - fdst(i,j))) +
     @           sigmaf(is-idir,js) *
     @            (fdst(i-idir,j+jdir) - fdst(i,j+jdir) +
     @             2.d0 * (fdst(i-idir,j) - fdst(i,j)))) +
     @   hym2 * (sigmaf(is,js-jdir) *
     @            (fdst(i+idir,j-jdir) - fdst(i+idir,j) +
     @             2.d0 * (fdst(i,j-jdir) - fdst(i,j))) +
     @           sigmaf(is,js) *
     @            (fdst(i+idir,j+jdir) - fdst(i+idir,j) +
     @             2.d0 * (fdst(i,j+jdir) - fdst(i,j))) +
     @           sigmaf(is-idir,js) *
     @            (fdst(i-idir,j+jdir) - fdst(i-idir,j) +
     @             2.d0 * (fdst(i,j+jdir) - fdst(i,j)))))
      if (irz .eq. 1) then
         r = (is + 0.5d0) * hx
         rfac0 = hx * hym2 / (2.d0 * (r - idir * hx))
         rfac1 = hx * hym2 / (2.d0 * r)
         res(i,j) = res(i,j) - idir * fac0 *
     @     ((rfac0 * sigmaf(is-idir,js) - rfac1 * sigmaf(is,js)) *
     @      (fdst(i,j+jdir) - fdst(i,j)) -
     @      rfac1 * sigmaf(is,js-jdir) * (fdst(i,j-jdir) - fdst(i,j)))
      endif
      fac0 = fac0 / ir
      do 10 m = idir, idir*(ir-1), idir
         fac1 = (ir-abs(m)) * fac0
         tmp = hxm2 *
     @     (sigmaf(i-m-1,js) *
     @       (fdst(i-m-1,j+jdir) - fdst(i-m,j+jdir) +
     @        2.d0 * (fdst(i-m-1,j) - fdst(i-m,j))) +
     @      sigmaf(i-m,js) *
     @       (fdst(i-m+1,j+jdir) - fdst(i-m,j+jdir) +
     @        2.d0 * (fdst(i-m+1,j) - fdst(i-m,j))))
         res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     @     (sigmaf(i-m-1,js) *
     @       (fdst(i-m-1,j+jdir) - fdst(i-m-1,j) +
     @        2.d0 * (fdst(i-m,j+jdir) - fdst(i-m,j))) +
     @	    sigmaf(i-m,js) *
     @       (fdst(i-m+1,j+jdir) - fdst(i-m+1,j) +
     @        2.d0 * (fdst(i-m,j+jdir) - fdst(i-m,j)))))
         if (irz .eq. 1) then
            r = (i - m + 0.5d0) * hx
            rfac0m = hx * hym2 / (2.d0 * (r - hx))
            rfac1m = hx * hym2 / (2.d0 * r)
            res(i,j) = res(i,j) - fac1 *
     @         (rfac0m * sigmaf(i-m-1,js) - rfac1m * sigmaf(i-m,js)) *
     @          (fdst(i-m,j+jdir) - fdst(i-m,j))
         endif
 10      continue
      fac0 = ir * fac0 / jr
      do 20 n = jdir, jdir*(jr-1), jdir
         fac1 = (jr-abs(n)) * fac0
         tmp = hxm2 *
     @     (sigmaf(is,j-n-1) *
     @       (fdst(i+idir,j-n-1) - fdst(i,j-n-1) +
     @        2.d0 * (fdst(i+idir,j-n) - fdst(i,j-n))) +
     @      sigmaf(is,j-n) *
     @       (fdst(i+idir,j-n+1) - fdst(i,j-n+1) +
     @        2.d0 * (fdst(i+idir,j-n) - fdst(i,j-n))))
         res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     @     (sigmaf(is,j-n-1) *
     @       (fdst(i+idir,j-n-1) - fdst(i+idir,j-n) + 
     @        2.d0 * (fdst(i,j-n-1) - fdst(i,j-n))) +
     @	    sigmaf(is,j-n) *
     @       (fdst(i+idir,j-n+1) - fdst(i+idir,j-n) +
     @        2.d0 * (fdst(i,j-n+1) - fdst(i,j-n)))))
         if (irz .eq. 1) then
            r = (is + 0.5d0) * hx
            rfac1 = -idir * hx * hym2 / (2.d0 * r)
            res(i,j) = res(i,j) - fac1 *
     @        rfac1 * (sigmaf(is,j-n-1) * (fdst(i,j-n-1) - fdst(i,j-n)) +
     @                 sigmaf(is,j-n)   * (fdst(i,j-n+1) - fdst(i,j-n)))
         endif
 20      continue
      end

c Note---assumes fdst linearly interpolated from cdst along edge
      subroutine hgdres(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ sigmaf, dims(sf),
     @ sigmac, dims(sc),
     @ dims(reg),
#ifdef TURNOFF
     @ hx, hy, ir, jr, jdir, irz1)
#else
     @ hx, hy, ir, jr, jdir, irz)
#endif
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(sf)
      integer dims(sc)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 sigmaf(dimdec(sf))
      real*8 sigmac(dimdec(sc))
      real*8 hx, hy
      integer ir, jr, jdir, irz
#ifdef TURNOFF
      irz = 0
#endif
      i = regl0
      j = regl1
      if (jdir .eq. 1) then
         js = j
      else
         js = j - 1
      endif
      hxm2 = 1.d0 / (ir * ir * hx * hx)
      hym2 = 1.d0 / (jr * jr * hy * hy)
      fac0 = (ir * jr) / (3.d0*ir*jr + 3.d0*ir + 3.d0*jr - 3.d0)
      res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     @  (hxm2 * (sigmac(i,js-jdir) *
     @            (cdst(i+1,j-jdir) - cdst(i,j-jdir) +
     @             2.d0 * (cdst(i+1,j) - cdst(i,j))) +
     @           sigmac(i-1,js) *
     @            (cdst(i-1,j+jdir) - cdst(i,j+jdir) +
     @             2.d0 * (cdst(i-1,j) - cdst(i,j)))) +
     @   hym2 * (sigmac(i,js-jdir) *
     @            (cdst(i+1,j-jdir) - cdst(i+1,j) +
     @             2.d0 * (cdst(i,j-jdir) - cdst(i,j))) +
     @           sigmac(i-1,js) *
     @            (cdst(i-1,j+jdir) - cdst(i-1,j) +
     @             2.d0 * (cdst(i,j+jdir) - cdst(i,j)))))
      if (irz .eq. 1) then
         r = (i + 0.5d0) * (hx * ir)
         rfac0 = ir * hx * hym2 / (2.d0 * (r - hx * ir))
         rfac1 = ir * hx * hym2 / (2.d0 * r)
         res(i*ir,j*jr) = res(i*ir,j*jr) - fac0 *
     @     (rfac0 * sigmac(i-1,js)    * (cdst(i,j+jdir) - cdst(i,j)) -
     @      rfac1 * sigmac(i,js-jdir) * (cdst(i,j-jdir) - cdst(i,j)))
      endif
      fac0 = fac0 / (ir * jr)
      hxm2 = ir * ir * hxm2
      hym2 = jr * jr * hym2
      i = i * ir
      j = j * jr
      if (jdir .eq. 1) then
         js = j
      else
         js = j - 1
      endif
      res(i,j) = res(i,j) - fac0 *
     @  (hxm2 * (sigmaf(i,js) *
     @            (fdst(i+1,j+jdir) - fdst(i,j+jdir) +
     @             2.d0 * (fdst(i+1,j) - fdst(i,j))) +
     @           sigmaf(i-1,js-jdir) *
     @            (fdst(i-1,j-jdir) - fdst(i,j-jdir) +
     @             2.d0 * (fdst(i-1,j) - fdst(i,j)))) +
     @   hym2 * (sigmaf(i,js) *
     @            (fdst(i+1,j+jdir) - fdst(i+1,j) +
     @             2.d0 * (fdst(i,j+jdir) - fdst(i,j))) +
     @           sigmaf(i-1,js-jdir) *
     @            (fdst(i-1,j-jdir) - fdst(i-1,j) +
     @             2.d0 * (fdst(i,j-jdir) - fdst(i,j)))))
      if (irz .eq. 1) then
         r = (i + 0.5d0) * hx
         rfac0 = hx * hym2 / (2.d0 * (r - hx))
         rfac1 = hx * hym2 / (2.d0 * r)
         res(i,j) = res(i,j) - fac0 *
     @     (rfac0 * sigmaf(i-1,js-jdir) * (fdst(i,j-jdir) - fdst(i,j)) -
     @      rfac1 * sigmaf(i,js)        * (fdst(i,j+jdir) - fdst(i,j)))
      endif
      fac0 = fac0 / ir
      do 10 m = 1, ir-1
         fac1 = (ir-m) * fac0
         tmp = hxm2 *
     @     (sigmaf(i+m-1,js) *
     @       (fdst(i+m-1,j+jdir) - fdst(i+m,j+jdir) +
     @        2.d0 * (fdst(i+m-1,j) - fdst(i+m,j))) +
     @      sigmaf(i+m,js) *
     @       (fdst(i+m+1,j+jdir) - fdst(i+m,j+jdir) +
     @        2.d0 * (fdst(i+m+1,j) - fdst(i+m,j))) +
     @      sigmaf(i-m-1,js-jdir) *
     @       (fdst(i-m-1,j-jdir) - fdst(i-m,j-jdir) +
     @        2.d0 * (fdst(i-m-1,j) - fdst(i-m,j))) +
     @      sigmaf(i-m,js-jdir) *
     @       (fdst(i-m+1,j-jdir) - fdst(i-m,j-jdir) +
     @        2.d0 * (fdst(i-m+1,j) - fdst(i-m,j))))
         res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     @     (sigmaf(i+m-1,js) *
     @       (fdst(i+m-1,j+jdir) - fdst(i+m-1,j) +
     @        2.d0 * (fdst(i+m,j+jdir) - fdst(i+m,j))) +
     @      sigmaf(i+m,js) *
     @       (fdst(i+m+1,j+jdir) - fdst(i+m+1,j) +
     @        2.d0 * (fdst(i+m,j+jdir) - fdst(i+m,j))) +
     @      sigmaf(i-m-1,js-jdir) *
     @       (fdst(i-m-1,j-jdir) - fdst(i-m-1,j) +
     @        2.d0 * (fdst(i-m,j-jdir) - fdst(i-m,j))) +
     @      sigmaf(i-m,js-jdir) *
     @       (fdst(i-m+1,j-jdir) - fdst(i-m+1,j) +
     @        2.d0 * (fdst(i-m,j-jdir) - fdst(i-m,j)))))
         if (irz .eq. 1) then
            r = (i - m + 0.5d0) * hx
            rfac0m = hx * hym2 / (2.d0 * (r - hx))
            rfac1m = hx * hym2 / (2.d0 * r)
            r = (i + m + 0.5d0) * hx
            rfac0p = hx * hym2 / (2.d0 * (r - hx))
            rfac1p = hx * hym2 / (2.d0 * r)
            res(i,j) = res(i,j) - fac1 *
     @       ((rfac0m * sigmaf(i-m-1,js-jdir) - rfac1m * sigmaf(i-m,js-jdir)) *
     @          (fdst(i-m,j-jdir) - fdst(i-m,j)) +
     @        (rfac0p * sigmaf(i+m-1,js) - rfac1p * sigmaf(i+m,js)) *
     @          (fdst(i+m,j+jdir) - fdst(i+m,j)))
         endif
 10      continue
      fac0 = ir * fac0 / jr
      do 20 n = jdir, jdir*(jr-1), jdir
         fac1 = (jr-abs(n)) * fac0
         tmp = hxm2 *
     @     (sigmaf(i,j+n-1) *
     @       (fdst(i+1,j+n-1) - fdst(i,j+n-1) +
     @        2.d0 * (fdst(i+1,j+n) - fdst(i,j+n))) +
     @      sigmaf(i,j+n) *
     @       (fdst(i+1,j+n+1) - fdst(i,j+n+1) +
     @        2.d0 * (fdst(i+1,j+n) - fdst(i,j+n))) +
     @      sigmaf(i-1,j-n-1) *
     @       (fdst(i-1,j-n-1) - fdst(i,j-n-1) +
     @        2.d0 * (fdst(i-1,j-n) - fdst(i,j-n))) +
     @      sigmaf(i-1,j-n) *
     @       (fdst(i-1,j-n+1) - fdst(i,j-n+1) +
     @        2.d0 * (fdst(i-1,j-n) - fdst(i,j-n))))
         res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     @     (sigmaf(i,j+n-1) *
     @       (fdst(i+1,j+n-1) - fdst(i+1,j+n) +
     @        2.d0 * (fdst(i,j+n-1) - fdst(i,j+n))) +
     @      sigmaf(i,j+n) *
     @       (fdst(i+1,j+n+1) - fdst(i+1,j+n) +
     @        2.d0 * (fdst(i,j+n+1) - fdst(i,j+n))) +
     @      sigmaf(i-1,j-n-1) *
     @       (fdst(i-1,j-n-1) - fdst(i-1,j-n) +
     @        2.d0 * (fdst(i,j-n-1) - fdst(i,j-n))) +
     @      sigmaf(i-1,j-n) *
     @       (fdst(i-1,j-n+1) - fdst(i-1,j-n) +
     @        2.d0 * (fdst(i,j-n+1) - fdst(i,j-n)))))
         if (irz .eq. 1) then
            r = (i + 0.5d0) * hx
            rfac0 = hx * hym2 / (2.d0 * (r - hx))
            rfac1 = hx * hym2 / (2.d0 * r)
            res(i,j) = res(i,j) - fac1 *
     @       (rfac0 * (sigmaf(i-1,j-n-1) * (fdst(i,j-n-1) - fdst(i,j-n)) +
     @                 sigmaf(i-1,j-n)   * (fdst(i,j-n+1) - fdst(i,j-n))) -
     @        rfac1 * (sigmaf(i,j+n-1) * (fdst(i,j+n-1) - fdst(i,j+n)) +
     @                 sigmaf(i,j+n)   * (fdst(i,j+n+1) - fdst(i,j+n))))
         endif
 20      continue
      end

#else

c constant density versions:

c Note---assumes fdst linearly interpolated from cdst along edge
      subroutine hgfres(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ dims(reg),
     @ hx, irat, idim, idir)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 hx
      integer irat, idim, idir
      ir = irat
      jr = irat
      hm2 = 1.d0 / (hx*hx)
      fac0 = 1.d0 / (3.d0 * (irat*irat+irat))
      if (idim .eq. 0) then
         i = regl0
         do 10 j = regl1, regh1
 10         res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 * hm2 *
     @        (2.d0 * (cdst(i-idir,j-1) + cdst(i-idir,j+1) +
     @                 cdst(i-idir,j)) +
     @         (cdst(i,j-1) + cdst(i,j+1)) -
     @         8.d0 * cdst(i,j))
         fac0 = fac0 / irat
         i = i * irat
         do 20 n = 0, irat-1
            fac1 = (irat-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do 20 j = irat*regl1, irat*regh1, irat
 20            res(i,j) = res(i,j) - fac1 * hm2 *
     @           (2.d0 * (fdst(i+idir,j-n-1) + fdst(i+idir,j-n+1) +
     @                    fdst(i+idir,j+n-1) + fdst(i+idir,j+n+1) +
     @                    fdst(i+idir,j-n) + fdst(i+idir,j+n)) +
     @            (fdst(i,j-n-1) + fdst(i,j-n+1) +
     @             fdst(i,j+n-1) + fdst(i,j+n+1)) -
     @            8.d0 * (fdst(i,j-n) + fdst(i,j+n)))
      else
         j = regl1
         do 30 i = regl0, regh0
 30         res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 * hm2 *
     @        (2.d0 * (cdst(i-1,j-idir) + cdst(i+1,j-idir) +
     @                 cdst(i,j-idir)) +
     @         (cdst(i-1,j) + cdst(i+1,j)) -
     @         8.d0 * cdst(i,j))
         fac0 = fac0 / irat
         j = j * irat
         do 40 m = 0, irat-1
            fac1 = (irat-m) * fac0
            if (m .eq. 0) fac1 = 0.5d0 * fac1
            do 40 i = irat*regl0, irat*regh0, irat
 40            res(i,j) = res(i,j) - fac1 * hm2 *
     @           (2.d0 * (fdst(i-m-1,j+idir) + fdst(i-m+1,j+idir) +
     @                    fdst(i+m-1,j+idir) + fdst(i+m+1,j+idir) +
     @                    fdst(i-m,j+idir) + fdst(i+m,j+idir)) +
     @            (fdst(i-m-1,j) + fdst(i-m+1,j) +
     @             fdst(i+m-1,j) + fdst(i+m+1,j)) -
     @            8.d0 * (fdst(i-m,j) + fdst(i+m,j)))
      endif
      end

c Note---assumes fdst linearly interpolated from cdst along edge
      subroutine hgores(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ dims(reg),
     @ hx, irat, idir, jdir)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 hx
      integer irat, idir, jdir
      ir = irat
      jr = irat
      i = regl0
      j = regl1
      hm2 = 1.d0 / (hx*hx)
      fac0 = 1.d0 / ((4.5d0*irat+3.d0)*irat-1.5d0)
      res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 * hm2 *
     @  (2.d0 * (cdst(i-idir,j-jdir) + cdst(i-idir,j+jdir) +
     @           cdst(i+idir,j-jdir) +
     @           cdst(i-idir,j) + cdst(i,j-jdir)) +
     @   (cdst(i+idir,j) + cdst(i,j+jdir)) -
     @   12.d0 * cdst(i,j))
      i = i * irat
      j = j * irat
      res(i,j) = res(i,j) - fac0 * hm2 *
     @  (2.d0 * fdst(i+idir,j+jdir) +
     @   (fdst(i+idir,j) + fdst(i,j+jdir)) -
     @   4.d0 * fdst(i,j))
      fac0 = fac0 / irat
      do 10 n = jdir, jdir*(irat-1), jdir
         fac1 = (irat-abs(n)) * fac0
         m = idir*abs(n)
 10      res(i,j) = res(i,j) - fac1 * hm2 *
     @     (2.d0 * (fdst(i+m-1,j+jdir) + fdst(i+m+1,j+jdir) +
     @              fdst(i+idir,j+n-1) + fdst(i+idir,j+n+1) +
     @              fdst(i+m,j+jdir) + fdst(i+idir,j+n)) +
     @      (fdst(i+m-1,j) + fdst(i+m+1,j) +
     @       fdst(i,j+n-1) + fdst(i,j+n+1)) -
     @      8.d0 * (fdst(i+m,j) + fdst(i,j+n)))
      end

c Note---assumes fdst linearly interpolated from cdst along edge
      subroutine hgires(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ dims(reg),
     @ hx, irat, idir, jdir)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 hx
      integer irat, idir, jdir
      ir = irat
      jr = irat
      i = regl0
      j = regl1
      hm2 = 1.d0 / (hx*hx)
      fac0 = 1.d0 / ((1.5d0*irat+3.d0)*irat+1.5d0)
      res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 * hm2 *
     @  (2.d0 * cdst(i-idir,j-jdir) +
     @   (cdst(i-idir,j) + cdst(i,j-jdir)) -
     @   4.d0 * cdst(i,j))
      i = i * irat
      j = j * irat
      res(i,j) = res(i,j) - fac0 * hm2 *
     @  (2.d0 * (fdst(i-idir,j+jdir) + fdst(i+idir,j-jdir) +
     @           fdst(i+idir,j+jdir) +
     @           fdst(i+idir,j) + fdst(i,j+jdir)) +
     @   (fdst(i-idir,j) + fdst(i,j-jdir)) -
     @   12.d0 * fdst(i,j))
      fac0 = fac0 / irat
      do 10 n = jdir, jdir*(irat-1), jdir
         fac1 = (irat-abs(n)) * fac0
         m = idir*abs(n)
 10      res(i,j) = res(i,j) - fac1 * hm2 *
     @     (2.d0 * (fdst(i-m-1,j+jdir) + fdst(i-m+1,j+jdir) +
     @              fdst(i+idir,j-n-1) + fdst(i+idir,j-n+1) +
     @              fdst(i-m,j+jdir) + fdst(i+idir,j-n)) +
     @      (fdst(i-m-1,j) + fdst(i-m+1,j) +
     @       fdst(i,j-n-1) + fdst(i,j-n+1)) -
     @      8.d0 * (fdst(i-m,j) + fdst(i,j-n)))
      end

c Note---assumes fdst linearly interpolated from cdst along edge
      subroutine hgdres(res, dims(res),
     @ src,  dims(src),
     @ fdst, dims(fdst),
     @ cdst, dims(cdst),
     @ dims(reg),
     @ hx, irat, jdir)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(fdst)
      integer dims(cdst)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 fdst(dimdec(fdst))
      real*8 cdst(dimdec(cdst))
      real*8 hx
      integer irat, jdir
      ir = irat
      jr = irat
      i = regl0
      j = regl1
      hm2 = 1.d0 / (hx*hx)
      fac0 = 1.d0 / ((3.d0*irat+6.d0)*irat-3.d0)
      res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 * hm2 *
     @  (2.d0 * (cdst(i-1,j+jdir) + cdst(i+1,j-jdir)) +
     @   (cdst(i-1,j) + cdst(i+1,j) + cdst(i,j-1) + cdst(i,j+1)) -
     @   8.d0 * cdst(i,j))
      i = i * irat
      j = j * irat
      res(i,j) = res(i,j) - fac0 * hm2 *
     @  (2.d0 * (fdst(i+1,j+jdir) + fdst(i-1,j-jdir)) +
     @   (fdst(i-1,j) + fdst(i+1,j) + fdst(i,j-1) + fdst(i,j+1)) -
     @   8.d0 * fdst(i,j))
      fac0 = fac0 / irat
      do 10 n = jdir, jdir*(irat-1), jdir
         fac1 = (irat-abs(n)) * fac0
         m = abs(n)
 10      res(i,j) = res(i,j) - fac1 * hm2 *
     @     (2.d0 * (fdst(i-1,j-n-1) + fdst(i-1,j-n+1) +
     @              fdst(i+1,j+n-1) + fdst(i+1,j+n+1) +
     @              fdst(i-m-1,j-jdir) + fdst(i-m+1,j-jdir) +
     @              fdst(i+m-1,j+jdir) + fdst(i+m+1,j+jdir) +
     @              fdst(i-1,j-n) + fdst(i+1,j+n) +
     @              fdst(i-m,j-jdir) + fdst(i+m,j+jdir)) +
     @      (fdst(i,j-n-1) + fdst(i,j-n+1) +
     @       fdst(i,j+n-1) + fdst(i,j+n+1) +
     @       fdst(i-m-1,j) + fdst(i-m+1,j) +
     @       fdst(i+m-1,j) + fdst(i+m+1,j)) -
     @      8.d0 * (fdst(i,j-n) + fdst(i,j+n) +
     @              fdst(i-m,j) + fdst(i+m,j)))
      end

#endif

#endif

#ifndef CONSTANT

c NODE-based data, factor of 2 only.
      subroutine hgints(dest, dims(dest),
     @ dims(reg),
     @ sigx, sigy, dims(sb),
     @ src,  dims(src),
     @ dims(bb), ir, jr)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(dest)
      integer dims(reg)
      integer dims(sb)
      integer dims(src)
      integer dims(bb)
      integer ir, jr
      real*8 dest(dimdec(dest))
      real*8 sigx(dimdec(sb))
      real*8 sigy(dimdec(sb))
      real*8 src(dimdec(src))
      do 10 jc = bbl1, bbh1
         do 10 ic = bbl0, bbh0
 10         dest(ir*ic,jr*jc) = src(ic,jc)
      if (ir .eq. 2) then
         do 20 jc = bbl1, bbh1
            do 20 ic = bbl0, bbh0-1
               i = ir * ic
               j = jr * jc
 20            dest(i+1,j) =
     @           ((sigx(i,j-1)+sigx(i,j)) * src(ic,jc) +
     @            (sigx(i+1,j-1)+sigx(i+1,j)) * src(ic+1,jc)) /
     @           (sigx(i,j-1)+sigx(i,j)+
     @            sigx(i+1,j-1)+sigx(i+1,j))
      endif
      if (jr .eq. 2) then
         do 30 jc = bbl1, bbh1-1
            do 30 ic = bbl0, bbh0
               i = ir * ic
               j = jr * jc
 30            dest(i,j+1) =
     @           ((sigy(i-1,j)+sigy(i,j)) * src(ic,jc) +
     @            (sigy(i-1,j+1)+sigy(i,j+1)) * src(ic,jc+1)) /
     @           (sigy(i-1,j)+sigy(i,j)+
     @            sigy(i-1,j+1)+sigy(i,j+1))
      endif
      if (ir .eq. 2 .and. jr .eq. 2) then
         do 40 jc = bbl1, bbh1-1
cdir$ ivdep
            do 40 ic = bbl0, bbh0-1
               i = ir * ic
               j = jr * jc
 40            dest(i+1,j+1) = ((sigx(i,j) + sigx(i,j+1)) *
     @                            dest(i,j+1) +
     @                          (sigx(i+1,j) + sigx(i+1,j+1)) *
     @                            dest(i+2,j+1) +
     @                          (sigy(i,j) + sigy(i+1,j)) *
     @                            dest(i+1,j) +
     @                          (sigy(i,j+1) + sigy(i+1,j+1)) *
     @                            dest(i+1,j+2)) /
     @                         (sigx(i,j) + sigx(i,j+1) +
     @                          sigx(i+1,j) + sigx(i+1,j+1) +
     @                          sigy(i,j) + sigy(i+1,j) +
     @                          sigy(i,j+1) + sigy(i+1,j+1))
      endif
      end

c CELL-based data only.
      subroutine hgsrst(destx, desty,
     @ dims(dest),
     @ dims(reg),
     @ srcx, srcy, dims(src),
     @ ir, jr)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer ir, jr
      real*8 destx(dimdec(dest))
      real*8 desty(dimdec(dest))
      real*8 srcx(dimdec(src))
      real*8 srcy(dimdec(src))
      if (ir .eq. 2 .and. jr .eq. 2) then
         do 10 j = regl1, regh1
            do 10 i = regl0, regh0
               i2 = 2 * i
               j2 = 2 * j
               destx(i,j) = 1.d0 /
     @                      (1.d0 / (srcx(i2,j2)   + srcx(i2,j2+1)) +
     @                       1.d0 / (srcx(i2+1,j2) + srcx(i2+1,j2+1)))
               desty(i,j) = 1.d0 /
     @                      (1.d0 / (srcy(i2,j2)   + srcy(i2+1,j2)) +
     @                       1.d0 / (srcy(i2,j2+1) + srcy(i2+1,j2+1)))
 10            continue
      else if (ir .eq. 2) then
         do 20 j = regl1, regh1
            do 20 i = regl0, regh0
               i2 = 2 * i
               destx(i,j) = 2.d0 /
     @                      (1.d0 / srcx(i2,j) + 1.d0 / srcx(i2+1,j))
               desty(i,j) = 0.5d0 * (srcy(i2,j) + srcy(i2+1,j))
 20            continue
      else
         do 30 j = regl1, regh1
            do 30 i = regl0, regh0
               j2 = 2 * j
               destx(i,j) = 0.5d0 * (srcx(i,j2) + srcx(i,j2+1))
               desty(i,j) = 2.d0 /
     @                      (1.d0 / srcy(i,j2) + 1.d0 / srcy(i,j2+1))
 30            continue
      endif
      end

#endif

#ifdef CROSS_STENCIL

#ifdef CONSTANT

c five-point constant stencils

#define AVG fac * (cor(i-1) + cor(i+1) + cor(i-idiff) + cor(i+idiff))

      subroutine hgrlxu(cor,
     @ res, dims(res),
     @ mask,
     @ dims(reg),
     @ hx)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(reg)
      real*8 cor(dimdecu(res))
      real*8 res(dimdecu(res))
      real*8 mask(dimdecu(res))
      real*8 hx
#if 0
      avg() = fac * (cor(i-1) + cor(i+1) +
     @               cor(i-idiff) + cor(i+idiff))
#endif
      hxm2 = 1.d0 / (hx*hx)
      fac = hxm2
      cen = 0.25d0 / hxm2
      idiff = resh0 - resl0 + 1
cdir$ ivdep
      do 10 i = (regl1 - resl1) * idiff + (regl0 - resl0) + 1,
     @          (regh1 - resl1) * idiff + (regh0 - resl0) + 1, 2
 10      cor(i) = cor(i) + mask(i) * ((AVG - res(i)) * cen - cor(i))
cdir$ ivdep
      do 20 i = (regl1 - resl1) * idiff + (regl0 - resl0) + 2,
     @          (regh1 - resl1) * idiff + (regh0 - resl0) + 1, 2
 20      cor(i) = cor(i) + mask(i) * ((AVG - res(i)) * cen - cor(i))
      end

      subroutine hgres(res, dims(res),
     @ src, dest,
     @ dims(reg),
     @ hx)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(res))
      real*8 dest(dimdec(res))
      real*8 hx
      hxm2 = 1.d0 / (hx*hx)
      fac = hxm2
      if (regh1 - regl1 .gt. regh0 - regl0) then
         do 10 i = regl0, regh0
            do 10 j = regl1, regh1
 10            res(i,j) = src(i,j) - fac *
     @           (dest(i-1,j) + dest(i+1,j) +
     @            dest(i,j-1) + dest(i,j+1) -
     @            4.d0 * dest(i,j))
      else
         do 20 j = regl1, regh1
            do 20 i = regl0, regh0
 20            res(i,j) = src(i,j) - fac *
     @           (dest(i-1,j) + dest(i+1,j) +
     @            dest(i,j-1) + dest(i,j+1) -
     @            4.d0 * dest(i,j))
      endif
      end

      subroutine hgresu(res, dims(res),
     @ src, dest,
     @ dims(reg),
     @ hx)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(reg)
      real*8 res(dimdecu(res))
      real*8 src(dimdecu(res))
      real*8 dest(dimdecu(res))
      real*8 hx
      hxm2 = 1.d0 / (hx*hx)
      fac = hxm2
      idiff = resh0 - resl0 + 1
      do 10 i = (regl1 - resl1) * idiff + (regl0 - resl0) + 1,
     @          (regh1 - resl1) * idiff + (regh0 - resl0) + 1
 10      res(i) = src(i) - fac *
     @     (dest(i-1) + dest(i+1) +
     @      dest(i-idiff) + dest(i+idiff) -
     @      4.d0 * dest(i))
      end

#else

c five-point variable stencils

      subroutine hgcen(cen, dims(cen),
     @ sigx, sigy, dims(sb),
     @ dims(reg),
     @ hx, hy)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(cen)
      integer dims(sb)
      integer dims(reg)
      real*8 cen(dimdec(cen))
      real*8 sigx(dimdec(sb))
      real*8 sigy(dimdec(sb))
      real*8 hx, hy
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      do 10 j = regl1, regh1
         do 10 i = regl0, regh0
 10         cen(i,j) = 2.d0 / (hxm2 * (sigx(i-1,j-1) + sigx(i-1,j) +
     @                                 sigx(i,j-1) + sigx(i,j)) +
     @                         hym2 * (sigy(i-1,j-1) + sigy(i-1,j) +
     @                                 sigy(i,j-1) + sigy(i,j)))
      end

#define AVG 0.5d0 * (hxm2 * \
               ((sigx(i-1,j-1) + sigx(i-1,j)) * cor(i-1,j) + \
                (sigx(i,j-1) + sigx(i,j)) * cor(i+1,j)) + \
                      hym2 * \
               ((sigy(i-1,j-1) + sigy(i,j-1)) * cor(i,j-1) + \
                (sigy(i-1,j) + sigy(i,j)) * cor(i,j+1)))

      subroutine hgrlx(cor, dims(cor),
     @ res, dims(res),
     @ sigx, sigy, dims(sf),
     @ cen, dims(cen),
     @ dims(reg),
     @ hx, hy, irz, imax)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(cor)
      integer dims(res)
      integer dims(sf)
      integer dims(cen)
      integer dims(reg)
      real*8 cor(dimdec(cor))
      real*8 res(dimdec(res))
      real*8 sigx(dimdec(sf))
      real*8 sigy(dimdec(sf))
      real*8 cen(dimdec(cen))
      real*8 hx, hy
      integer irz, imax
#if 0
      avg() = 0.5d0 * (hxm2 *
     @          ((sigx(i-1,j-1) + sigx(i-1,j)) * cor(i-1,j) +
     @           (sigx(i,j-1) + sigx(i,j)) * cor(i+1,j)) +
     @                 hym2 *
     @          ((sigy(i-1,j-1) + sigy(i,j-1)) * cor(i,j-1) +
     @           (sigy(i-1,j) + sigy(i,j)) * cor(i,j+1)))
#endif
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      if (regh1 - regl1 .gt. regh0 - regl0) then
         do 10 ipass = 1, 0, -1
            ipar = ipass
            do 10 i = regl0, regh0
               ipar = 1 - ipar
cdir$ ivdep
               do 10 j = regl1 + ipar, regh1, 2
                  cor(i,j) = (AVG - res(i,j)) * cen(i,j)
 10               continue
      else
         do 20 ipass = 1, 0, -1
            ipar = ipass
            do 20 j = regl1, regh1
               ipar = 1 - ipar
cdir$ ivdep
               do 20 i = regl0 + ipar, regh0, 2
                  cor(i,j) = (AVG - res(i,j)) * cen(i,j)
 20               continue
      endif
      end

#define RHSL (res(i,j) - 0.5d0 * hxm2 * \
               ((sigx(i-1,j-1) + sigx(i-1,j)) * cor(i-1,j) + \
                (sigx(i,j-1) + sigx(i,j)) * cor(i+1,j)))

      subroutine hgrlxl(cor, dims(cor),
     @ res, dims(res),
     @ sigx, sigy, dims(sf),
     @ cen, dims(cen),
     @ dims(reg),
     @ dims(dom),
     @ hx, hy, irz, imax, lsd)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(cor)
      integer dims(res)
      integer dims(sf)
      integer dims(cen)
      integer dims(reg)
      integer dims(dom)
      real*8 cor(dimdec(cor))
      real*8 res(dimdec(res))
      real*8 sigx(dimdec(sf))
      real*8 sigy(dimdec(sf))
      real*8 cen(dimdec(cen))
      real*8 hx, hy
      real*8 wrk(256)
      integer irz, imax, lsd
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      if (lsd .eq. 1) then
         do 20 ipass = 0, 1
         do 20 i = regl0 + ipass, regh0, 2
            j = regl1
            betm = -cen(i,j)
            if (betm .eq. 0.d0) then
c dirichlet bdy:
               cor(i,j) = 0.d0
               aj = 0.d0
               wrk(1) = 0.d0
            else if (regl1 .eq. doml1) then
c neumann bdy:
               cor(i,j) = RHSL * betm
               aj = 0.5d0 * hym2 * (sigy(i-1,j) + sigy(i,j))
               wrk(1) = 2.d0 * aj * betm
            else
c interface to grid at same level:
               aj = 0.5d0 * hym2 * (sigy(i-1,j-1) + sigy(i,j-1))
               cor(i,j) = (RHSL - aj * cor(i,j-1)) * betm
               aj = 0.5d0 * hym2 * (sigy(i-1,j) + sigy(i,j))
               wrk(1) = aj * betm
            endif
c forward solve loop:
            do 10 j = regl1 + 1, regh1 - 1
               jw = j - regl1
               if (cen(i,j) .eq. 0.d0) then
                  betm = 0.d0
               else
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(jw))
               endif
               cor(i,j) = (RHSL - aj * cor(i,j-1)) * betm
               aj = 0.5d0 * hym2 * (sigy(i-1,j) + sigy(i,j))
               wrk(jw + 1) = aj * betm
 10            continue
            j = regh1
            jw = j - regl1
            if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
               cor(i,j) = 0.d0
            else if (regh1 .eq. domh1) then
c neumann bdy:
               aj = 2.d0 * aj
               betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(jw))
               cor(i,j) = (RHSL - aj * cor(i,j-1)) * betm
            else if (jw .gt. 0) then
c interface to grid at same level:
               betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(jw))
               cor(i,j) = RHSL - aj * cor(i,j-1)
               aj = 0.5d0 * hym2 * (sigy(i-1,j) + sigy(i,j))
               cor(i,j) = (cor(i,j) - aj * cor(i,j+1)) * betm
            endif
c back substitution loop:
            do 20 j = regh1 - 1, regl1, -1
               jw = j - regl1
               cor(i,j) = cor(i,j) - wrk(jw + 1) * cor(i,j+1)
 20            continue
      else
         print *, "Line solve not implemented in dimension", lsd
      endif
      end

      subroutine hgrlnf(cor, dims(cor),
     @ res, dims(res),
     @ wrk, dims(wrk),
     @ sigx, sigy, dims(sf),
     @ cen, dims(cen),
     @ dims(reg),
     @ dims(dom),
     @ hx, hy, irz, imax, lsd, ipass)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(cor)
      integer dims(res)
      integer dims(wrk)
      integer dims(sf)
      integer dims(cen)
      integer dims(reg)
      integer dims(dom)
      real*8 cor(dimdec(cor))
      real*8 res(dimdec(res))
      real*8 wrk(dimdec(wrk))
      real*8 sigx(dimdec(sf))
      real*8 sigy(dimdec(sf))
      real*8 cen(dimdec(cen))
      real*8 hx, hy
      integer irz, imax, lsd, ipass
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      if (lsd .eq. 1) then
         if (mod(regl0,2) .eq. 0) then
            ioff = ipass
         else
            ioff = 1 - ipass
         endif
         j = regl1
         aj = 0.d0
         do 10 i = regl0 + ioff, regh0, 2
            if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
               wrk(i,j) = 0.d0
            else if (regl1 .eq. doml1) then
c neumann bdy:
               betm = -cen(i,j)
               cor(i,j) = RHSL * betm
               aj = 0.5d0 * hym2 * (sigy(i-1,j) + sigy(i,j))
               wrk(i,j) = 2.d0 * aj * betm
            else
c interface to grid at same level:
               aj = 0.5d0 * hym2 * (sigy(i-1,j) + sigy(i,j))
            endif
            wrk(i,j+1) = aj
 10         continue
c forward solve loop:
         do 20 j = regl1 + 1, regh1 - 1
            do 20 i = regl0 + ioff, regh0, 2
               if (cen(i,j) .eq. 0.d0) then
                  betm = 0.d0
               else
                  aj = wrk(i,j)
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
               endif
               cor(i,j) = (RHSL - aj * cor(i,j-1)) * betm
               aj = 0.5d0 * hym2 * (sigy(i-1,j) + sigy(i,j))
               wrk(i,j+1) = aj
 20            wrk(i,j) = aj * betm
         j = regh1
         do 30 i = regl0 + ioff, regh0, 2
            if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
            else if (regh1 .eq. domh1) then
c neumann bdy:
               aj = 2.d0 * wrk(i,j)
               betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
               cor(i,j) = (RHSL - aj * cor(i,j-1)) * betm
            else if (j .gt. regl1) then
c interface to grid at same level:
               aj = wrk(i,j)
               betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
               cor(i,j) = (RHSL - aj * cor(i,j-1)) * betm
               aj = 0.5d0 * hym2 * (sigy(i-1,j) + sigy(i,j))
               wrk(i,j) = aj * betm
            endif
 30         continue
      else
         print *, "Line solve not implemented in dimension", lsd
      endif
      end

      subroutine hgres(res, dims(res),
     @ src,  dims(src),
     @ dest, dims(dest),
     @ sigx, sigy, dims(sf),
     @ dims(reg),
     @ hx, hy)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(dest)
      integer dims(sf)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 dest(dimdec(dest))
      real*8 sigx(dimdec(sf))
      real*8 sigy(dimdec(sf))
      real*8 hx, hy
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      if (regh1 - regl1 .gt. regh0 - regl0) then
         do 10 i = regl0, regh0
            do 10 j = regl1, regh1
 10            res(i,j) = src(i,j) - 0.5d0 * (hxm2 *
     @           ((sigx(i-1,j-1) + sigx(i-1,j)) *
     @             (dest(i-1,j) - dest(i,j)) +
     @            (sigx(i,j-1) + sigx(i,j)) *
     @             (dest(i+1,j) - dest(i,j))) +
     @                                        hym2 *
     @           ((sigy(i-1,j-1) + sigy(i,j-1)) *
     @             (dest(i,j-1) - dest(i,j)) +
     @            (sigy(i-1,j) + sigy(i,j)) *
     @             (dest(i,j+1) - dest(i,j))))
      else
         do 20 j = regl1, regh1
            do 20 i = regl0, regh0
 20            res(i,j) = src(i,j) - 0.5d0 * (hxm2 *
     @           ((sigx(i-1,j-1) + sigx(i-1,j)) *
     @             (dest(i-1,j) - dest(i,j)) +
     @            (sigx(i,j-1) + sigx(i,j)) *
     @             (dest(i+1,j) - dest(i,j))) +
     @                                        hym2 *
     @           ((sigy(i-1,j-1) + sigy(i,j-1)) *
     @             (dest(i,j-1) - dest(i,j)) +
     @            (sigy(i-1,j) + sigy(i,j)) *
     @             (dest(i,j+1) - dest(i,j))))
      endif
      end

#endif

#else

#ifdef CONSTANT

c nine-point constant stencils

#define AVG fac * (cor(i-1,j-1) + cor(i-1,j) + cor(i-1,j+1) + \
                    cor(i,j-1) + cor(i,j+1) + \
                    cor(i+1,j-1) + cor(i+1,j) + cor(i+1,j+1))

      subroutine hgrlx(cor,
     @ res, dims(res),
     @ mask,
     @ dims(reg),
     @ hx)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(reg)
      real*8 cor(dimdec(res))
      real*8 res(dimdec(res))
      real*8 mask(dimdec(res))
      real*8 hx
#if 0
      avg() = fac * (cor(i-1,j-1) + cor(i-1,j) + cor(i-1,j+1) +
     @               cor(i,j-1) + cor(i,j+1) +
     @               cor(i+1,j-1) + cor(i+1,j) + cor(i+1,j+1))
#endif
      hxm2 = 1.d0 / (hx*hx)
      fac = hxm2 / 3.d0
      cen = 0.375d0 / hxm2
      if (regh1 - regl1 .gt. regh0 - regl0) then
         do 10 i = regl0, regh0
cdir$ ivdep
            do 10 j = regl1, regh1
 10            cor(i,j) = cor(i,j) +
     @           mask(i,j) * ((AVG - res(i,j)) * cen - cor(i,j))
      else
         do 20 j = regl1, regh1
cdir$ ivdep
            do 20 i = regl0, regh0
 20            cor(i,j) = cor(i,j) +
     @           mask(i,j) * ((AVG - res(i,j)) * cen - cor(i,j))
      endif
      end

      subroutine hgres(res, dims(res),
     @ src, dest,
     @ dims(reg),
     @ hx)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(res))
      real*8 dest(dimdec(res))
      real*8 hx
      hxm2 = 1.d0 / (hx*hx)
      fac = hxm2 / 3.d0
      if (regh1 - regl1 .gt. regh0 - regl0) then
         do 10 i = regl0, regh0
            do 10 j = regl1, regh1
 10            res(i,j) = src(i,j) - fac *
     @           (dest(i-1,j-1) + dest(i-1,j) + dest(i-1,j+1) +
     @            dest(i,j-1) + dest(i,j+1) +
     @            dest(i+1,j-1) + dest(i+1,j) + dest(i+1,j+1) -
     @            8.d0 * dest(i,j))
      else
         do 20 j = regl1, regh1
            do 20 i = regl0, regh0
 20            res(i,j) = src(i,j) - fac *
     @           (dest(i-1,j-1) + dest(i-1,j) + dest(i-1,j+1) +
     @            dest(i,j-1) + dest(i,j+1) +
     @            dest(i+1,j-1) + dest(i+1,j) + dest(i+1,j+1) -
     @            8.d0 * dest(i,j))
      endif
      end

      subroutine hgresu(res, dims(res),
     @ src, dest,
     @ dims(reg),
     @ hx)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(reg)
      real*8 res(dimdecu(res))
      real*8 src(dimdecu(res))
      real*8 dest(dimdecu(res))
      real*8 hx
      hxm2 = 1.d0 / (hx*hx)
      fac = hxm2 / 3.d0
      idiff = resh0 - resl0 + 1
      do 10 i = (regl1 - resl1) * idiff + (regl0 - resl0) + 1,
     @          (regh1 - resl1) * idiff + (regh0 - resl0) + 1
 10      res(i) = src(i) - fac *
     @     (dest(i-1-idiff) + dest(i-1) + dest(i-1+idiff) +
     @      dest(i-idiff) + dest(i+idiff) +
     @      dest(i+1-idiff) + dest(i+1) + dest(i+1+idiff) -
     @      8.d0 * dest(i))
      end

#else

c nine-point variable stencils

      subroutine hgcen(cen, dims(cen),
     @ sigx, sigy, dims(sb),
     @ dims(reg),
#ifdef TURNOFF
     @ hx, hy, irz1, imax)
#else
     @ hx, hy, irz, imax)
#endif
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(cen)
      integer dims(sb)
      integer dims(reg)
      real*8 cen(dimdec(cen))
      real*8 sigx(dimdec(sb))
      real*8 sigy(dimdec(sb))
      real*8 hx, hy
      integer irz, imax
#ifdef TURNOFF
      irz = 0
#endif
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      if (irz .eq. 0) then
         do 10 j = regl1, regh1
            do 10 i = regl0, regh0
 10            cen(i,j) = 3.d0 / (hxm2 * (sigx(i-1,j-1) + sigx(i-1,j) +
     @                                    sigx(i,j-1) + sigx(i,j)) +
     @                            hym2 * (sigy(i-1,j-1) + sigy(i-1,j) +
     @                                    sigy(i,j-1) + sigy(i,j)))
      else
         facrz = 0.25d0 * hx * hym2
         do 20 i = regl0, regh0
            r1 = (i + 0.5d0) * hx
            r0 = r1 - hx
            if (i .eq. imax) then
               r1 = -r0
            endif
            do 20 j = regl1, regh1
 20            cen(i,j) = 3.d0 / (hxm2 * (sigx(i-1,j-1) + sigx(i-1,j) +
     @                                    sigx(i,j-1) + sigx(i,j)) +
     @                            hym2 * (sigy(i-1,j-1) + sigy(i-1,j) +
     @                                    sigy(i,j-1) + sigy(i,j)) +
     @                            facrz * ((sigy(i-1,j-1) + sigy(i-1,j)) / r0 -
     @                                     (sigy(i,j-1)   + sigy(i,j))   / r1))
      endif
      end

#define AVG fac * (hxm2 * \
               (sigx(i-1,j-1) * \
                 (cor(i-1,j-1) - cor(i,j-1) + 2.d0 * cor(i-1,j)) + \
                sigx(i-1,j) * \
                 (cor(i-1,j+1) - cor(i,j+1) + 2.d0 * cor(i-1,j)) + \
                sigx(i,j-1) * \
                 (cor(i+1,j-1) - cor(i,j-1) + 2.d0 * cor(i+1,j)) + \
                sigx(i,j) * \
                 (cor(i+1,j+1) - cor(i,j+1) + 2.d0 * cor(i+1,j))) + \
                    hym2 * \
               (sigy(i-1,j-1) * \
                 (cor(i-1,j-1) - cor(i-1,j) + 2.d0 * cor(i,j-1)) + \
                sigy(i-1,j) * \
                 (cor(i-1,j+1) - cor(i-1,j) + 2.d0 * cor(i,j+1)) + \
                sigy(i,j-1) * \
                 (cor(i+1,j-1) - cor(i+1,j) + 2.d0 * cor(i,j-1)) + \
                sigy(i,j) * \
                 (cor(i+1,j+1) - cor(i+1,j) + 2.d0 * cor(i,j+1))))

#define AVGRZ AVG + facrz * \
                 ((sigy(i-1,j-1) / r0 - sigy(i,j-1) / r1) * cor(i,j-1) + \
                  (sigy(i-1,j)   / r0 - sigy(i,j)   / r1) * cor(i,j+1))

      subroutine hgrlx(cor, dims(cor),
     @ res, dims(res),
     @ sigx, sigy, dims(sf),
     @ cen, dims(cen),
     @ dims(reg),
#ifdef TURNOFF
     @ hx, hy, irz1, imax)
#else
     @ hx, hy, irz, imax)
#endif
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(cor)
      integer dims(res)
      integer dims(sf)
      integer dims(cen)
      integer dims(reg)
      real*8 cor(dimdec(cor))
      real*8 res(dimdec(res))
      real*8 sigx(dimdec(sf))
      real*8 sigy(dimdec(sf))
      real*8 cen(dimdec(cen))
      real*8 hx, hy
      integer irz, imax
#if 0
      avg() = fac * (hxm2 *
     @          (sigx(i-1,j-1) *
     @            (cor(i-1,j-1) - cor(i,j-1) + 2.d0 * cor(i-1,j)) +
     @           sigx(i-1,j) *
     @            (cor(i-1,j+1) - cor(i,j+1) + 2.d0 * cor(i-1,j)) +
     @           sigx(i,j-1) *
     @            (cor(i+1,j-1) - cor(i,j-1) + 2.d0 * cor(i+1,j)) +
     @           sigx(i,j) *
     @            (cor(i+1,j+1) - cor(i,j+1) + 2.d0 * cor(i+1,j))) +
     @               hym2 *
     @          (sigy(i-1,j-1) *
     @            (cor(i-1,j-1) - cor(i-1,j) + 2.d0 * cor(i,j-1)) +
     @           sigy(i-1,j) *
     @            (cor(i-1,j+1) - cor(i-1,j) + 2.d0 * cor(i,j+1)) +
     @           sigy(i,j-1) *
     @            (cor(i+1,j-1) - cor(i+1,j) + 2.d0 * cor(i,j-1)) +
     @           sigy(i,j) *
     @            (cor(i+1,j+1) - cor(i+1,j) + 2.d0 * cor(i,j+1))))
      avgrz() = avg() + facrz *
     @            ((sigy(i-1,j-1) / r0 - sigy(i,j-1) / r1) * cor(i,j-1) +
     @             (sigy(i-1,j)   / r0 - sigy(i,j)   / r1) * cor(i,j+1))
#endif
#ifdef TURNOFF
      irz = 0
#endif
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      fac = 1.d0 / 6.d0
      if (irz .eq. 0) then
         if (regh1 - regl1 .gt. regh0 - regl0) then
            do 10 i = regl0, regh0
cdir$ ivdep
               do 10 j = regl1, regh1
                  cor(i,j) = (AVG - res(i,j)) * cen(i,j)
 10               continue
         else
            do 20 j = regl1, regh1
cdir$ ivdep
               do 20 i = regl0, regh0
                  cor(i,j) = (AVG - res(i,j)) * cen(i,j)
 20               continue
         endif
      else
         facrz = hx * hym2 / 12.d0
         if (regh1 - regl1 .gt. regh0 - regl0) then
            do 30 i = regl0, regh0
               r1 = (i + 0.5d0) * hx
               r0 = r1 - hx
               if (i .eq. imax) then
                  r1 = -r0
               endif
cdir$ ivdep
               do 30 j = regl1, regh1
                  cor(i,j) = (AVGRZ - res(i,j)) * cen(i,j)
 30               continue
         else
            do 40 j = regl1, regh1
cdir$ ivdep
               do 40 i = regl0, regh0
                  r1 = (i + 0.5d0) * hx
                  r0 = r1 - hx
                  if (i .eq. imax) then
                     r1 = -r0
                  endif
                  cor(i,j) = (AVGRZ - res(i,j)) * cen(i,j)
 40               continue
         endif
      endif
      end

#define RHSL (res(i,j) - fac * (hxm2 * \
               (sigx(i-1,j-1) * \
                 (cor(i-1,j-1) + 2.d0 * cor(i-1,j)) + \
                sigx(i-1,j) * \
                 (cor(i-1,j+1) + 2.d0 * cor(i-1,j)) + \
                sigx(i,j-1) * \
                 (cor(i+1,j-1) + 2.d0 * cor(i+1,j)) + \
                sigx(i,j) * \
                 (cor(i+1,j+1) + 2.d0 * cor(i+1,j))) + \
                                hym2 * \
               (sigy(i-1,j-1) * \
                 (cor(i-1,j-1) - cor(i-1,j)) + \
                sigy(i-1,j) * \
                 (cor(i-1,j+1) - cor(i-1,j)) + \
                sigy(i,j-1) * \
                 (cor(i+1,j-1) - cor(i+1,j)) + \
                sigy(i,j) * \
                 (cor(i+1,j+1) - cor(i+1,j)))))

#define AVGRZL AVGL + facrz * \
                 ((sigy(i-1,j-1) / r0 - sigy(i,j-1) / r1) * cor(i,j-1) + \
                  (sigy(i-1,j)   / r0 - sigy(i,j)   / r1) * cor(i,j+1))

      subroutine hgrlxl(cor, dims(cor),
     @ res, dims(res),
     @ sigx, sigy, dims(sf),
     @ cen, dims(cen),
     @ dims(reg),
     @ dims(dom),
#ifdef TURNOFF
     @ hx, hy, irz1, imax, lsd)
#else
     @ hx, hy, irz, imax, lsd)
#endif
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(cor)
      integer dims(res)
      integer dims(sf)
      integer dims(cen)
      integer dims(reg)
      integer dims(dom)
      real*8 cor(dimdec(cor))
      real*8 res(dimdec(res))
      real*8 sigx(dimdec(sf))
      real*8 sigy(dimdec(sf))
      real*8 cen(dimdec(cen))
      real*8 hx, hy
      real*8 wrk(256)
      integer irz, imax, lsd
#ifdef TURNOFF
      irz = 0
#endif
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      fac = 1.d0 / 6.d0
      if (irz .eq. 0) then
         if (lsd .eq. 1) then
            do 30 ipass = 0, 1
            do 30 i = regl0 + ipass, regh0, 2
               j = regl1
               betm = -cen(i,j)
               if (betm .eq. 0.d0) then
c dirichlet bdy:
                  cor(i,j) = 0.d0
                  wrk(1) = 0.d0
               else if (regl1 .eq. doml1) then
c neumann bdy:
                  cor(i,j) = RHSL * betm
                  aj = fac *
     @              (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     @               hxm2 * (sigx(i-1,j) + sigx(i,j)))
                  wrk(1) = 2.d0 * aj * betm
               else
c interface to grid at same level:
                  aj = fac *
     @              (hym2 * 2.d0 * (sigy(i-1,j-1) + sigy(i,j-1)) -
     @               hxm2 * (sigx(i-1,j-1) + sigx(i,j-1)))
                  cor(i,j) = (RHSL - aj * cor(i,j-1)) * betm
                  aj = fac *
     @              (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     @               hxm2 * (sigx(i-1,j) + sigx(i,j)))
                  wrk(1) = aj * betm
               endif
c forward solve loop:
               do 10 j = regl1 + 1, regh1 - 1
                  jw = j - regl1
                  if (cen(i,j) .eq. 0.d0) then
                     betm = 0.d0
                  else
                     betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(jw))
                  endif
                  cor(i,j) = (RHSL - aj * cor(i,j-1)) * betm
                  aj = fac *
     @              (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     @               hxm2 * (sigx(i-1,j) + sigx(i,j)))
                  wrk(jw + 1) = aj * betm
 10               continue
               j = regh1
               jw = j - regl1
               if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
                  cor(i,j) = 0.d0
               else if (regh1 .eq. domh1) then
c neumann bdy:
                  aj = 2.d0 * aj
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(jw))
                  cor(i,j) = (RHSL - aj * cor(i,j-1)) * betm
               else if (jw .gt. 0) then
c interface to grid at same level:
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(jw))
                  cor(i,j) = RHSL - aj * cor(i,j-1)
                  aj = fac *
     @              (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     @               hxm2 * (sigx(i-1,j) + sigx(i,j)))
                  cor(i,j) = (cor(i,j) - aj * cor(i,j+1)) * betm
               endif
c back substitution loop:
               do 20 j = regh1 - 1, regl1, -1
                  jw = j - regl1
                  cor(i,j) = cor(i,j) - wrk(jw + 1) * cor(i,j+1)
 20               continue
c update neumann bdys:
               if (regl1 .eq. doml1) then
                  cor(i,regl1-1) = cor(i,regl1+1)
               endif
               if (regh1 .eq. domh1) then
                  cor(i,regh1+1) = cor(i,regh1-1)
               endif
 30            continue
         else
            print *, "Line solve not implemented in dimension", lsd
         endif
      else
         print *, "Line solve not implemented for rz:"
         print *, "   Using Gauss-Seidel instead."
         facrz = hx * hym2 / 12.d0
         if (regh1 - regl1 .gt. regh0 - regl0) then
            do 40 i = regl0, regh0
               r1 = (i + 0.5d0) * hx
               r0 = r1 - hx
               if (i .eq. imax) then
                  r1 = -r0
               endif
cdir$ ivdep
               do 40 j = regl1, regh1
                  cor(i,j) = (AVGRZ - res(i,j)) * cen(i,j)
 40               continue
         else
            do 50 j = regl1, regh1
cdir$ ivdep
               do 50 i = regl0, regh0
                  r1 = (i + 0.5d0) * hx
                  r0 = r1 - hx
                  if (i .eq. imax) then
                     r1 = -r0
                  endif
                  cor(i,j) = (AVGRZ - res(i,j)) * cen(i,j)
 50               continue
         endif
      endif
      end

      subroutine hgrlnf(cor, dims(cor),
     @ res, dims(res),
     @ wrk, dims(wrk),
     @ sigx, sigy, dims(sf),
     @ cen, dims(cen),
     @ dims(reg),
     @ dims(dom),
#ifdef TURNOFF
     @ hx, hy, irz1, imax, lsd, ipass)
#else
     @ hx, hy, irz, imax, lsd, ipass)
#endif
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(cor)
      integer dims(res)
      integer dims(wrk)
      integer dims(sf)
      integer dims(cen)
      integer dims(reg)
      integer dims(dom)
      real*8 cor(dimdec(cor))
      real*8 res(dimdec(res))
      real*8 wrk(dimdec(wrk))
      real*8 sigx(dimdec(sf))
      real*8 sigy(dimdec(sf))
      real*8 cen(dimdec(cen))
      real*8 hx, hy
      integer irz, imax, lsd, ipass
#ifdef TURNOFF
      irz = 0
#endif
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      fac = 1.d0 / 6.d0
      if (irz .eq. 0) then
         if (lsd .eq. 1) then
            if (mod(regl0,2) .eq. 0) then
               ioff = ipass
            else
               ioff = 1 - ipass
            endif
            j = regl1
            aj = 0.d0
            do 10 i = regl0 + ioff, regh0, 2
               if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
                  wrk(i,j) = 0.d0
               else if (regl1 .eq. doml1) then
c neumann bdy:
                  betm = -cen(i,j)
                  cor(i,j) = RHSL * betm
                  aj = fac *
     @              (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     @               hxm2 * (sigx(i-1,j) + sigx(i,j)))
                  wrk(i,j) = 2.d0 * aj * betm
               else
c interface to grid at same level:
                  aj = fac *
     @              (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     @               hxm2 * (sigx(i-1,j) + sigx(i,j)))
               endif
               wrk(i,j+1) = aj
 10            continue
c forward solve loop:
            do 20 j = regl1 + 1, regh1 - 1
               do 20 i = regl0 + ioff, regh0, 2
                  if (cen(i,j) .eq. 0.d0) then
                     betm = 0.d0
                  else
                     aj = wrk(i,j)
                     betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
                  endif
                  cor(i,j) = (RHSL - aj * cor(i,j-1)) * betm
                  aj = fac *
     @              (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     @               hxm2 * (sigx(i-1,j) + sigx(i,j)))
                  wrk(i,j+1) = aj
 20               wrk(i,j) = aj * betm
            j = regh1
            do 30 i = regl0 + ioff, regh0, 2
               if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
               else if (regh1 .eq. domh1) then
c neumann bdy:
                  aj = 2.d0 * wrk(i,j)
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
                  cor(i,j) = (RHSL - aj * cor(i,j-1)) * betm
               else if (j .gt. regl1) then
c interface to grid at same level:
                  aj = wrk(i,j)
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
                  cor(i,j) = (RHSL - aj * cor(i,j-1)) * betm
                  aj = fac *
     @              (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     @               hxm2 * (sigx(i-1,j) + sigx(i,j)))
                  wrk(i,j) = aj * betm
               endif
 30            continue
         else
            print *, "Line solve not implemented in dimension", lsd
         endif
      else
         print *, "Line solve not implemented for rz:"
      endif
      end

      subroutine hgres(res, dims(res),
     @ src,  dims(src),
     @ dest, dims(dest),
     @ sigx, sigy, dims(sf),
     @ dims(reg),
#ifdef TURNOFF
     @ hx, hy, irz1, imax)
#else
     @ hx, hy, irz, imax)
#endif
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(res)
      integer dims(src)
      integer dims(dest)
      integer dims(sf)
      integer dims(reg)
      real*8 res(dimdec(res))
      real*8 src(dimdec(src))
      real*8 dest(dimdec(dest))
      real*8 sigx(dimdec(sf))
      real*8 sigy(dimdec(sf))
      real*8 hx, hy
      integer irz, imax
#ifdef TURNOFF
      irz = 0
#endif
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      fac = 1.d0 / 6.d0
      if (regh1 - regl1 .gt. regh0 - regl0) then
         do 10 i = regl0, regh0
            do 10 j = regl1, regh1
               tmp = hxm2 *
     @           (sigx(i-1,j-1) *
     @             (dest(i-1,j-1) - dest(i,j-1) +
     @              2.d0 * (dest(i-1,j) - dest(i,j))) +
     @            sigx(i-1,j) *
     @             (dest(i-1,j+1) - dest(i,j+1) +
     @              2.d0 * (dest(i-1,j) - dest(i,j))) +
     @            sigx(i,j-1) *
     @             (dest(i+1,j-1) - dest(i,j-1) +
     @              2.d0 * (dest(i+1,j) - dest(i,j))) +
     @            sigx(i,j) *
     @             (dest(i+1,j+1) - dest(i,j+1) +
     @              2.d0 * (dest(i+1,j) - dest(i,j))))
 10            res(i,j) = src(i,j) - fac * (tmp + hym2 *
     @           (sigy(i-1,j-1) *
     @             (dest(i-1,j-1) - dest(i-1,j) +
     @              2.d0 * (dest(i,j-1) - dest(i,j))) +
     @            sigy(i-1,j) *
     @             (dest(i-1,j+1) - dest(i-1,j) +
     @              2.d0 * (dest(i,j+1) - dest(i,j))) +
     @            sigy(i,j-1) *
     @             (dest(i+1,j-1) - dest(i+1,j) +
     @              2.d0 * (dest(i,j-1) - dest(i,j))) +
     @            sigy(i,j) *
     @             (dest(i+1,j+1) - dest(i+1,j) +
     @              2.d0 * (dest(i,j+1) - dest(i,j)))))
      else
         do 20 j = regl1, regh1
            do 20 i = regl0, regh0
               tmp = hxm2 *
     @           (sigx(i-1,j-1) *
     @             (dest(i-1,j-1) - dest(i,j-1) +
     @              2.d0 * (dest(i-1,j) - dest(i,j))) +
     @            sigx(i-1,j) *
     @             (dest(i-1,j+1) - dest(i,j+1) +
     @              2.d0 * (dest(i-1,j) - dest(i,j))) +
     @            sigx(i,j-1) *
     @             (dest(i+1,j-1) - dest(i,j-1) +
     @              2.d0 * (dest(i+1,j) - dest(i,j))) +
     @            sigx(i,j) *
     @             (dest(i+1,j+1) - dest(i,j+1) +
     @              2.d0 * (dest(i+1,j) - dest(i,j))))
 20            res(i,j) = src(i,j) - fac * (tmp + hym2 *
     @           (sigy(i-1,j-1) *
     @             (dest(i-1,j-1) - dest(i-1,j) +
     @              2.d0 * (dest(i,j-1) - dest(i,j))) +
     @            sigy(i-1,j) *
     @             (dest(i-1,j+1) - dest(i-1,j) +
     @              2.d0 * (dest(i,j+1) - dest(i,j))) +
     @            sigy(i,j-1) *
     @             (dest(i+1,j-1) - dest(i+1,j) +
     @              2.d0 * (dest(i,j-1) - dest(i,j))) +
     @            sigy(i,j) *
     @             (dest(i+1,j+1) - dest(i+1,j) +
     @              2.d0 * (dest(i,j+1) - dest(i,j)))))
      endif
      if (irz .eq. 1) then
         fac = hx * hym2 / 12.d0
         do 30 i = regl0, regh0
            r1 = (i + 0.5d0) * hx
            r0 = r1 - hx
            if (i .eq. imax) then
               r1 = -r0
            endif
            do 30 j = regl1, regh1
 30            res(i,j) = res(i,j) - fac *
     @          ((sigy(i-1,j-1) * (dest(i,j-1) - dest(i,j)) +
     @            sigy(i-1,j)   * (dest(i,j+1) - dest(i,j))) / r0 -
     @           (sigy(i,j-1)   * (dest(i,j-1) - dest(i,j)) +
     @            sigy(i,j)     * (dest(i,j+1) - dest(i,j))) / r1)
      endif
      end

#endif

#endif

      subroutine hgrlnb(cor, dims(cor),
     @ wrk, dims(wrk),
     @ dims(reg),
     @ lsd, ipass)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(cor)
      integer dims(wrk)
      integer dims(reg)
      real*8 cor(dimdec(cor))
      real*8 wrk(dimdec(wrk))
      integer lsd, ipass
      if (lsd .eq. 1) then
         if (mod(regl0,2) .eq. 0) then
            ioff = ipass
         else
            ioff = 1 - ipass
         endif
c back substitution loop:
         do 20 j = regh1 - 1, regl1, -1
            do 20 i = regl0 + ioff, regh0, 2
               cor(i,j) = cor(i,j) - wrk(i,j) * cor(i,j+1)
 20            continue
      else
         print *, "Line solve not implemented in dimension", lsd
      endif
      end

#ifdef CONSTANT
#  define CGOPT 2
#else
#  define CGOPT 1
#endif

#if (CGOPT == 1)

c Unrolled indexing in these 3 routines uses the fact that each array
c has a border of width 1

c Works for NODE-based data.
      subroutine hgip(v0, v1, mask,
     @ dims(reg), sum)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(reg)
      real*8 v0(dimdecu(reg))
      real*8 v1(dimdecu(reg))
      real*8 mask(dimdecu(reg))
      real*8 sum
c      do 10 i = 1, (regh0 - regl0 + 1) * (regh1 - regl1 + 1)
      idiff = regh0 - regl0 + 1
      do 10 i = idiff + 2, idiff * (regh1 - regl1) - 1
 10      sum = sum + mask(i) * v0(i) * v1(i)
      end

      subroutine hgcg1(r, p, z, x, w, c, mask,
     @ dims(reg),
     @ alpha, rho)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(reg)
      real*8 r(dimdecu(reg))
      real*8 p(dimdecu(reg))
      real*8 z(dimdecu(reg))
      real*8 x(dimdecu(reg))
      real*8 w(dimdecu(reg))
      real*8 c(dimdecu(reg))
      real*8 mask(dimdecu(reg))
      real*8 alpha, rho
c      do 10 i = 1, (regh0 - regl0 + 1) * (regh1 - regl1 + 1)
      idiff = regh0 - regl0 + 1
      do 10 i = idiff + 2, idiff * (regh1 - regl1) - 1
         r(i) = r(i) - alpha * w(i)
         x(i) = x(i) + alpha * p(i)
         z(i) = r(i) * c(i)
 10      rho = rho + mask(i) * z(i) * r(i)
      end

      subroutine hgcg2(p, z,
     @ dims(reg),
     @ alpha)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(reg)
      real*8 p(dimdecu(reg))
      real*8 z(dimdecu(reg))
      real*8 alpha
c      do 10 i = 1, (regh0 - regl0 + 1) * (regh1 - regl1 + 1)
      idiff = regh0 - regl0 + 1
      do 10 i = idiff + 2, idiff * (regh1 - regl1) - 1
 10      p(i) = alpha * p(i) + z(i)
      end

#endif
#if (CGOPT == 2)

      subroutine hgcg(r, p, z, x, w, c, mask,
     @ ngrids, strid, nvals,
     @ rs, ps, zs, xs, ws, cs, ms,
     @ nsets, edge, nevals, dstart, sstart, dstrid, sstrid,
     @ hx, alpha, rho, it, pcode)
      implicit real*8(a-h,o-z), integer(i-n)
      real*8 r(0:1), p(0:1), z(0:1)
      real*8 x(0:1), w(0:1), c(0:1), mask(0:1)
      integer ngrids
      integer strid(ngrids), nvals(ngrids)
      integer rs(ngrids), ps(ngrids), zs(ngrids)
      integer xs(ngrids), ws(ngrids), cs(ngrids), ms(ngrids)
      real*8 edge(0:1)
      integer nsets
      integer nevals(nsets), dstart(nsets), sstart(nsets)
      integer dstrid(nsets), sstrid(nsets)
      real*8 hx, alpha, rho
      integer it, pcode
      hxm2 = 1.d0 / (hx*hx)
#ifdef CROSS_STENCIL
      facm = -hxm2
#else
      facm = -hxm2 / 3.d0
#endif
      rho = 0.d0
      do 10 k = 1, ngrids
c must do this loop over entire array to initialize borders of p
         do 10 i = 0, nvals(k) - 1
            z(zs(k)+i) = r(rs(k)+i) * c(cs(k)+i)
            p(ps(k)+i) = z(zs(k)+i)
 10         rho = rho + mask(ms(k)+i) * z(zs(k)+i) * r(rs(k)+i)
      tol = 1.d-3 * rho
      if (tol .le. 0.d0) then
         return
      endif
 100  continue
      it = it + 1
      if (it .gt. 250) then
         print *, "FORTRAN hgcg---conjugate-gradient iteration failed"
      endif
      rhoold = rho
      alpha = 0.d0
      do 110 k = 1, nsets
         do 110 i = 0, nevals(k) - 1
 110        edge(dstart(k) + i * dstrid(k)) =
     @      edge(sstart(k) + i * sstrid(k))
      do 120 k = 1, ngrids
         do 120 i = strid(k) + 1, nvals(k) - strid(k) - 2
#ifdef CROSS_STENCIL
            w(ws(k)+i) = facm *
     @        (p(ps(k)+i-1) + p(ps(k)+i+1) +
     @         p(ps(k)+i-strid(k)) + p(ps(k)+i+strid(k)) -
     @         4.d0 * p(ps(k)+i))
#else
            w(ws(k)+i) = facm *
     @        (p(ps(k)+i-1-strid(k)) + p(ps(k)+i-1) +
     @         p(ps(k)+i-1+strid(k)) +
     @         p(ps(k)+i-strid(k)) + p(ps(k)+i+strid(k)) +
     @         p(ps(k)+i+1-strid(k)) + p(ps(k)+i+1) +
     @         p(ps(k)+i+1+strid(k)) -
     @         8.d0 * p(ps(k)+i))
#endif
 120        alpha = alpha + mask(ms(k)+i) * p(ps(k)+i) * w(ws(k)+i)
      alpha = rho / alpha
      rho = 0.d0
      do 130 k = 1, ngrids
         do 130 i = strid(k) + 1, nvals(k) - strid(k) - 2
            r(rs(k)+i) = r(rs(k)+i) - alpha * w(ws(k)+i)
            x(xs(k)+i) = x(xs(k)+i) + alpha * p(ps(k)+i)
            z(zs(k)+i) = r(rs(k)+i) * c(cs(k)+i)
 130        rho = rho + mask(ms(k)+i) * z(zs(k)+i) * r(rs(k)+i)
      if (pcode .ge. 3) then
         print *, it, rho
      endif
      if (rho .le. tol .or. it .gt. 250) then
         return
      endif
      alpha = rho / rhoold
      do 140 k = 1, ngrids
         do 140 i = strid(k) + 1, nvals(k) - strid(k) - 2
 140        p(ps(k)+i) = alpha * p(ps(k)+i) + z(zs(k)+i)
      goto 100
      end

#endif
