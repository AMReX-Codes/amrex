c variable density versions:

c Note---assumes fdst linearly interpolated from cdst along edge
      subroutine hgfres(res, resl0, resh0, resl1, resh1,
     & src,  srcl0, srch0, srcl1, srch1,
     & fdst, fdstl0, fdsth0, fdstl1, fdsth1,
     & cdst, cdstl0, cdsth0, cdstl1, cdsth1,
     & sigmaf, sfl0, sfh0, sfl1, sfh1,
     & sigmac, scl0, sch0, scl1, sch1,
     & regl0, regh0, regl1, regh1,
     & hx, hy, ir, jr, idim, idir)
      integer resl0, resh0, resl1, resh1
      integer srcl0, srch0, srcl1, srch1
      integer fdstl0, fdsth0, fdstl1, fdsth1
      integer cdstl0, cdsth0, cdstl1, cdsth1
      integer sfl0, sfh0, sfl1, sfh1
      integer scl0, sch0, scl1, sch1
      integer regl0, regh0, regl1, regh1
      double precision res(resl0:resh0,resl1:resh1)
      double precision src(srcl0:srch0,srcl1:srch1)
      double precision fdst(fdstl0:fdsth0,fdstl1:fdsth1)
      double precision cdst(cdstl0:cdsth0,cdstl1:cdsth1)
      double precision sigmaf(sfl0:sfh0,sfl1:sfh1)
      double precision sigmac(scl0:sch0,scl1:sch1)
      double precision hx, hy
      integer ir, jr, idim, idir
      double precision hxm2, hym2, fac0, fac1, tmp
      integer i, is, j, js, m, n
      if (idim .eq. 0) then
         i = regl0
         if (idir .eq. 1) then
            is = i - 1
         else
            is = i
         end if
         fac0 = ir / (ir + 1.d0)
         hxm2 = 1.d0 / (ir * ir * hx * hx)
         hym2 = 1.d0 / (jr * jr * hy * hy)
         do 10 j = regl1, regh1
 10         res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     &        (hxm2 *
     &          ((sigmac(is,j-1) + sigmac(is,j)) *
     &            (cdst(i-idir,j) - cdst(i,j))) +
     &         hym2 *
     &          (sigmac(is,j-1) *
     &            (cdst(i,j-1) - cdst(i,j)) +
     &           sigmac(is,j) *
     &            (cdst(i,j+1) - cdst(i,j))))
         fac0 = fac0 / (ir * jr * jr)
         hxm2 = ir * ir * hxm2
         hym2 = jr * jr * hym2
         i = i * ir
         if (idir .eq. 1) then
            is = i
         else
            is = i - 1
         end if
         do 20 n = 0, jr-1
            fac1 = (jr-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do 20 j = jr*regl1, jr*regh1, jr
               tmp = hxm2 *
     &           ((sigmaf(is,j-n-1) + sigmaf(is,j-n)) *
     &             (fdst(i+idir,j-n) - fdst(i,j-n)) +
     &		  (sigmaf(is,j+n-1) + sigmaf(is,j+n)) *
     &             (fdst(i+idir,j+n) - fdst(i,j+n)))
 20            res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     &		 (sigmaf(is,j-n-1) *
     &             (fdst(i,j-n-1) - fdst(i,j-n)) +
     &		  sigmaf(is,j-n) *
     &             (fdst(i,j-n+1) - fdst(i,j-n)) +
     &		  sigmaf(is,j+n-1) *
     &             (fdst(i,j+n-1) - fdst(i,j+n)) +
     &		  sigmaf(is,j+n) *
     &             (fdst(i,j+n+1) - fdst(i,j+n))))
      else
         j = regl1
         if (idir .eq. 1) then
            js = j - 1
         else
            js = j
         end if
         fac0 = jr / (jr + 1.d0)
         hxm2 = 1.d0 / (ir * ir * hx * hx)
         hym2 = 1.d0 / (jr * jr * hy * hy)
         do 30 i = regl0, regh0
 30         res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     &        (hxm2 *
     &          (sigmac(i-1,js) *
     &            (cdst(i-1,j) - cdst(i,j)) +
     &           sigmac(i,js) *
     &            (cdst(i+1,j) - cdst(i,j))) +
     &         hym2 *
     &          ((sigmac(i-1,js) + sigmac(i,js)) *
     &            (cdst(i,j-idir) - cdst(i,j))))
         fac0 = fac0 / (ir * ir * jr)
         hxm2 = ir * ir * hxm2
         hym2 = jr * jr * hym2
         j = j * jr
         if (idir .eq. 1) then
            js = j
         else
            js = j - 1
         end if
         do 40 m = 0, ir-1
            fac1 = (ir-m) * fac0
            if (m .eq. 0) fac1 = 0.5d0 * fac1
            do 40 i = ir*regl0, ir*regh0, ir
               tmp = hxm2 *
     &		 (sigmaf(i-m-1,js) *
     &             (fdst(i-m-1,j) - fdst(i-m,j)) +
     &		  sigmaf(i-m,js) *
     &             (fdst(i-m+1,j) - fdst(i-m,j)) +
     &		  sigmaf(i+m-1,js) *
     &             (fdst(i+m-1,j) - fdst(i+m,j)) +
     &		  sigmaf(i+m,js) *
     &             (fdst(i+m+1,j) - fdst(i+m,j)))
 40            res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     &           ((sigmaf(i-m-1,js) + sigmaf(i-m,js)) *
     &             (fdst(i-m,j+idir) - fdst(i-m,j)) +
     &		  (sigmaf(i+m-1,js) + sigmaf(i+m,js)) *
     &             (fdst(i+m,j+idir) - fdst(i+m,j))))
      end if
      end

c Note---assumes fdst linearly interpolated from cdst along edges
      subroutine hgcres(res, resl0, resh0, resl1, resh1,
     & src,  srcl0, srch0, srcl1, srch1,
     & fdst, fdstl0, fdsth0, fdstl1, fdsth1,
     & cdst, cdstl0, cdsth0, cdstl1, cdsth1,
     & sigmaf, sfl0, sfh0, sfl1, sfh1,
     & sigmac, scl0, sch0, scl1, sch1,
     & regl0, regh0, regl1, regh1,
     & hx, hy, ir, jr, ga)
      integer resl0, resh0, resl1, resh1
      integer srcl0, srch0, srcl1, srch1
      integer fdstl0, fdsth0, fdstl1, fdsth1
      integer cdstl0, cdsth0, cdstl1, cdsth1
      integer sfl0, sfh0, sfl1, sfh1
      integer scl0, sch0, scl1, sch1
      integer regl0, regh0, regl1, regh1
      double precision res(resl0:resh0,resl1:resh1)
      double precision src(srcl0:srch0,srcl1:srch1)
      double precision fdst(fdstl0:fdsth0,fdstl1:fdsth1)
      double precision cdst(cdstl0:cdsth0,cdstl1:cdsth1)
      double precision sigmaf(sfl0:sfh0,sfl1:sfh1)
      double precision sigmac(scl0:sch0,scl1:sch1)
      double precision hx, hy
      integer ir, jr, ga(0:1,0:1)
      double precision hxm2, hym2, hxm2c, hym2c, sum, center, 
     &	ffac, cfac, fac, fac1
      integer ic, jc, if, jf, ii, ji, idir, jdir, m, n
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
     &           (hxm2 * (fdst(if+idir,jf) - fdst(if,jf)) +
     &            hym2 * (fdst(if,jf+jdir) - fdst(if,jf)))
            else
               center = center + cfac
               sum = sum + ir * jr * sigmac(ic+ii-1,jc+ji-1) *
     &           (hxm2c * (cdst(ic+idir,jc) - cdst(ic,jc)) +
     &            hym2c * (cdst(ic,jc+jdir) - cdst(ic,jc)))
            end if
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
     &              (hxm2 * (sigmaf(if+m-1,jf+ji-1) *
     &                        (fdst(if+m-1,jf) - fdst(if+m,jf)) +
     &                       sigmaf(if+m,jf+ji-1) *
     &                        (fdst(if+m+1,jf) - fdst(if+m,jf))) +
     &               hym2 *
     &                 (sigmaf(if+m-1,jf+ji-1) + sigmaf(if+m,jf+ji-1)) *
     &                 (fdst(if+m,jf+jdir) - fdst(if+m,jf)))
            end if
            if (ga(ii,ji) - ga(1-ii,ji) .eq. 1) then
               fac1 = 1.d0 / jr
               ffac = 0.5d0 * (jr-1)
               center = center + ffac
               do 30 n = jdir, jdir*(jr-1), jdir
                  fac = (jr-abs(n)) * fac1
 30               sum = sum + fac *
     &              (hxm2 *
     &                 (sigmaf(if+ii-1,jf+n-1) + sigmaf(if+ii-1,jf+n)) *
     &                 (fdst(if+idir,jf+n) - fdst(if,jf+n)) +
     &               hym2 * (sigmaf(if+ii-1,jf+n-1) *
     &                        (fdst(if,jf+n-1) - fdst(if,jf+n)) +
     &                       sigmaf(if+ii-1,jf+n) *
     &                        (fdst(if,jf+n+1) - fdst(if,jf+n))))
            end if
 40         continue
c weighting
      res(if,jf) = src(if,jf) - sum / center
      end


c 9-point stencil versions:


c 9-point variable density stencils:
c-----------------------------------------------------------------------
c Note---assumes fdst linearly interpolated from cdst along edge
      subroutine hgfres_full_stencil(res, resl0, resh0, resl1, resh1,
     & src,  srcl0, srch0, srcl1, srch1,
     & fdst, fdstl0, fdsth0, fdstl1, fdsth1,
     & cdst, cdstl0, cdsth0, cdstl1, cdsth1,
     & sigmaf, sfl0, sfh0, sfl1, sfh1,
     & sigmac, scl0, sch0, scl1, sch1,
     & regl0, regh0, regl1, regh1,
     & hx, hy, ir, jr, idim, idir, irz, imax)
      integer resl0, resh0, resl1, resh1
      integer srcl0, srch0, srcl1, srch1
      integer fdstl0, fdsth0, fdstl1, fdsth1
      integer cdstl0, cdsth0, cdstl1, cdsth1
      integer sfl0, sfh0, sfl1, sfh1
      integer scl0, sch0, scl1, sch1
      integer regl0, regh0, regl1, regh1
      double precision res(resl0:resh0,resl1:resh1)
      double precision src(srcl0:srch0,srcl1:srch1)
      double precision fdst(fdstl0:fdsth0,fdstl1:fdsth1)
      double precision cdst(cdstl0:cdsth0,cdstl1:cdsth1)
      double precision sigmaf(sfl0:sfh0,sfl1:sfh1)
      double precision sigmac(scl0:sch0,scl1:sch1)
      double precision hx, hy
      integer ir, jr, idim, idir, irz, imax
      double precision hxm2, hym2, fac0, fac1, r, rfac,
     &		rfac0, rfac1, tmp
      double precision rfac0m, rfac1m, rfac0p, rfac1p
      integer i, j, is, js, m, n
      if (idim .eq. 0) then
         i = regl0
         if (idir .eq. 1) then
            is = i - 1
         else
            is = i
         end if
         fac0 = ir / (3.d0 * (ir + 1.d0))
         hxm2 = 1.d0 / (ir * ir * hx * hx)
         hym2 = 1.d0 / (jr * jr * hy * hy)
         do 10 j = regl1, regh1
 10         res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     &        (hxm2 *
     &          (sigmac(is,j-1) *
     &            (cdst(i-idir,j-1) - cdst(i,j-1) +
     &             2.d0 * (cdst(i-idir,j) - cdst(i,j))) +
     &           sigmac(is,j) *
     &            (cdst(i-idir,j+1) - cdst(i,j+1) +
     &             2.d0 * (cdst(i-idir,j) - cdst(i,j)))) +
     &         hym2 *
     &          (sigmac(is,j-1) *
     &            (cdst(i-idir,j-1) - cdst(i-idir,j) +
     &             2.d0 * (cdst(i,j-1) - cdst(i,j))) +
     &           sigmac(is,j) *
     &            (cdst(i-idir,j+1) - cdst(i-idir,j) +
     &             2.d0 * (cdst(i,j+1) - cdst(i,j)))))
         if (irz .eq. 1) then
            r = (is + 0.5d0) * (hx * ir)
            rfac = idir * ir * hx * hym2 / (2.d0 * r)
            do 20 j = regl1, regh1
 20            res(i*ir,j*jr) = res(i*ir,j*jr) - fac0 *
     &           (rfac * (sigmac(is,j-1) * (cdst(i,j-1) - cdst(i,j)) +
     &                    sigmac(is,j)   * (cdst(i,j+1) - cdst(i,j))))
         end if
         fac0 = fac0 / (ir * jr * jr)
         hxm2 = ir * ir * hxm2
         hym2 = jr * jr * hym2
         i = i * ir
         if (idir .eq. 1) then
            is = i
         else
            is = i - 1
         end if
         if (irz .eq. 1) then
            r = (is + 0.5d0) * hx
            rfac = idir * hx * hym2 / (2.d0 * r)
         end if
         do 50 n = 0, jr-1
            fac1 = (jr-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do 30 j = jr*regl1, jr*regh1, jr
               tmp = hxm2 *
     &           (sigmaf(is,j-n-1) *
     &             (fdst(i+idir,j-n-1) - fdst(i,j-n-1) +
     &              2.d0 * (fdst(i+idir,j-n) - fdst(i,j-n))) +
     &		  sigmaf(is,j-n) *
     &             (fdst(i+idir,j-n+1) - fdst(i,j-n+1) +
     &              2.d0 * (fdst(i+idir,j-n) - fdst(i,j-n))) +
     &		  sigmaf(is,j+n-1) *
     &             (fdst(i+idir,j+n-1) - fdst(i,j+n-1) +
     &              2.d0 * (fdst(i+idir,j+n) - fdst(i,j+n))) +
     &		  sigmaf(is,j+n) *
     &             (fdst(i+idir,j+n+1) - fdst(i,j+n+1) +
     &              2.d0 * (fdst(i+idir,j+n) - fdst(i,j+n))))
 30            res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     &		 (sigmaf(is,j-n-1) *
     &             (fdst(i+idir,j-n-1) - fdst(i+idir,j-n) +
     &              2.d0 * (fdst(i,j-n-1) - fdst(i,j-n))) +
     &		  sigmaf(is,j-n) *
     &             (fdst(i+idir,j-n+1) - fdst(i+idir,j-n) +
     &              2.d0 * (fdst(i,j-n+1) - fdst(i,j-n))) +
     &		  sigmaf(is,j+n-1) *
     &             (fdst(i+idir,j+n-1) - fdst(i+idir,j+n) +
     &              2.d0 * (fdst(i,j+n-1) - fdst(i,j+n))) +
     &		  sigmaf(is,j+n) *
     &             (fdst(i+idir,j+n+1) - fdst(i+idir,j+n) +
     &              2.d0 * (fdst(i,j+n+1) - fdst(i,j+n)))))
            if (irz .eq. 1) then
               do 40 j = jr*regl1, jr*regh1, jr
 40               res(i,j) = res(i,j) + fac1 *
     &       (rfac * (sigmaf(is,j-n-1) * (fdst(i,j-n-1) - fdst(i,j-n)) +
     &                sigmaf(is,j-n)   * (fdst(i,j-n+1) - fdst(i,j-n)) +
     &                sigmaf(is,j+n-1) * (fdst(i,j+n-1) - fdst(i,j+n)) +
     &                sigmaf(is,j+n)   * (fdst(i,j+n+1) - fdst(i,j+n))))
            end if
 50         continue
      else
         j = regl1
         if (idir .eq. 1) then
            js = j - 1
         else
            js = j
         end if
         fac0 = jr / (3.d0 * (jr + 1.d0))
         hxm2 = 1.d0 / (ir * ir * hx * hx)
         hym2 = 1.d0 / (jr * jr * hy * hy)
         do 60 i = regl0, regh0
 60         res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     &        (hxm2 *
     &          (sigmac(i-1,js) *
     &            (cdst(i-1,j-idir) - cdst(i,j-idir) +
     &             2.d0 * (cdst(i-1,j) - cdst(i,j))) +
     &           sigmac(i,js) *
     &            (cdst(i+1,j-idir) - cdst(i,j-idir) +
     &             2.d0 * (cdst(i+1,j) - cdst(i,j)))) +
     &         hym2 *
     &          (sigmac(i-1,js) *
     &            (cdst(i-1,j-idir) - cdst(i-1,j) +
     &             2.d0 * (cdst(i,j-idir) - cdst(i,j))) +
     &           sigmac(i,js) *
     &            (cdst(i+1,j-idir) - cdst(i+1,j) +
     &             2.d0 * (cdst(i,j-idir) - cdst(i,j)))))
         if (irz .eq. 1 .and. regh0 .lt. imax) then
            do 70 i = regl0, regh0
               r = (i + 0.5d0) * (hx * ir)
               rfac0 = ir * hx * hym2 / (2.d0 * (r - hx * ir))
               rfac1 = ir * hx * hym2 / (2.d0 * r)
 70            res(i*ir,j*jr) = res(i*ir,j*jr) - fac0 *
     &           (rfac0 * sigmac(i-1,js) - rfac1 * sigmac(i,js)) *
     &           (cdst(i,j-idir) - cdst(i,j))
         else if (irz .eq. 1) then
c This should only occur with a corner at the outer boundary:
            i = regh0
            r = (i - 0.5d0) * (hx * ir)
            rfac0 = ir * hx * hym2 / (2.d0 * r)
            rfac1 = -rfac0
            res(i*ir,j*jr) = res(i*ir,j*jr) - fac0 *
     &           (rfac0 * sigmac(i-1,js) - rfac1 * sigmac(i,js)) *
     &           (cdst(i,j-idir) - cdst(i,j))
         end if
         fac0 = fac0 / (ir * ir * jr)
         hxm2 = ir * ir * hxm2
         hym2 = jr * jr * hym2
         j = j * jr
         if (idir .eq. 1) then
            js = j
         else
            js = j - 1
         end if
         do 100 m = 0, ir-1
            fac1 = (ir-m) * fac0
            if (m .eq. 0) fac1 = 0.5d0 * fac1
            do 80 i = ir*regl0, ir*regh0, ir
               tmp = hxm2 *
     &		 (sigmaf(i-m-1,js) *
     &             (fdst(i-m-1,j+idir) - fdst(i-m,j+idir) +
     &              2.d0 * (fdst(i-m-1,j) - fdst(i-m,j))) +
     &		  sigmaf(i-m,js) *
     &             (fdst(i-m+1,j+idir) - fdst(i-m,j+idir) +
     &              2.d0 * (fdst(i-m+1,j) - fdst(i-m,j))) +
     &		  sigmaf(i+m-1,js) *
     &             (fdst(i+m-1,j+idir) - fdst(i+m,j+idir) +
     &              2.d0 * (fdst(i+m-1,j) - fdst(i+m,j))) +
     &		  sigmaf(i+m,js) *
     &             (fdst(i+m+1,j+idir) - fdst(i+m,j+idir) +
     &              2.d0 * (fdst(i+m+1,j) - fdst(i+m,j))))
 80            res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     &           (sigmaf(i-m-1,js) *
     &             (fdst(i-m-1,j+idir) - fdst(i-m-1,j) +
     &              2.d0 * (fdst(i-m,j+idir) - fdst(i-m,j))) +
     &		  sigmaf(i-m,js) *
     &             (fdst(i-m+1,j+idir) - fdst(i-m+1,j) +
     &              2.d0 * (fdst(i-m,j+idir) - fdst(i-m,j))) +
     &		  sigmaf(i+m-1,js) *
     &             (fdst(i+m-1,j+idir) - fdst(i+m-1,j) +
     &              2.d0 * (fdst(i+m,j+idir) - fdst(i+m,j))) +
     &		  sigmaf(i+m,js) *
     &             (fdst(i+m+1,j+idir) - fdst(i+m+1,j) +
     &              2.d0 * (fdst(i+m,j+idir) - fdst(i+m,j)))))
         if (irz .eq. 1 .and. regh0 .lt. imax) then
            do 90 i = ir*regl0, ir*regh0, ir
               r = (i + 0.5d0) * hx
               rfac0m = hx * hym2 / (2.d0 * (r - (m + 1) * hx))
               rfac1m = hx * hym2 / (2.d0 * (r - m * hx))
               rfac0p = hx * hym2 / (2.d0 * (r + (m - 1) * hx))
               rfac1p = hx * hym2 / (2.d0 * (r + m * hx))
 90            res(i,j) = res(i,j) - fac1 *
     &          ((rfac0m * sigmaf(i-m-1,js) - rfac1m * sigmaf(i-m,js)) *
     &           (fdst(i-m,j+idir) - fdst(i-m,j)) +
     &           (rfac0p * sigmaf(i+m-1,js) - rfac1p * sigmaf(i+m,js)) *
     &           (fdst(i+m,j+idir) - fdst(i+m,j)))
         else if (irz .eq. 1) then
c This should only occur with a corner at the outer boundary:
            i = ir * regh0
            r = (i + 0.5d0) * hx
            rfac0m = hx * hym2 / (2.d0 * (r - (m + 1) * hx))
            if (m .eq. 0) then
               rfac1m = -rfac0m
            else
               rfac1m = hx * hym2 / (2.d0 * (r - m * hx))
            end if
            rfac0p = -rfac1m
            rfac1p = -rfac0m
            res(i,j) = res(i,j) - fac1 *
     &          ((rfac0m * sigmaf(i-m-1,js) - rfac1m * sigmaf(i-m,js)) *
     &           (fdst(i-m,j+idir) - fdst(i-m,j)) +
     &           (rfac0p * sigmaf(i+m-1,js) - rfac1p * sigmaf(i+m,js)) *
     &           (fdst(i+m,j+idir) - fdst(i+m,j)))
         end if
 100     continue
      end if
      end
c-----------------------------------------------------------------------
c Note---assumes fdst linearly interpolated from cdst along edge
      subroutine hgores_full_stencil(res, resl0, resh0, resl1, resh1,
     & src,  srcl0, srch0, srcl1, srch1,
     & fdst, fdstl0, fdsth0, fdstl1, fdsth1,
     & cdst, cdstl0, cdsth0, cdstl1, cdsth1,
     & sigmaf, sfl0, sfh0, sfl1, sfh1,
     & sigmac, scl0, sch0, scl1, sch1,
     & regl0, regh0, regl1, regh1,
     & hx, hy, ir, jr, idir, jdir, irz)
      integer resl0, resh0, resl1, resh1
      integer srcl0, srch0, srcl1, srch1
      integer fdstl0, fdsth0, fdstl1, fdsth1
      integer cdstl0, cdsth0, cdstl1, cdsth1
      integer sfl0, sfh0, sfl1, sfh1
      integer scl0, sch0, scl1, sch1
      integer regl0, regh0, regl1, regh1
      double precision res(resl0:resh0,resl1:resh1)
      double precision src(srcl0:srch0,srcl1:srch1)
      double precision fdst(fdstl0:fdsth0,fdstl1:fdsth1)
      double precision cdst(cdstl0:cdsth0,cdstl1:cdsth1)
      double precision sigmaf(sfl0:sfh0,sfl1:sfh1)
      double precision sigmac(scl0:sch0,scl1:sch1)
      double precision hx, hy
      integer ir, jr, idir, jdir, irz
      double precision hxm2, hym2, fac0, fac1, r, rfac,
     &		rfac0, rfac1, tmp
      double precision rfac0p, rfac1p
      integer i, j, is, js, m, n
      i = regl0
      j = regl1
      if (idir .eq. 1) then
         is = i - 1
      else
         is = i
      end if
      if (jdir .eq. 1) then
         js = j - 1
      else
         js = j
      end if
      hxm2 = 1.d0 / (ir * ir * hx * hx)
      hym2 = 1.d0 / (jr * jr * hy * hy)
      fac0 = (ir * jr) / (4.5d0*ir*jr + 1.5d0*ir + 1.5d0*jr - 1.5d0)
      res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     &  (hxm2 * (sigmac(is,js+jdir) *
     &            (cdst(i-idir,j+jdir) - cdst(i,j+jdir) +
     &             2.d0 * (cdst(i-idir,j) - cdst(i,j))) +
     &           sigmac(is,js) *
     &            (cdst(i-idir,j-jdir) - cdst(i,j-jdir) +
     &             2.d0 * (cdst(i-idir,j) - cdst(i,j))) +
     &           sigmac(is+idir,js) *
     &            (cdst(i+idir,j-jdir) - cdst(i,j-jdir) +
     &             2.d0 * (cdst(i+idir,j) - cdst(i,j)))) +
     &   hym2 * (sigmac(is,js+jdir) *
     &            (cdst(i-idir,j+jdir) - cdst(i-idir,j) +
     &             2.d0 * (cdst(i,j+jdir) - cdst(i,j))) +
     &           sigmac(is,js) *
     &            (cdst(i-idir,j-jdir) - cdst(i-idir,j) +
     &             2.d0 * (cdst(i,j-jdir) - cdst(i,j))) +
     &           sigmac(is+idir,js) *
     &            (cdst(i+idir,j-jdir) - cdst(i+idir,j) +
     &             2.d0 * (cdst(i,j-jdir) - cdst(i,j)))))
      if (irz .eq. 1) then
         r = (is + 0.5d0) * (hx * ir)
         rfac0 = ir * hx * hym2 / (2.d0 * (r + idir * hx * ir))
         rfac1 = ir * hx * hym2 / (2.d0 * r)
         res(i*ir,j*jr) = res(i*ir,j*jr) + idir * fac0 *
     &     ((rfac0 * sigmac(is+idir,js) - rfac1 * sigmac(is,js)) *
     &      (cdst(i,j-jdir) - cdst(i,j)) -
     &      rfac1 * sigmac(is,js+jdir) * (cdst(i,j+jdir) - cdst(i,j)))
      end if
      fac0 = fac0 / (ir * jr)
      hxm2 = ir * ir * hxm2
      hym2 = jr * jr * hym2
      i = i * ir
      j = j * jr
      if (idir .eq. 1) then
         is = i
      else
         is = i - 1
      end if
      if (jdir .eq. 1) then
         js = j
      else
         js = j - 1
      end if
      res(i,j) = res(i,j) - fac0 * sigmaf(is,js) *
     &  (hxm2 * (fdst(i+idir,j+jdir) - fdst(i,j+jdir) +
     &           2.d0 * (fdst(i+idir,j) - fdst(i,j))) +
     &   hym2 * (fdst(i+idir,j+jdir) - fdst(i+idir,j) +
     &           2.d0 * (fdst(i,j+jdir) - fdst(i,j))))
      if (irz .eq. 1) then
         r = (is + 0.5d0) * hx
         rfac = hx * hym2 / (2.d0 * r)
         res(i,j) = res(i,j) + idir * fac0 *
     &      rfac * sigmaf(is,js) * (fdst(i,j+jdir) - fdst(i,j))
      end if
      fac0 = fac0 / ir
      do 10 m = idir, idir*(ir-1), idir
         fac1 = (ir-abs(m)) * fac0
         tmp = hxm2 *
     &     (sigmaf(i+m-1,js) *
     &       (fdst(i+m-1,j+jdir) - fdst(i+m,j+jdir) +
     &        2.d0 * (fdst(i+m-1,j) - fdst(i+m,j))) +
     &      sigmaf(i+m,js) *
     &       (fdst(i+m+1,j+jdir) - fdst(i+m,j+jdir) +
     &        2.d0 * (fdst(i+m+1,j) - fdst(i+m,j))))
         res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     &     (sigmaf(i+m-1,js) *
     &       (fdst(i+m-1,j+jdir) - fdst(i+m-1,j) +
     &        2.d0 * (fdst(i+m,j+jdir) - fdst(i+m,j))) +
     &	    sigmaf(i+m,js) *
     &       (fdst(i+m+1,j+jdir) - fdst(i+m+1,j) +
     &        2.d0 * (fdst(i+m,j+jdir) - fdst(i+m,j)))))
         if (irz .eq. 1) then
            r = (i + m + 0.5d0) * hx
            rfac0p = hx * hym2 / (2.d0 * (r - hx))
            rfac1p = hx * hym2 / (2.d0 * r)
            res(i,j) = res(i,j) - fac1 *
     &         (rfac0p * sigmaf(i+m-1,js) - rfac1p * sigmaf(i+m,js)) *
     &          (fdst(i+m,j+jdir) - fdst(i+m,j))
         end if
 10      continue
      fac0 = ir * fac0 / jr
      do 20 n = jdir, jdir*(jr-1), jdir
         fac1 = (jr-abs(n)) * fac0
         tmp = hxm2 *
     &     (sigmaf(is,j+n-1) *
     &       (fdst(i+idir,j+n-1) - fdst(i,j+n-1) +
     &        2.d0 * (fdst(i+idir,j+n) - fdst(i,j+n))) +
     &      sigmaf(is,j+n) *
     &       (fdst(i+idir,j+n+1) - fdst(i,j+n+1) +
     &        2.d0 * (fdst(i+idir,j+n) - fdst(i,j+n))))
         res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     &     (sigmaf(is,j+n-1) *
     &       (fdst(i+idir,j+n-1) - fdst(i+idir,j+n) + 
     &        2.d0 * (fdst(i,j+n-1) - fdst(i,j+n))) +
     &	    sigmaf(is,j+n) *
     &       (fdst(i+idir,j+n+1) - fdst(i+idir,j+n) +
     &        2.d0 * (fdst(i,j+n+1) - fdst(i,j+n)))))
         if (irz .eq. 1) then
            r = (is + 0.5d0) * hx
            rfac1 = -idir * hx * hym2 / (2.d0 * r)
            res(i,j) = res(i,j) - fac1 *
     &        rfac1 *(sigmaf(is,j+n-1) * (fdst(i,j+n-1) - fdst(i,j+n)) +
     &                sigmaf(is,j+n)   * (fdst(i,j+n+1) - fdst(i,j+n)))
         end if
 20      continue
      end
c-----------------------------------------------------------------------
c Note---assumes fdst linearly interpolated from cdst along edge
      subroutine hgires_full_stencil(res, resl0, resh0, resl1, resh1,
     & src,  srcl0, srch0, srcl1, srch1,
     & fdst, fdstl0, fdsth0, fdstl1, fdsth1,
     & cdst, cdstl0, cdsth0, cdstl1, cdsth1,
     & sigmaf, sfl0, sfh0, sfl1, sfh1,
     & sigmac, scl0, sch0, scl1, sch1,
     & regl0, regh0, regl1, regh1,
     & hx, hy, ir, jr, idir, jdir, irz)
      integer resl0, resh0, resl1, resh1
      integer srcl0, srch0, srcl1, srch1
      integer fdstl0, fdsth0, fdstl1, fdsth1
      integer cdstl0, cdsth0, cdstl1, cdsth1
      integer sfl0, sfh0, sfl1, sfh1
      integer scl0, sch0, scl1, sch1
      integer regl0, regh0, regl1, regh1
      double precision res(resl0:resh0,resl1:resh1)
      double precision src(srcl0:srch0,srcl1:srch1)
      double precision fdst(fdstl0:fdsth0,fdstl1:fdsth1)
      double precision cdst(cdstl0:cdsth0,cdstl1:cdsth1)
      double precision sigmaf(sfl0:sfh0,sfl1:sfh1)
      double precision sigmac(scl0:sch0,scl1:sch1)
      double precision hx, hy
      integer ir, jr, idir, jdir, irz
      double precision hxm2, hym2, fac0, fac1, r, rfac,
     &		rfac0, rfac1, tmp
      double precision rfac0m, rfac1m
      integer i, j, is, js, m, n
      i = regl0
      j = regl1
      if (idir .eq. 1) then
         is = i - 1
      else
         is = i
      end if
      if (jdir .eq. 1) then
         js = j - 1
      else
         js = j
      end if
      hxm2 = 1.d0 / (ir * ir * hx * hx)
      hym2 = 1.d0 / (jr * jr * hy * hy)
      fac0 = (ir * jr) / (1.5d0*ir*jr + 1.5d0*ir + 1.5d0*jr + 1.5d0)
      res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 * sigmac(is,js) *
     &  (hxm2 * (cdst(i-idir,j-jdir) - cdst(i,j-jdir) +
     &           2.d0 * (cdst(i-idir,j) - cdst(i,j))) +
     &   hym2 * (cdst(i-idir,j-jdir) - cdst(i-idir,j) +
     &           2.d0 * (cdst(i,j-jdir) - cdst(i,j))))
      if (irz .eq. 1) then
         r = (is + 0.5d0) * (hx * ir)
         rfac = ir * hx * hym2 / (2.d0 * r)
         res(i*ir,j*jr) = res(i*ir,j*jr) - idir * fac0 *
     &     rfac * sigmac(is,js) * (cdst(i,j-jdir) - cdst(i,j))
      end if
      fac0 = fac0 / (ir * jr)
      hxm2 = ir * ir * hxm2
      hym2 = jr * jr * hym2
      i = i * ir
      j = j * jr
      if (idir .eq. 1) then
         is = i
      else
         is = i - 1
      end if
      if (jdir .eq. 1) then
         js = j
      else
         js = j - 1
      end if
      res(i,j) = res(i,j) - fac0 *
     &  (hxm2 * (sigmaf(is,js-jdir) *
     &            (fdst(i+idir,j-jdir) - fdst(i,j-jdir) +
     &             2.d0 * (fdst(i+idir,j) - fdst(i,j))) +
     &           sigmaf(is,js) *
     &            (fdst(i+idir,j+jdir) - fdst(i,j+jdir) +
     &             2.d0 * (fdst(i+idir,j) - fdst(i,j))) +
     &           sigmaf(is-idir,js) *
     &            (fdst(i-idir,j+jdir) - fdst(i,j+jdir) +
     &             2.d0 * (fdst(i-idir,j) - fdst(i,j)))) +
     &   hym2 * (sigmaf(is,js-jdir) *
     &            (fdst(i+idir,j-jdir) - fdst(i+idir,j) +
     &             2.d0 * (fdst(i,j-jdir) - fdst(i,j))) +
     &           sigmaf(is,js) *
     &            (fdst(i+idir,j+jdir) - fdst(i+idir,j) +
     &             2.d0 * (fdst(i,j+jdir) - fdst(i,j))) +
     &           sigmaf(is-idir,js) *
     &            (fdst(i-idir,j+jdir) - fdst(i-idir,j) +
     &             2.d0 * (fdst(i,j+jdir) - fdst(i,j)))))
      if (irz .eq. 1) then
         r = (is + 0.5d0) * hx
         rfac0 = hx * hym2 / (2.d0 * (r - idir * hx))
         rfac1 = hx * hym2 / (2.d0 * r)
         res(i,j) = res(i,j) - idir * fac0 *
     &     ((rfac0 * sigmaf(is-idir,js) - rfac1 * sigmaf(is,js)) *
     &      (fdst(i,j+jdir) - fdst(i,j)) -
     &      rfac1 * sigmaf(is,js-jdir) * (fdst(i,j-jdir) - fdst(i,j)))
      end if
      fac0 = fac0 / ir
      do 10 m = idir, idir*(ir-1), idir
         fac1 = (ir-abs(m)) * fac0
         tmp = hxm2 *
     &     (sigmaf(i-m-1,js) *
     &       (fdst(i-m-1,j+jdir) - fdst(i-m,j+jdir) +
     &        2.d0 * (fdst(i-m-1,j) - fdst(i-m,j))) +
     &      sigmaf(i-m,js) *
     &       (fdst(i-m+1,j+jdir) - fdst(i-m,j+jdir) +
     &        2.d0 * (fdst(i-m+1,j) - fdst(i-m,j))))
         res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     &     (sigmaf(i-m-1,js) *
     &       (fdst(i-m-1,j+jdir) - fdst(i-m-1,j) +
     &        2.d0 * (fdst(i-m,j+jdir) - fdst(i-m,j))) +
     &	    sigmaf(i-m,js) *
     &       (fdst(i-m+1,j+jdir) - fdst(i-m+1,j) +
     &        2.d0 * (fdst(i-m,j+jdir) - fdst(i-m,j)))))
         if (irz .eq. 1) then
            r = (i - m + 0.5d0) * hx
            rfac0m = hx * hym2 / (2.d0 * (r - hx))
            rfac1m = hx * hym2 / (2.d0 * r)
            res(i,j) = res(i,j) - fac1 *
     &         (rfac0m * sigmaf(i-m-1,js) - rfac1m * sigmaf(i-m,js)) *
     &          (fdst(i-m,j+jdir) - fdst(i-m,j))
         end if
 10      continue
      fac0 = ir * fac0 / jr
      do 20 n = jdir, jdir*(jr-1), jdir
         fac1 = (jr-abs(n)) * fac0
         tmp = hxm2 *
     &     (sigmaf(is,j-n-1) *
     &       (fdst(i+idir,j-n-1) - fdst(i,j-n-1) +
     &        2.d0 * (fdst(i+idir,j-n) - fdst(i,j-n))) +
     &      sigmaf(is,j-n) *
     &       (fdst(i+idir,j-n+1) - fdst(i,j-n+1) +
     &        2.d0 * (fdst(i+idir,j-n) - fdst(i,j-n))))
         res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     &     (sigmaf(is,j-n-1) *
     &       (fdst(i+idir,j-n-1) - fdst(i+idir,j-n) + 
     &        2.d0 * (fdst(i,j-n-1) - fdst(i,j-n))) +
     &	    sigmaf(is,j-n) *
     &       (fdst(i+idir,j-n+1) - fdst(i+idir,j-n) +
     &        2.d0 * (fdst(i,j-n+1) - fdst(i,j-n)))))
         if (irz .eq. 1) then
            r = (is + 0.5d0) * hx
            rfac1 = -idir * hx * hym2 / (2.d0 * r)
            res(i,j) = res(i,j) - fac1 *
     &        rfac1 *(sigmaf(is,j-n-1) * (fdst(i,j-n-1) - fdst(i,j-n)) +
     &                sigmaf(is,j-n)   * (fdst(i,j-n+1) - fdst(i,j-n)))
         end if
 20      continue
      end
c-----------------------------------------------------------------------
c Note---assumes fdst linearly interpolated from cdst along edge
      subroutine hgdres_full_stencil(res, resl0, resh0, resl1, resh1,
     & src,  srcl0, srch0, srcl1, srch1,
     & fdst, fdstl0, fdsth0, fdstl1, fdsth1,
     & cdst, cdstl0, cdsth0, cdstl1, cdsth1,
     & sigmaf, sfl0, sfh0, sfl1, sfh1,
     & sigmac, scl0, sch0, scl1, sch1,
     & regl0, regh0, regl1, regh1,
     & hx, hy, ir, jr, jdir, irz)
      integer resl0, resh0, resl1, resh1
      integer srcl0, srch0, srcl1, srch1
      integer fdstl0, fdsth0, fdstl1, fdsth1
      integer cdstl0, cdsth0, cdstl1, cdsth1
      integer sfl0, sfh0, sfl1, sfh1
      integer scl0, sch0, scl1, sch1
      integer regl0, regh0, regl1, regh1
      double precision res(resl0:resh0,resl1:resh1)
      double precision src(srcl0:srch0,srcl1:srch1)
      double precision fdst(fdstl0:fdsth0,fdstl1:fdsth1)
      double precision cdst(cdstl0:cdsth0,cdstl1:cdsth1)
      double precision sigmaf(sfl0:sfh0,sfl1:sfh1)
      double precision sigmac(scl0:sch0,scl1:sch1)
      double precision hx, hy
      integer ir, jr, jdir, irz
      double precision hxm2, hym2, fac0, fac1, r, rfac0, rfac1, tmp
      double precision rfac0m, rfac1m, rfac0p, rfac1p
      integer i, j, js, m, n
      i = regl0
      j = regl1
      if (jdir .eq. 1) then
         js = j
      else
         js = j - 1
      end if
      hxm2 = 1.d0 / (ir * ir * hx * hx)
      hym2 = 1.d0 / (jr * jr * hy * hy)
      fac0 = (ir * jr) / (3.d0*ir*jr + 3.d0*ir + 3.d0*jr - 3.d0)
      res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     &  (hxm2 * (sigmac(i,js-jdir) *
     &            (cdst(i+1,j-jdir) - cdst(i,j-jdir) +
     &             2.d0 * (cdst(i+1,j) - cdst(i,j))) +
     &           sigmac(i-1,js) *
     &            (cdst(i-1,j+jdir) - cdst(i,j+jdir) +
     &             2.d0 * (cdst(i-1,j) - cdst(i,j)))) +
     &   hym2 * (sigmac(i,js-jdir) *
     &            (cdst(i+1,j-jdir) - cdst(i+1,j) +
     &             2.d0 * (cdst(i,j-jdir) - cdst(i,j))) +
     &           sigmac(i-1,js) *
     &            (cdst(i-1,j+jdir) - cdst(i-1,j) +
     &             2.d0 * (cdst(i,j+jdir) - cdst(i,j)))))
      if (irz .eq. 1) then
         r = (i + 0.5d0) * (hx * ir)
         rfac0 = ir * hx * hym2 / (2.d0 * (r - hx * ir))
         rfac1 = ir * hx * hym2 / (2.d0 * r)
         res(i*ir,j*jr) = res(i*ir,j*jr) - fac0 *
     &     (rfac0 * sigmac(i-1,js)    * (cdst(i,j+jdir) - cdst(i,j)) -
     &      rfac1 * sigmac(i,js-jdir) * (cdst(i,j-jdir) - cdst(i,j)))
      end if
      fac0 = fac0 / (ir * jr)
      hxm2 = ir * ir * hxm2
      hym2 = jr * jr * hym2
      i = i * ir
      j = j * jr
      if (jdir .eq. 1) then
         js = j
      else
         js = j - 1
      end if
      res(i,j) = res(i,j) - fac0 *
     &  (hxm2 * (sigmaf(i,js) *
     &            (fdst(i+1,j+jdir) - fdst(i,j+jdir) +
     &             2.d0 * (fdst(i+1,j) - fdst(i,j))) +
     &           sigmaf(i-1,js-jdir) *
     &            (fdst(i-1,j-jdir) - fdst(i,j-jdir) +
     &             2.d0 * (fdst(i-1,j) - fdst(i,j)))) +
     &   hym2 * (sigmaf(i,js) *
     &            (fdst(i+1,j+jdir) - fdst(i+1,j) +
     &             2.d0 * (fdst(i,j+jdir) - fdst(i,j))) +
     &           sigmaf(i-1,js-jdir) *
     &            (fdst(i-1,j-jdir) - fdst(i-1,j) +
     &             2.d0 * (fdst(i,j-jdir) - fdst(i,j)))))
      if (irz .eq. 1) then
         r = (i + 0.5d0) * hx
         rfac0 = hx * hym2 / (2.d0 * (r - hx))
         rfac1 = hx * hym2 / (2.d0 * r)
         res(i,j) = res(i,j) - fac0 *
     &     (rfac0 * sigmaf(i-1,js-jdir) * (fdst(i,j-jdir) - fdst(i,j)) -
     &      rfac1 * sigmaf(i,js)        * (fdst(i,j+jdir) - fdst(i,j)))
      end if
      fac0 = fac0 / ir
      do 10 m = 1, ir-1
         fac1 = (ir-m) * fac0
         tmp = hxm2 *
     &     (sigmaf(i+m-1,js) *
     &       (fdst(i+m-1,j+jdir) - fdst(i+m,j+jdir) +
     &        2.d0 * (fdst(i+m-1,j) - fdst(i+m,j))) +
     &      sigmaf(i+m,js) *
     &       (fdst(i+m+1,j+jdir) - fdst(i+m,j+jdir) +
     &        2.d0 * (fdst(i+m+1,j) - fdst(i+m,j))) +
     &      sigmaf(i-m-1,js-jdir) *
     &       (fdst(i-m-1,j-jdir) - fdst(i-m,j-jdir) +
     &        2.d0 * (fdst(i-m-1,j) - fdst(i-m,j))) +
     &      sigmaf(i-m,js-jdir) *
     &       (fdst(i-m+1,j-jdir) - fdst(i-m,j-jdir) +
     &        2.d0 * (fdst(i-m+1,j) - fdst(i-m,j))))
         res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     &     (sigmaf(i+m-1,js) *
     &       (fdst(i+m-1,j+jdir) - fdst(i+m-1,j) +
     &        2.d0 * (fdst(i+m,j+jdir) - fdst(i+m,j))) +
     &      sigmaf(i+m,js) *
     &       (fdst(i+m+1,j+jdir) - fdst(i+m+1,j) +
     &        2.d0 * (fdst(i+m,j+jdir) - fdst(i+m,j))) +
     &      sigmaf(i-m-1,js-jdir) *
     &       (fdst(i-m-1,j-jdir) - fdst(i-m-1,j) +
     &        2.d0 * (fdst(i-m,j-jdir) - fdst(i-m,j))) +
     &      sigmaf(i-m,js-jdir) *
     &       (fdst(i-m+1,j-jdir) - fdst(i-m+1,j) +
     &        2.d0 * (fdst(i-m,j-jdir) - fdst(i-m,j)))))
         if (irz .eq. 1) then
            r = (i - m + 0.5d0) * hx
            rfac0m = hx * hym2 / (2.d0 * (r - hx))
            rfac1m = hx * hym2 / (2.d0 * r)
            r = (i + m + 0.5d0) * hx
            rfac0p = hx * hym2 / (2.d0 * (r - hx))
            rfac1p = hx * hym2 / (2.d0 * r)
            res(i,j) = res(i,j) - fac1 *
     &       ((rfac0m * sigmaf(i-m-1,js-jdir) 
     &			- rfac1m * sigmaf(i-m,js-jdir)) *
     &         (fdst(i-m,j-jdir) - fdst(i-m,j)) +
     &        (rfac0p * sigmaf(i+m-1,js) - rfac1p * sigmaf(i+m,js)) *
     &          (fdst(i+m,j+jdir) - fdst(i+m,j)))
         end if
 10      continue
      fac0 = ir * fac0 / jr
      do 20 n = jdir, jdir*(jr-1), jdir
         fac1 = (jr-abs(n)) * fac0
         tmp = hxm2 *
     &     (sigmaf(i,j+n-1) *
     &       (fdst(i+1,j+n-1) - fdst(i,j+n-1) +
     &        2.d0 * (fdst(i+1,j+n) - fdst(i,j+n))) +
     &      sigmaf(i,j+n) *
     &       (fdst(i+1,j+n+1) - fdst(i,j+n+1) +
     &        2.d0 * (fdst(i+1,j+n) - fdst(i,j+n))) +
     &      sigmaf(i-1,j-n-1) *
     &       (fdst(i-1,j-n-1) - fdst(i,j-n-1) +
     &        2.d0 * (fdst(i-1,j-n) - fdst(i,j-n))) +
     &      sigmaf(i-1,j-n) *
     &       (fdst(i-1,j-n+1) - fdst(i,j-n+1) +
     &        2.d0 * (fdst(i-1,j-n) - fdst(i,j-n))))
         res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     &     (sigmaf(i,j+n-1) *
     &       (fdst(i+1,j+n-1) - fdst(i+1,j+n) +
     &        2.d0 * (fdst(i,j+n-1) - fdst(i,j+n))) +
     &      sigmaf(i,j+n) *
     &       (fdst(i+1,j+n+1) - fdst(i+1,j+n) +
     &        2.d0 * (fdst(i,j+n+1) - fdst(i,j+n))) +
     &      sigmaf(i-1,j-n-1) *
     &       (fdst(i-1,j-n-1) - fdst(i-1,j-n) +
     &        2.d0 * (fdst(i,j-n-1) - fdst(i,j-n))) +
     &      sigmaf(i-1,j-n) *
     &       (fdst(i-1,j-n+1) - fdst(i-1,j-n) +
     &        2.d0 * (fdst(i,j-n+1) - fdst(i,j-n)))))
         if (irz .eq. 1) then
            r = (i + 0.5d0) * hx
            rfac0 = hx * hym2 / (2.d0 * (r - hx))
            rfac1 = hx * hym2 / (2.d0 * r)
            res(i,j) = res(i,j) - fac1 *
     &       (rfac0 * (sigmaf(i-1,j-n-1) *
     &				(fdst(i,j-n-1) - fdst(i,j-n)) +
     &                 sigmaf(i-1,j-n)   * 
     &				(fdst(i,j-n+1) - fdst(i,j-n))) -
     &        rfac1 * (sigmaf(i,j+n-1) * (fdst(i,j+n-1) - fdst(i,j+n)) +
     &                 sigmaf(i,j+n)   * (fdst(i,j+n+1) - fdst(i,j+n))))
         end if
 20      continue
      end

      
      subroutine hgfres_terrain(res, resl0, resh0, resl1, resh1,
     @ src,  srcl0, srch0, srcl1, srch1,
     @ fdst, fdstl0, fdsth0, fdstl1, fdsth1,
     @ cdst, cdstl0, cdsth0, cdstl1, cdsth1,
     @ sigmaf, sfl0, sfh0, sfl1, sfh1,
     @ sigmac, scl0, sch0, scl1, sch1,
     @ regl0, regh0, regl1, regh1,
     @ ir, jr, idim, idir)
      integer resl0, resh0, resl1, resh1
      integer srcl0, srch0, srcl1, srch1
      integer fdstl0, fdsth0, fdstl1, fdsth1
      integer cdstl0, cdsth0, cdstl1, cdsth1
      integer sfl0, sfh0, sfl1, sfh1
      integer scl0, sch0, scl1, sch1
      integer regl0, regh0, regl1, regh1
      double precision res(resl0:resh0,resl1:resh1)
      double precision src(srcl0:srch0,srcl1:srch1)
      double precision fdst(fdstl0:fdsth0,fdstl1:fdsth1)
      double precision cdst(cdstl0:cdsth0,cdstl1:cdsth1)
      double precision sigmaf(sfl0:sfh0,sfl1:sfh1, 3)
      double precision sigmac(scl0:sch0,scl1:sch1, 3)
      integer ir, jr, idim, idir
      double precision fac0, fac1, tmp
      integer i, j, is, js, m, n
      if (idim .eq. 0) then
         i = regl0
         if (idir .eq. 1) then
            is = i - 1
         else
            is = i
         end if
         fac0 = 1.d0 / 6.d0
         do 10 j = regl1, regh1
 10         res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     @          (sigmac(is,j-1,1) *
     @            (cdst(i-idir,j-1) - cdst(i,j-1) +
     @             2.d0 * (cdst(i-idir,j) - cdst(i,j))) +
     @           sigmac(is,j,1) *
     @            (cdst(i-idir,j+1) - cdst(i,j+1) +
     @             2.d0 * (cdst(i-idir,j) - cdst(i,j))) +
     @           sigmac(is,j-1,2) *
     @            (cdst(i-idir,j-1) - cdst(i-idir,j) +
     @             2.d0 * (cdst(i,j-1) - cdst(i,j))) +
     @           sigmac(is,j,2) *
     @            (cdst(i-idir,j+1) - cdst(i-idir,j) +
     @             2.d0 * (cdst(i,j+1) - cdst(i,j))) -
     @           3.d0 * (sigmac(is,j-(i-is),3) *
     @                     (cdst(i-idir,j-idir) - cdst(i,j)) -
     @                   sigmac(is,j-(1-(i-is)),3) *
     @                     (cdst(i-idir,j+idir) - cdst(i,j))))
         fac0 = fac0 / jr
         i = i * ir
         if (idir .eq. 1) then
            is = i
         else
            is = i - 1
         end if
         do 50 n = 0, jr-1
            fac1 = (jr-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do 50 j = jr*regl1, jr*regh1, jr
               tmp =
     @            sigmaf(is,j-n-1,1) *
     @             (fdst(i+idir,j-n-1) - fdst(i,j-n-1) +
     @              2.d0 * (fdst(i+idir,j-n) - fdst(i,j-n))) +
     @		  sigmaf(is,j-n,1) *
     @             (fdst(i+idir,j-n+1) - fdst(i,j-n+1) +
     @              2.d0 * (fdst(i+idir,j-n) - fdst(i,j-n))) +
     @		  sigmaf(is,j+n-1,1) *
     @             (fdst(i+idir,j+n-1) - fdst(i,j+n-1) +
     @              2.d0 * (fdst(i+idir,j+n) - fdst(i,j+n))) +
     @		  sigmaf(is,j+n,1) *
     @             (fdst(i+idir,j+n+1) - fdst(i,j+n+1) +
     @              2.d0 * (fdst(i+idir,j+n) - fdst(i,j+n)))
               res(i,j) = res(i,j) - fac1 * (tmp +
     @		  sigmaf(is,j-n-1,2) *
     @             (fdst(i+idir,j-n-1) - fdst(i+idir,j-n) +
     @              2.d0 * (fdst(i,j-n-1) - fdst(i,j-n))) +
     @		  sigmaf(is,j-n,2) *
     @             (fdst(i+idir,j-n+1) - fdst(i+idir,j-n) +
     @              2.d0 * (fdst(i,j-n+1) - fdst(i,j-n))) +
     @		  sigmaf(is,j+n-1,2) *
     @             (fdst(i+idir,j+n-1) - fdst(i+idir,j+n) +
     @              2.d0 * (fdst(i,j+n-1) - fdst(i,j+n))) +
     @		  sigmaf(is,j+n,2) *
     @             (fdst(i+idir,j+n+1) - fdst(i+idir,j+n) +
     @              2.d0 * (fdst(i,j+n+1) - fdst(i,j+n))) -
     @		  3.d0 * (sigmaf(is,j-n-(i-is),3) *
     @             (fdst(i+idir,j-n+idir) - fdst(i+idir,j-n)) -
     @		          sigmaf(is,j-n-(1-(i-is)),3) *
     @             (fdst(i+idir,j-n-idir) - fdst(i+idir,j-n)) +
     @		          sigmaf(is,j+n-(i-is),3) *
     @             (fdst(i+idir,j+n+idir) - fdst(i+idir,j+n)) -
     @		          sigmaf(is,j+n-(1-(i-is)),3) *
     @             (fdst(i+idir,j+n-idir) - fdst(i+idir,j+n))))
 50         continue
      else
         j = regl1
         if (idir .eq. 1) then
            js = j - 1
         else
            js = j
         end if
         fac0 = 1.d0 / 6.d0
         do 60 i = regl0, regh0
 60         res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     @          (sigmac(i-1,js,1) *
     @            (cdst(i-1,j-idir) - cdst(i,j-idir) +
     @             2.d0 * (cdst(i-1,j) - cdst(i,j))) +
     @           sigmac(i,js,1) *
     @            (cdst(i+1,j-idir) - cdst(i,j-idir) +
     @             2.d0 * (cdst(i+1,j) - cdst(i,j))) +
     @           sigmac(i-1,js,2) *
     @            (cdst(i-1,j-idir) - cdst(i-1,j) +
     @             2.d0 * (cdst(i,j-idir) - cdst(i,j))) +
     @           sigmac(i,js,2) *
     @            (cdst(i+1,j-idir) - cdst(i+1,j) +
     @             2.d0 * (cdst(i,j-idir) - cdst(i,j))) -
     @           3.d0 * (sigmac(i-(j-js),js,3) *
     @                     (cdst(i-idir,j-idir) - cdst(i,j)) -
     @                   sigmac(i-(1-(j-js)),js,3) *
     @                     (cdst(i+idir,j-idir) - cdst(i,j))))
         fac0 = fac0 / ir
         j = j * jr
         if (idir .eq. 1) then
            js = j
         else
            js = j - 1
         end if
         do 100 m = 0, ir-1
            fac1 = (ir-m) * fac0
            if (m .eq. 0) fac1 = 0.5d0 * fac1
            do 100 i = ir*regl0, ir*regh0, ir
               tmp =
     @		  sigmaf(i-m-1,js,1) *
     @             (fdst(i-m-1,j+idir) - fdst(i-m,j+idir) +
     @              2.d0 * (fdst(i-m-1,j) - fdst(i-m,j))) +
     @		  sigmaf(i-m,js,1) *
     @             (fdst(i-m+1,j+idir) - fdst(i-m,j+idir) +
     @              2.d0 * (fdst(i-m+1,j) - fdst(i-m,j))) +
     @		  sigmaf(i+m-1,js,1) *
     @             (fdst(i+m-1,j+idir) - fdst(i+m,j+idir) +
     @              2.d0 * (fdst(i+m-1,j) - fdst(i+m,j))) +
     @		  sigmaf(i+m,js,1) *
     @             (fdst(i+m+1,j+idir) - fdst(i+m,j+idir) +
     @              2.d0 * (fdst(i+m+1,j) - fdst(i+m,j)))
               res(i,j) = res(i,j) - fac1 * (tmp +
     @            sigmaf(i-m-1,js,2) *
     @             (fdst(i-m-1,j+idir) - fdst(i-m-1,j) +
     @              2.d0 * (fdst(i-m,j+idir) - fdst(i-m,j))) +
     @		  sigmaf(i-m,js,2) *
     @             (fdst(i-m+1,j+idir) - fdst(i-m+1,j) +
     @              2.d0 * (fdst(i-m,j+idir) - fdst(i-m,j))) +
     @		  sigmaf(i+m-1,js,2) *
     @             (fdst(i+m-1,j+idir) - fdst(i+m-1,j) +
     @              2.d0 * (fdst(i+m,j+idir) - fdst(i+m,j))) +
     @		  sigmaf(i+m,js,2) *
     @             (fdst(i+m+1,j+idir) - fdst(i+m+1,j) +
     @              2.d0 * (fdst(i+m,j+idir) - fdst(i+m,j))) -
     @		  3.d0 * (sigmaf(i-m-(j-js),js,3) *
     @             (fdst(i-m+idir,j+idir) - fdst(i-m,j+idir)) -
     @		          sigmaf(i-m-(1-(j-js)),js,3) *
     @             (fdst(i-m-idir,j+idir) - fdst(i-m,j+idir)) +
     @		          sigmaf(i+m-(j-js),js,3) *
     @             (fdst(i+m+idir,j+idir) - fdst(i+m,j+idir)) -
     @		          sigmaf(i+m-(1-(j-js)),js,3) *
     @             (fdst(i+m-idir,j+idir) - fdst(i+m,j+idir))))
 100     continue
      end if
      end
      
      subroutine hgcres_terrain(res, resl0, resh0, resl1, resh1,
     @ src,  srcl0, srch0, srcl1, srch1,
     @ fdst, fdstl0, fdsth0, fdstl1, fdsth1,
     @ cdst, cdstl0, cdsth0, cdstl1, cdsth1,
     @ sigmaf, sfl0, sfh0, sfl1, sfh1,
     @ sigmac, scl0, sch0, scl1, sch1,
     @ regl0, regh0, regl1, regh1,
     @ ir, jr, ga)
      integer resl0, resh0, resl1, resh1
      integer srcl0, srch0, srcl1, srch1
      integer fdstl0, fdsth0, fdstl1, fdsth1
      integer cdstl0, cdsth0, cdstl1, cdsth1
      integer sfl0, sfh0, sfl1, sfh1
      integer scl0, sch0, scl1, sch1
      integer regl0, regh0, regl1, regh1
      double precision res(resl0:resh0,resl1:resh1)
      double precision src(srcl0:srch0,srcl1:srch1)
      double precision fdst(fdstl0:fdsth0,fdstl1:fdsth1)
      double precision cdst(cdstl0:cdsth0,cdstl1:cdsth1)
      double precision sigmaf(sfl0:sfh0,sfl1:sfh1, 3)
      double precision sigmac(scl0:sch0,scl1:sch1, 3)
      integer ir, jr, ga(0:1,0:1)
      double precision sum, fac, fac1
      integer ic, jc, if, jf, ii, ji, idir, jdir, m, n
      ic = regl0
      jc = regl1
      if = ic * ir
      jf = jc * jr
      sum = 0.d0
c quadrants
      do 10 ji = 0, 1
         jdir = 2 * ji - 1
         do 10 ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,ji) .eq. 1) then
               sum = sum + sigmaf(if+ii-1,jf+ji-1,1) *
     @            (fdst(if+idir,jf+jdir) - fdst(if,jf+jdir) +
     @             2.d0 * (fdst(if+idir,jf) - fdst(if,jf))) +
     @                     sigmaf(if+ii-1,jf+ji-1,2) *
     @            (fdst(if+idir,jf+jdir) - fdst(if+idir,jf) +
     @             2.d0 * (fdst(if,jf+jdir) - fdst(if,jf))) -
     @                     sigmaf(if+ii-1,jf+ji-1,3) * idir * jdir *
     @            3.d0 * (fdst(if+idir,jf+jdir) - fdst(if,jf))
            else
               sum = sum + sigmac(ic+ii-1,jc+ji-1,1) *
     @            (cdst(ic+idir,jc+jdir) - cdst(ic,jc+jdir) +
     @             2.d0 * (cdst(ic+idir,jc) - cdst(ic,jc))) +
     @                     sigmac(ic+ii-1,jc+ji-1,2) *
     @            (cdst(ic+idir,jc+jdir) - cdst(ic+idir,jc) +
     @             2.d0 * (cdst(ic,jc+jdir) - cdst(ic,jc))) -
     @                     sigmac(ic+ii-1,jc+ji-1,3) * idir * jdir *
     @            3.d0 * (cdst(ic+idir,jc+jdir) - cdst(ic,jc))
            endif
 10         continue
c edges
      do 40 ji = 0, 1
         jdir = 2 * ji - 1
         do 40 ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,ji) - ga(ii,1-ji) .eq. 1) then
               fac1 = 1.d0 / ir
               do 20 m = idir, idir*(ir-1), idir
                  fac = (ir-abs(m)) * fac1
 20               sum = sum + fac * (sigmaf(if+m-1,jf+ji-1,1) *
     @               (fdst(if+m-1,jf+jdir) - fdst(if+m,jf+jdir) +
     @                2.d0 * (fdst(if+m-1,jf) - fdst(if+m,jf))) +
     @                               sigmaf(if+m-1,jf+ji-1,2) *
     @               (fdst(if+m-1,jf+jdir) - fdst(if+m-1,jf) +
     @                2.d0 * (fdst(if+m,jf+jdir) - fdst(if+m,jf))) +
     @                               sigmaf(if+m-1,jf+ji-1,3) * jdir *
     @               3.d0 * (fdst(if+m-1,jf+jdir) - fdst(if+m,jf)) +
     @                               sigmaf(if+m,jf+ji-1,1) *
     @               (fdst(if+m+1,jf+jdir) - fdst(if+m,jf+jdir) +
     @                2.d0 * (fdst(if+m+1,jf) - fdst(if+m,jf))) +
     @                               sigmaf(if+m,jf+ji-1,2) *
     @               (fdst(if+m+1,jf+jdir) - fdst(if+m+1,jf) +
     @                2.d0 * (fdst(if+m,jf+jdir) - fdst(if+m,jf))) -
     @                               sigmaf(if+m,jf+ji-1,3) * jdir *
     @               3.d0 * (fdst(if+m+1,jf+jdir) - fdst(if+m,jf)))
            endif
            if (ga(ii,ji) - ga(1-ii,ji) .eq. 1) then
               fac1 = 1.d0 / jr
               do 30 n = jdir, jdir*(jr-1), jdir
                  fac = (jr-abs(n)) * fac1
 30               sum = sum + fac * (sigmaf(if+ii-1,jf+n-1,1) *
     @               (fdst(if+idir,jf+n-1) - fdst(if,jf+n-1) +
     @                2.d0 * (fdst(if+idir,jf+n) - fdst(if,jf+n))) +
     @                               sigmaf(if+ii-1,jf+n-1,2) *
     @               (fdst(if+idir,jf+n-1) - fdst(if+idir,jf+n) +
     @                2.d0 * (fdst(if,jf+n-1) - fdst(if,jf+n))) +
     @                               sigmaf(if+ii-1,jf+n-1,3) * idir *
     @               3.d0 * (fdst(if+idir,jf+n-1) - fdst(if,jf+n)) +
     @                               sigmaf(if+ii-1,jf+n,1) *
     @               (fdst(if+idir,jf+n+1) - fdst(if,jf+n+1) +
     @                2.d0 * (fdst(if+idir,jf+n) - fdst(if,jf+n))) +
     @                               sigmaf(if+ii-1,jf+n,2) *
     @               (fdst(if+idir,jf+n+1) - fdst(if+idir,jf+n) +
     @                2.d0 * (fdst(if,jf+n+1) - fdst(if,jf+n))) -
     @                               sigmaf(if+ii-1,jf+n,3) * idir *
     @               3.d0 * (fdst(if+idir,jf+n+1) - fdst(if,jf+n)))
            endif
 40         continue
c weighting
      res(if,jf) = src(if,jf) - sum / 6.d0
      end

      
c-----------------------------------------------------------------------
c NODE-based data, factor of 2 only.
      subroutine hgints_no_sigma_node(
     & dest, destl0, desth0, destl1, desth1,
     &       regl0, regh0, regl1, regh1,
     & sigx, sigy,
     &       sbl0, sbh0, sbl1, sbh1,
     & src,  srcl0, srch0, srcl1, srch1,
     &        bbl0, bbh0, bbl1, bbh1,
     & ir, jr)
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer sbl0, sbh0, sbl1, sbh1
      integer srcl0, srch0, srcl1, srch1
      integer bbl0, bbh0, bbl1, bbh1
      integer ir, jr
      double precision dest(destl0:desth0,destl1:desth1)
      double precision sigx(sbl0:sbh0,sbl1:sbh1)
      double precision sigy(sbl0:sbh0,sbl1:sbh1)
      double precision src(srcl0:srch0,srcl1:srch1)
      integer i, j, ic, jc
      do 10 jc = bbl1, bbh1
         do 10 ic = bbl0, bbh0
 10         dest(ir*ic,jr*jc) = src(ic,jc)
      if (ir .eq. 2) then
         do 20 jc = bbl1, bbh1
            do 20 ic = bbl0, bbh0-1
               i = ir * ic
               j = jr * jc
 20            dest(i+1,j) =
     &           ((sigx(i,j-1)+sigx(i,j)) * src(ic,jc) +
     &            (sigx(i+1,j-1)+sigx(i+1,j)) * src(ic+1,jc)) /
     &           (sigx(i,j-1)+sigx(i,j)+
     &            sigx(i+1,j-1)+sigx(i+1,j))
      end if
      if (jr .eq. 2) then
         do 30 jc = bbl1, bbh1-1
            do 30 ic = bbl0, bbh0
               i = ir * ic
               j = jr * jc
 30            dest(i,j+1) =
     &           ((sigy(i-1,j)+sigy(i,j)) * src(ic,jc) +
     &            (sigy(i-1,j+1)+sigy(i,j+1)) * src(ic,jc+1)) /
     &           (sigy(i-1,j)+sigy(i,j)+
     &            sigy(i-1,j+1)+sigy(i,j+1))
      end if
      if (ir .eq. 2 .and. jr .eq. 2) then
         do 40 jc = bbl1, bbh1-1
cdir$ ivdep
            do 40 ic = bbl0, bbh0-1
               i = ir * ic
               j = jr * jc
 40            dest(i+1,j+1) = ((sigx(i,j) + sigx(i,j+1)) *
     &                            dest(i,j+1) +
     &                          (sigx(i+1,j) + sigx(i+1,j+1)) *
     &                            dest(i+2,j+1) +
     &                          (sigy(i,j) + sigy(i+1,j)) *
     &                            dest(i+1,j) +
     &                          (sigy(i,j+1) + sigy(i+1,j+1)) *
     &                            dest(i+1,j+2)) /
     &                         (sigx(i,j) + sigx(i,j+1) +
     &                          sigx(i+1,j) + sigx(i+1,j+1) +
     &                          sigy(i,j) + sigy(i+1,j) +
     &                          sigy(i,j+1) + sigy(i+1,j+1))
      end if
      end
c-----------------------------------------------------------------------
c NODE-based data, factor of 2 only.
      subroutine hgints(
     & dest, destl0, desth0, destl1, desth1,
     &       regl0, regh0, regl1, regh1,
     & signd, sbl0, sbh0, sbl1, sbh1,
     & src,   srcl0, srch0, srcl1, srch1,
     &        bbl0, bbh0, bbl1, bbh1,
     & ir, jr)
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer sbl0, sbh0, sbl1, sbh1
      integer srcl0, srch0, srcl1, srch1
      integer bbl0, bbh0, bbl1, bbh1
      integer ir, jr
      double precision dest(destl0:desth0,destl1:desth1)
      double precision signd(sbl0:sbh0,sbl1:sbh1, 2)
      double precision src(srcl0:srch0,srcl1:srch1)
      stop 'no code for hgints (sigma node) in 2D'
      end
      
c-----------------------------------------------------------------------
c CELL-based data only.
      subroutine hgsrst(destx, desty,
     & destl0, desth0, destl1, desth1,
     & regl0, regh0, regl1, regh1,
     & srcx, srcy, srcl0, srch0, srcl1, srch1,
     & ir, jr)
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer srcl0, srch0, srcl1, srch1
      integer ir, jr
      double precision destx(destl0:desth0,destl1:desth1)
      double precision desty(destl0:desth0,destl1:desth1)
      double precision srcx(srcl0:srch0,srcl1:srch1)
      double precision srcy(srcl0:srch0,srcl1:srch1)
      integer i, j, i2, j2
      if (ir .eq. 2 .and. jr .eq. 2) then
         do 10 j = regl1, regh1
            do 10 i = regl0, regh0
               i2 = 2 * i
               j2 = 2 * j
               destx(i,j) = 1.d0 /
     &                      (1.d0 / (srcx(i2,j2)   + srcx(i2,j2+1)) +
     &                       1.d0 / (srcx(i2+1,j2) + srcx(i2+1,j2+1)))
               desty(i,j) = 1.d0 /
     &                      (1.d0 / (srcy(i2,j2)   + srcy(i2+1,j2)) +
     &                       1.d0 / (srcy(i2,j2+1) + srcy(i2+1,j2+1)))
 10            continue
      else if (ir .eq. 2) then
         do 20 j = regl1, regh1
            do 20 i = regl0, regh0
               i2 = 2 * i
               destx(i,j) = 2.d0 /
     &                      (1.d0 / srcx(i2,j) + 1.d0 / srcx(i2+1,j))
               desty(i,j) = 0.5d0 * (srcy(i2,j) + srcy(i2+1,j))
 20            continue
      else
         do 30 j = regl1, regh1
            do 30 i = regl0, regh0
               j2 = 2 * j
               destx(i,j) = 0.5d0 * (srcx(i,j2) + srcx(i,j2+1))
               desty(i,j) = 2.d0 /
     &                      (1.d0 / srcy(i,j2) + 1.d0 / srcy(i,j2+1))
 30            continue
      end if
      end
      
c five-point variable stencils
c-----------------------------------------------------------------------
      subroutine hgrlx(cor, corl0, corh0, corl1, corh1,
     & res, resl0, resh0, resl1, resh1,
     & sigx, sigy, sfl0, sfh0, sfl1, sfh1,
     & cen, cenl0, cenh0, cenl1, cenh1,
     & regl0, regh0, regl1, regh1,
     & hx, hy, irz, imax)
      integer corl0, corh0, corl1, corh1
      integer resl0, resh0, resl1, resh1
      integer sfl0, sfh0, sfl1, sfh1
      integer cenl0, cenh0, cenl1, cenh1
      integer regl0, regh0, regl1, regh1
      double precision cor(corl0:corh0,corl1:corh1)
      double precision res(resl0:resh0,resl1:resh1)
      double precision sigx(sfl0:sfh0,sfl1:sfh1)
      double precision sigy(sfl0:sfh0,sfl1:sfh1)
      double precision cen(cenl0:cenh0,cenl1:cenh1)
      double precision hx, hy
      integer irz, imax
      double precision hxm2, hym2
      integer i, j, ipass, ipar
      double precision AVG
      AVG() = 0.5d0 * (hxm2 *
     &          ((sigx(i-1,j-1) + sigx(i-1,j)) * cor(i-1,j) +
     &           (sigx(i,j-1) + sigx(i,j)) * cor(i+1,j)) +
     &                 hym2 *
     &          ((sigy(i-1,j-1) + sigy(i,j-1)) * cor(i,j-1) +
     &           (sigy(i-1,j) + sigy(i,j)) * cor(i,j+1)))
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      if (regh1 - regl1 .gt. regh0 - regl0) then
         do 10 ipass = 1, 0, -1
            ipar = ipass
            do 10 i = regl0, regh0
               ipar = 1 - ipar
cdir$ ivdep
               do 10 j = regl1 + ipar, regh1, 2
                  cor(i,j) = (AVG() - res(i,j)) * cen(i,j)
 10               continue
      else
         do 20 ipass = 1, 0, -1
            ipar = ipass
            do 20 j = regl1, regh1
               ipar = 1 - ipar
cdir$ ivdep
               do 20 i = regl0 + ipar, regh0, 2
                  cor(i,j) = (AVG() - res(i,j)) * cen(i,j)
 20               continue
      end if
      end
c-----------------------------------------------------------------------
      subroutine hgrlxl(cor, corl0, corh0, corl1, corh1,
     & res, resl0, resh0, resl1, resh1,
     & sigx, sigy, sfl0, sfh0, sfl1, sfh1,
     & cen, cenl0, cenh0, cenl1, cenh1,
     & regl0, regh0, regl1, regh1,
     & doml0, domh0, doml1, domh1,
     & hx, hy, irz, imax, lsd)
      integer corl0, corh0, corl1, corh1
      integer resl0, resh0, resl1, resh1
      integer sfl0, sfh0, sfl1, sfh1
      integer cenl0, cenh0, cenl1, cenh1
      integer regl0, regh0, regl1, regh1
      integer doml0, domh0, doml1, domh1
      double precision cor(corl0:corh0,corl1:corh1)
      double precision res(resl0:resh0,resl1:resh1)
      double precision sigx(sfl0:sfh0,sfl1:sfh1)
      double precision sigy(sfl0:sfh0,sfl1:sfh1)
      double precision cen(cenl0:cenh0,cenl1:cenh1)
      double precision hx, hy
      double precision wrk(256)
      integer irz, imax, lsd
      double precision hxm2, hym2, betm, aj
      double precision RHSL0, RHSL1
      integer i, j, jw, ipass
      RHSL0()=(res(i,j) - 0.5d0 * hym2 * 
     &          ((sigy(i-1,j-1) + sigy(i,j-1)) * cor(i,j-1) + 
     &           (sigy(i-1,j) + sigy(i,j)) * cor(i,j+1)))

      RHSL1()=(res(i,j) - 0.5d0 * hxm2 * 
     &          ((sigx(i-1,j-1) + sigx(i-1,j)) * cor(i-1,j) + 
     &           (sigx(i,j-1) + sigx(i,j)) * cor(i+1,j)))

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
               cor(i,j) = RHSL1()* betm
               aj = 0.5d0 * hym2 * (sigy(i-1,j) + sigy(i,j))
               wrk(1) = 2.d0 * aj * betm
            else
c interface to grid at same level:
               aj = 0.5d0 * hym2 * (sigy(i-1,j-1) + sigy(i,j-1))
               cor(i,j) = (RHSL1() - aj * cor(i,j-1)) * betm
               aj = 0.5d0 * hym2 * (sigy(i-1,j) + sigy(i,j))
               wrk(1) = aj * betm
            end if
c forward solve loop:
            do 10 j = regl1 + 1, regh1 - 1
               jw = j - regl1
               if (cen(i,j) .eq. 0.d0) then
                  betm = 0.d0
               else
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(jw))
               end if
               cor(i,j) = (RHSL1() - aj * cor(i,j-1)) * betm
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
               cor(i,j) = (RHSL1() - aj * cor(i,j-1)) * betm
            else if (jw .gt. 0) then
c interface to grid at same level:
               betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(jw))
               cor(i,j) = RHSL1() - aj * cor(i,j-1)
               aj = 0.5d0 * hym2 * (sigy(i-1,j) + sigy(i,j))
               cor(i,j) = (cor(i,j) - aj * cor(i,j+1)) * betm
            end if
c back substitution loop:
            do 20 j = regh1 - 1, regl1, -1
               jw = j - regl1
               cor(i,j) = cor(i,j) - wrk(jw + 1) * cor(i,j+1)
 20            continue
      else
         print *, "Line solve not implemented in dimension", lsd
	    stop
      end if
      end
c-----------------------------------------------------------------------
      subroutine hgrlnf(cor, corl0, corh0, corl1, corh1,
     & res, resl0, resh0, resl1, resh1,
     & wrk, wrkl0, wrkh0, wrkl1, wrkh1,
     & sigx, sigy, sfl0, sfh0, sfl1, sfh1,
     & cen, cenl0, cenh0, cenl1, cenh1,
     & regl0, regh0, regl1, regh1,
     & doml0, domh0, doml1, domh1,
     & hx, hy, irz, imax, lsd, ipass)
      integer corl0, corh0, corl1, corh1
      integer resl0, resh0, resl1, resh1
      integer wrkl0, wrkh0, wrkl1, wrkh1
      integer sfl0, sfh0, sfl1, sfh1
      integer cenl0, cenh0, cenl1, cenh1
      integer regl0, regh0, regl1, regh1
      integer doml0, domh0, doml1, domh1
      double precision cor(corl0:corh0,corl1:corh1)
      double precision res(resl0:resh0,resl1:resh1)
      double precision wrk(wrkl0:wrkh0,wrkl1:wrkh1)
      double precision sigx(sfl0:sfh0,sfl1:sfh1)
      double precision sigy(sfl0:sfh0,sfl1:sfh1)
      double precision cen(cenl0:cenh0,cenl1:cenh1)
      double precision hx, hy
      integer irz, imax, lsd, ipass
      double precision hxm2, hym2, betm, aj
      double precision RHSL0, RHSL1
      integer i, j, ioff
      RHSL0()=(res(i,j) - 0.5d0 * hym2 * 
     &          ((sigy(i-1,j-1) + sigy(i,j-1)) * cor(i,j-1) + 
     &           (sigy(i-1,j) + sigy(i,j)) * cor(i,j+1)))

      RHSL1()=(res(i,j) - 0.5d0 * hxm2 * 
     &          ((sigx(i-1,j-1) + sigx(i-1,j)) * cor(i-1,j) + 
     &           (sigx(i,j-1) + sigx(i,j)) * cor(i+1,j)))
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      if (lsd .eq. 0) then
         if (mod(regl1,2) .eq. 0) then
            ioff = ipass
         else
            ioff = 1 - ipass
         end if
         i = regl0
         do 10 j = regl1 + ioff, regh1, 2
            aj = 0.5d0 * hxm2 * (sigx(i,j-1) + sigx(i,j))
            if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
               wrk(i,j) = 0.d0
            else if (regl0 .eq. doml0) then
c neumann bdy:
               betm = -cen(i,j)
               cor(i,j) = RHSL0() * betm
               wrk(i,j) = 2.d0 * aj * betm
            end if
            wrk(i+1,j) = aj
 10         continue
c forward solve loop:
         do 20 i = regl0 + 1, regh0 - 1
            do 20 j = regl1 + ioff, regh1, 2
               aj = wrk(i,j)
               if (cen(i,j) .eq. 0.d0) then
                  betm = 0.d0
               else
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i-1,j))
               end if
               cor(i,j) = (RHSL0() - aj * cor(i-1,j)) * betm
               aj = 0.5d0 * hxm2 * (sigx(i,j-1) + sigx(i,j))
               wrk(i+1,j) = aj
 20            wrk(i,j) = aj * betm
         i = regh0
         do 30 j = regl1 + ioff, regh1, 2
            if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
            else if (regh0 .eq. domh0) then
c neumann bdy:
               aj = 2.d0 * wrk(i,j)
               betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i-1,j))
               cor(i,j) = (RHSL0() - aj * cor(i-1,j)) * betm
            else if (i .gt. regl0) then
c interface to grid at same level:
               aj = wrk(i,j)
               betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i-1,j))
               cor(i,j) = (RHSL0() - aj * cor(i-1,j)) * betm
               aj = 0.5d0 * hxm2 * (sigx(i,j-1) + sigx(i,j))
               wrk(i,j) = aj * betm
            end if
 30         continue
      else
         if (mod(regl0,2) .eq. 0) then
            ioff = ipass
         else
            ioff = 1 - ipass
         end if
         j = regl1
         do 110 i = regl0 + ioff, regh0, 2
            aj = 0.5d0 * hym2 * (sigy(i-1,j) + sigy(i,j))
            if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
               wrk(i,j) = 0.d0
            else if (regl1 .eq. doml1) then
c neumann bdy:
               betm = -cen(i,j)
               cor(i,j) = RHSL1() * betm
               wrk(i,j) = 2.d0 * aj * betm
            end if
            wrk(i,j+1) = aj
 110        continue
c forward solve loop:
         do 120 j = regl1 + 1, regh1 - 1
            do 120 i = regl0 + ioff, regh0, 2
               aj = wrk(i,j)
               if (cen(i,j) .eq. 0.d0) then
                  betm = 0.d0
               else
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
               end if
               cor(i,j) = (RHSL1() - aj * cor(i,j-1)) * betm
               aj = 0.5d0 * hym2 * (sigy(i-1,j) + sigy(i,j))
               wrk(i,j+1) = aj
 120           wrk(i,j) = aj * betm
         j = regh1
         do 130 i = regl0 + ioff, regh0, 2
            if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
            else if (regh1 .eq. domh1) then
c neumann bdy:
               aj = 2.d0 * wrk(i,j)
               betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
               cor(i,j) = (RHSL1() - aj * cor(i,j-1)) * betm
            else if (j .gt. regl1) then
c interface to grid at same level:
               aj = wrk(i,j)
               betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
               cor(i,j) = (RHSL1() - aj * cor(i,j-1)) * betm
               aj = 0.5d0 * hym2 * (sigy(i-1,j) + sigy(i,j))
               wrk(i,j) = aj * betm
            end if
 130        continue
      end if
      end
c-----------------------------------------------------------------------
      subroutine hgres(res, resl0, resh0, resl1, resh1,
     & src,  srcl0, srch0, srcl1, srch1,
     & dest, destl0, desth0, destl1, desth1,
     & sigx, sigy, sfl0, sfh0, sfl1, sfh1,
     & regl0, regh0, regl1, regh1,
     & hx, hy)
      integer resl0, resh0, resl1, resh1
      integer srcl0, srch0, srcl1, srch1
      integer destl0, desth0, destl1, desth1
      integer sfl0, sfh0, sfl1, sfh1
      integer regl0, regh0, regl1, regh1
      double precision res(resl0:resh0,resl1:resh1)
      double precision src(srcl0:srch0,srcl1:srch1)
      double precision dest(destl0:desth0,destl1:desth1)
      double precision sigx(sfl0:sfh0,sfl1:sfh1)
      double precision sigy(sfl0:sfh0,sfl1:sfh1)
      double precision hx, hy
      double precision hxm2, hym2
      integer i, j
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      if (regh1 - regl1 .gt. regh0 - regl0) then
         do 10 i = regl0, regh0
            do 10 j = regl1, regh1
 10            res(i,j) = src(i,j) - 0.5d0 * (hxm2 *
     &           ((sigx(i-1,j-1) + sigx(i-1,j)) *
     &             (dest(i-1,j) - dest(i,j)) +
     &            (sigx(i,j-1) + sigx(i,j)) *
     &             (dest(i+1,j) - dest(i,j))) +
     &                                        hym2 *
     &           ((sigy(i-1,j-1) + sigy(i,j-1)) *
     &             (dest(i,j-1) - dest(i,j)) +
     &            (sigy(i-1,j) + sigy(i,j)) *
     &             (dest(i,j+1) - dest(i,j))))
      else
         do 20 j = regl1, regh1
            do 20 i = regl0, regh0
 20            res(i,j) = src(i,j) - 0.5d0 * (hxm2 *
     &           ((sigx(i-1,j-1) + sigx(i-1,j)) *
     &             (dest(i-1,j) - dest(i,j)) +
     &            (sigx(i,j-1) + sigx(i,j)) *
     &             (dest(i+1,j) - dest(i,j))) +
     &                                        hym2 *
     &           ((sigy(i-1,j-1) + sigy(i,j-1)) *
     &             (dest(i,j-1) - dest(i,j)) +
     &            (sigy(i-1,j) + sigy(i,j)) *
     &             (dest(i,j+1) - dest(i,j))))
      end if
      end
      
c five-point variable stencils
      subroutine hgcen_no_sigma_node(
     @ cen, cenl0,cenh0,cenl1,cenh1,
     @ sigx, sigy,
     @       sbl0,sbh0,sbl1,sbh1,
     @       regl0,regh0,regl1,regh1,
     @ hx, hy)
      integer cenl0,cenh0,cenl1,cenh1
      integer sbl0,sbh0,sbl1,sbh1
      integer regl0,regh0,regl1,regh1
      double precision cen(cenl0:cenh0,cenl1:cenh1)
      double precision sigx(sbl0:sbh0,sbl1:sbh1)
      double precision sigy(sbl0:sbh0,sbl1:sbh1)
      double precision hx, hy
      double precision hxm2, hym2
      integer i, j
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      do 10 j = regl1, regh1
         do 10 i = regl0, regh0
 10         cen(i,j) = 2.d0 / (hxm2 * (sigx(i-1,j-1) + sigx(i-1,j) +
     @                                 sigx(i,j-1) + sigx(i,j)) +
     @                         hym2 * (sigy(i-1,j-1) + sigy(i-1,j) +
     @                                 sigy(i,j-1) + sigy(i,j)))
      end
      
      subroutine hgcen(
     @ cen, cenl0,cenh0,cenl1,cenh1,
     @ signd,
     @       sbl0,sbh0,sbl1,sbh1,
     @       regl0,regh0,regl1,regh1,
     @ hx, hy)
      integer cenl0,cenh0,cenl1,cenh1
      integer sbl0,sbh0,sbl1,sbh1
      integer regl0,regh0,regl1,regh1
      double precision cen(cenl0:cenh0,cenl1:cenh1)
      double precision signd(sbl0:sbh0,sbl1:sbh1, 2)
      double precision hx, hy
      double precision hxm2, hym2
      integer i, j
      stop 'no code for hgcen (sigma_node) in 2D'
      end
      
c nine-point terrain stencils
c-----------------------------------------------------------------------
      subroutine hgcen_terrain(cen, cenl0, cenh0, cenl1, cenh1,
     & sig, sbl0, sbh0, sbl1, sbh1,
     & regl0, regh0, regl1, regh1)
      integer cenl0, cenh0, cenl1, cenh1
      integer sbl0, sbh0, sbl1, sbh1
      integer regl0, regh0, regl1, regh1
      double precision cen(cenl0:cenh0,cenl1:cenh1)
      double precision sig(sbl0:sbh0,sbl1:sbh1, 3)
      integer i, j
      do 10 j = regl1, regh1
         do 10 i = regl0, regh0
 10         cen(i,j) = 6.d0 / (2.d0 * (sig(i-1,j-1,1) + sig(i-1,j,1) +
     &                                 sig(i,j-1,1) + sig(i,j,1) +
     &                                 sig(i-1,j-1,2) + sig(i-1,j,2) +
     &                                 sig(i,j-1,2) + sig(i,j,2)) +
     &                         3.d0 * (sig(i-1,j,3) - sig(i-1,j-1,3) +
     &                                 sig(i,j-1,3) - sig(i,j,3)))
      end
c-----------------------------------------------------------------------
      subroutine hgrlx_terrain(cor, corl0, corh0, corl1, corh1,
     & res, resl0, resh0, resl1, resh1,
     & sig, sfl0, sfh0, sfl1, sfh1,
     & cen, cenl0, cenh0, cenl1, cenh1,
     & regl0, regh0, regl1, regh1)
      integer corl0, corh0, corl1, corh1
      integer resl0, resh0, resl1, resh1
      integer sfl0, sfh0, sfl1, sfh1
      integer cenl0, cenh0, cenl1, cenh1
      integer regl0, regh0, regl1, regh1
      double precision cor(corl0:corh0,corl1:corh1)
      double precision res(resl0:resh0,resl1:resh1)
      double precision sig(sfl0:sfh0,sfl1:sfh1, 3)
      double precision cen(cenl0:cenh0,cenl1:cenh1)
      double precision fac
      double precision AVG
      integer i, j
      AVG()= fac * (
     &           sig(i-1,j-1,1) *
     &           (cor(i-1,j-1) - cor(i,j-1) + 2.d0 * cor(i-1,j)) +
     &           sig(i-1,j,1) *
     &            (cor(i-1,j+1) - cor(i,j+1) + 2.d0 * cor(i-1,j)) +
     &           sig(i,j-1,1) *
     &            (cor(i+1,j-1) - cor(i,j-1) + 2.d0 * cor(i+1,j)) +
     &           sig(i,j,1) *
     &            (cor(i+1,j+1) - cor(i,j+1) + 2.d0 * cor(i+1,j)) +
     &           sig(i-1,j-1,2) *
     &            (cor(i-1,j-1) - cor(i-1,j) + 2.d0 * cor(i,j-1)) +
     &           sig(i-1,j,2) *
     &            (cor(i-1,j+1) - cor(i-1,j) + 2.d0 * cor(i,j+1)) +
     &           sig(i,j-1,2) *
     &            (cor(i+1,j-1) - cor(i+1,j) + 2.d0 * cor(i,j-1)) +
     &           sig(i,j,2) *
     &            (cor(i+1,j+1) - cor(i+1,j) + 2.d0 * cor(i,j+1)) -
     &           3.d0 * (sig(i-1,j-1,3) * cor(i-1,j-1) -
     &                   sig(i-1,j,3)   * cor(i-1,j+1) -
     &                   sig(i,j-1,3)   * cor(i+1,j-1) +
     &                   sig(i,j,3)     * cor(i+1,j+1)))

      fac = 1.d0 / 6.d0
      if (regh1 - regl1 .gt. regh0 - regl0) then
         do 10 i = regl0, regh0
cdir$ ivdep
            do 10 j = regl1, regh1
               cor(i,j) = (AVG() - res(i,j)) * cen(i,j)
 10            continue
      else
         do 20 j = regl1, regh1
cdir$ ivdep
            do 20 i = regl0, regh0
               cor(i,j) = (AVG() - res(i,j)) * cen(i,j)
 20            continue
      end if
      end
c-----------------------------------------------------------------------

      subroutine hgrlnf_terrain(cor, corl0, corh0, corl1, corh1,
     & res, resl0, resh0, resl1, resh1,
     & wrk, wrkl0, wrkh0, wrkl1, wrkh1,
     & sig, sfl0, sfh0, sfl1, sfh1,
     & cen, cenl0, cenh0, cenl1, cenh1,
     & regl0, regh0, regl1, regh1,
     & doml0, domh0, doml1, domh1,
     & lsd, ipass)
      integer corl0, corh0, corl1, corh1
      integer resl0, resh0, resl1, resh1
      integer wrkl0, wrkh0, wrkl1, wrkh1
      integer sfl0, sfh0, sfl1, sfh1
      integer cenl0, cenh0, cenl1, cenh1
      integer regl0, regh0, regl1, regh1
      integer doml0, domh0, doml1, domh1
      double precision cor(corl0:corh0,corl1:corh1)
      double precision res(resl0:resh0,resl1:resh1)
      double precision wrk(wrkl0:wrkh0,wrkl1:wrkh1)
      double precision sig(sfl0:sfh0,sfl1:sfh1, 3)
      double precision cen(cenl0:cenh0,cenl1:cenh1)
      integer lsd, ipass
      double precision fac, betm, aj
      double precision RHSL0,RHSL1
      integer i, j, ioff
      RHSL0() = (res(i,j) - fac * ( 
     &           sig(i-1,j-1,1) * 
     &             (cor(i-1,j-1) - cor(i,j-1)) + 
     &            sig(i-1,j,1) * 
     &             (cor(i-1,j+1) - cor(i,j+1)) + 
     &            sig(i,j-1,1) * 
     &             (cor(i+1,j-1) - cor(i,j-1)) + 
     &            sig(i,j,1) * 
     &             (cor(i+1,j+1) - cor(i,j+1)) + 
     &            sig(i-1,j-1,2) * 
     &             (cor(i-1,j-1) + 2.d0 * cor(i,j-1)) + 
     &            sig(i-1,j,2) * 
     &             (cor(i-1,j+1) + 2.d0 * cor(i,j+1)) + 
     &            sig(i,j-1,2) * 
     &             (cor(i+1,j-1) + 2.d0 * cor(i,j-1)) + 
     &            sig(i,j,2) * 
     &             (cor(i+1,j+1) + 2.d0 * cor(i,j+1)) - 
     &            3.d0 * (sig(i-1,j-1,3) * cor(i-1,j-1) - 
     &                    sig(i-1,j,3)   * cor(i-1,j+1) - 
     &                    sig(i,j-1,3)   * cor(i+1,j-1) + 
     &                    sig(i,j,3)     * cor(i+1,j+1))))

      RHSL1() = (res(i,j) - fac * ( 
     &           sig(i-1,j-1,1) * 
     &             (cor(i-1,j-1) + 2.d0 * cor(i-1,j)) + 
     &            sig(i-1,j,1) * 
     &             (cor(i-1,j+1) + 2.d0 * cor(i-1,j)) + 
     &            sig(i,j-1,1) * 
     &             (cor(i+1,j-1) + 2.d0 * cor(i+1,j)) + 
     &            sig(i,j,1) * 
     &             (cor(i+1,j+1) + 2.d0 * cor(i+1,j)) + 
     &            sig(i-1,j-1,2) * 
     &             (cor(i-1,j-1) - cor(i-1,j)) + 
     &            sig(i-1,j,2) * 
     &             (cor(i-1,j+1) - cor(i-1,j)) + 
     &            sig(i,j-1,2) * 
     &            (cor(i+1,j-1) - cor(i+1,j)) + 
     &            sig(i,j,2) * 
     &             (cor(i+1,j+1) - cor(i+1,j)) - 
     &            3.d0 * (sig(i-1,j-1,3) * cor(i-1,j-1) - 
     &                   sig(i-1,j,3)   * cor(i-1,j+1) - 
     &                    sig(i,j-1,3)   * cor(i+1,j-1) + 
     &                    sig(i,j,3)     * cor(i+1,j+1))))
      fac = 1.d0 / 6.d0
      if (lsd .eq. 0) then
         if (mod(regl1,2) .eq. 0) then
            ioff = ipass
         else
            ioff = 1 - ipass
         end if
         i = regl0
         do 10 j = regl1 + ioff, regh1, 2
            aj = fac *
     &        (2.d0 * (sig(i,j-1,1) + sig(i,j,1)) -
     &         (sig(i,j-1,2) + sig(i,j,2)))
            if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
               wrk(i,j) = 0.d0
            else if (regl0 .eq. doml0) then
c neumann bdy:
               betm = -cen(i,j)
               cor(i,j) = RHSL0() * betm
               wrk(i,j) = 2.d0 * aj * betm
            end if
            wrk(i+1,j) = aj
 10         continue
c forward solve loop:
         do 20 i = regl0 + 1, regh0 - 1
            do 20 j = regl1 + ioff, regh1, 2
               aj = wrk(i,j)
               if (cen(i,j) .eq. 0.d0) then
                  betm = 0.d0
               else
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i-1,j))
               end if
               cor(i,j) = (RHSL0() - aj * cor(i-1,j)) * betm
               aj = fac *
     &           (2.d0 * (sig(i,j-1,1) + sig(i,j,1)) -
     &            (sig(i,j-1,2) + sig(i,j,2)))
               wrk(i+1,j) = aj
 20            wrk(i,j) = aj * betm
         i = regh0
         do 30 j = regl1 + ioff, regh1, 2
            if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
            else if (regh0 .eq. domh0) then
c neumann bdy:
               aj = 2.d0 * wrk(i,j)
               betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i-1,j))
               cor(i,j) = (RHSL0() - aj * cor(i-1,j)) * betm
            else if (i .gt. regl0) then
c interface to grid at same level:
               aj = wrk(i,j)
               betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i-1,j))
               cor(i,j) = (RHSL0() - aj * cor(i-1,j)) * betm
               aj = fac *
     &           (2.d0 * (sig(i,j-1,1) + sig(i,j,1)) -
     &            (sig(i,j-1,2) + sig(i,j,2)))
               wrk(i,j) = aj * betm
            end if
 30         continue
      else
         if (mod(regl0,2) .eq. 0) then
            ioff = ipass
         else
            ioff = 1 - ipass
         end if
         j = regl1
         do 110 i = regl0 + ioff, regh0, 2
            aj = fac *
     &        (2.d0 * (sig(i-1,j,2) + sig(i,j,2)) -
     &         (sig(i-1,j,1) + sig(i,j,1)))
            if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
               wrk(i,j) = 0.d0
            else if (regl1 .eq. doml1) then
c neumann bdy:
               betm = -cen(i,j)
               cor(i,j) = RHSL1() * betm
               wrk(i,j) = 2.d0 * aj * betm
            end if
            wrk(i,j+1) = aj
 110        continue
c forward solve loop:
         do 120 j = regl1 + 1, regh1 - 1
            do 120 i = regl0 + ioff, regh0, 2
               aj = wrk(i,j)
               if (cen(i,j) .eq. 0.d0) then
                  betm = 0.d0
               else
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
               end if
               cor(i,j) = (RHSL1() - aj * cor(i,j-1)) * betm
               aj = fac *
     &           (2.d0 * (sig(i-1,j,2) + sig(i,j,2)) -
     &            (sig(i-1,j,1) + sig(i,j,1)))
               wrk(i,j+1) = aj
 120           wrk(i,j) = aj * betm
         j = regh1
         do 130 i = regl0 + ioff, regh0, 2
            if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
            else if (regh1 .eq. domh1) then
c neumann bdy:
               aj = 2.d0 * wrk(i,j)
               betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
               cor(i,j) = (RHSL1() - aj * cor(i,j-1)) * betm
            else if (j .gt. regl1) then
c interface to grid at same level:
               aj = wrk(i,j)
               betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
               cor(i,j) = (RHSL1() - aj * cor(i,j-1)) * betm
               aj = fac *
     &           (2.d0 * (sig(i-1,j,2) + sig(i,j,2)) -
     &            (sig(i-1,j,1) + sig(i,j,1)))
               wrk(i,j) = aj * betm
            end if
 130        continue
      end if
      end
c-----------------------------------------------------------------------
      subroutine hgres_terrain(res, resl0, resh0, resl1, resh1,
     & src,  srcl0, srch0, srcl1, srch1,
     & dest, destl0, desth0, destl1, desth1,
     & sig, sfl0, sfh0, sfl1, sfh1,
     & cen, cenl0, cenh0, cenl1, cenh1,
     & regl0, regh0, regl1, regh1)
      integer resl0, resh0, resl1, resh1
      integer srcl0, srch0, srcl1, srch1
      integer destl0, desth0, destl1, desth1
      integer sfl0, sfh0, sfl1, sfh1
      integer cenl0, cenh0, cenl1, cenh1
      integer regl0, regh0, regl1, regh1
      double precision res(resl0:resh0,resl1:resh1)
      double precision src(srcl0:srch0,srcl1:srch1)
      double precision dest(destl0:desth0,destl1:desth1)
      double precision sig(sfl0:sfh0,sfl1:sfh1, 3)
      double precision cen(cenl0:cenh0,cenl1:cenh1)
      double precision fac, tmp
      integer i, j
      fac = 1.d0 / 6.d0
      if (regh1 - regl1 .gt. regh0 - regl0) then
         do 10 i = regl0, regh0
            do 10 j = regl1, regh1
               tmp =
     &           (sig(i-1,j-1,1) *
     &             (dest(i-1,j-1) - dest(i,j-1) +
     &              2.d0 * (dest(i-1,j) - dest(i,j))) +
     &            sig(i-1,j,1) *
     &             (dest(i-1,j+1) - dest(i,j+1) +
     &              2.d0 * (dest(i-1,j) - dest(i,j))) +
     &            sig(i,j-1,1) *
     &             (dest(i+1,j-1) - dest(i,j-1) +
     &              2.d0 * (dest(i+1,j) - dest(i,j))) +
     &            sig(i,j,1) *
     &             (dest(i+1,j+1) - dest(i,j+1) +
     &              2.d0 * (dest(i+1,j) - dest(i,j))))
               tmp = src(i,j) - fac * (tmp +
     &           (sig(i-1,j-1,2) *
     &             (dest(i-1,j-1) - dest(i-1,j) +
     &              2.d0 * (dest(i,j-1) - dest(i,j))) +
     &            sig(i-1,j,2) *
     &             (dest(i-1,j+1) - dest(i-1,j) +
     &              2.d0 * (dest(i,j+1) - dest(i,j))) +
     &            sig(i,j-1,2) *
     &             (dest(i+1,j-1) - dest(i+1,j) +
     &              2.d0 * (dest(i,j-1) - dest(i,j))) +
     &            sig(i,j,2) *
     &             (dest(i+1,j+1) - dest(i+1,j) +
     &              2.d0 * (dest(i,j+1) - dest(i,j)))))
 10            res(i,j) = tmp + 0.5d0 *
     &           (sig(i-1,j-1,3) *
     &             (dest(i-1,j-1) - dest(i,j)) -
     &            sig(i-1,j,3) *
     &             (dest(i-1,j+1) - dest(i,j)) -
     &            sig(i,j-1,3) *
     &             (dest(i+1,j-1) - dest(i,j)) +
     &            sig(i,j,3) *
     &             (dest(i+1,j+1) - dest(i,j)))
      else
         do 20 j = regl1, regh1
            do 20 i = regl0, regh0
               tmp =
     &           (sig(i-1,j-1,1) *
     &             (dest(i-1,j-1) - dest(i,j-1) +
     &              2.d0 * (dest(i-1,j) - dest(i,j))) +
     &            sig(i-1,j,1) *
     &             (dest(i-1,j+1) - dest(i,j+1) +
     &              2.d0 * (dest(i-1,j) - dest(i,j))) +
     &            sig(i,j-1,1) *
     &             (dest(i+1,j-1) - dest(i,j-1) +
     &              2.d0 * (dest(i+1,j) - dest(i,j))) +
     &            sig(i,j,1) *
     &             (dest(i+1,j+1) - dest(i,j+1) +
     &              2.d0 * (dest(i+1,j) - dest(i,j))))
               tmp = src(i,j) - fac * (tmp +
     &           (sig(i-1,j-1,2) *
     &             (dest(i-1,j-1) - dest(i-1,j) +
     &              2.d0 * (dest(i,j-1) - dest(i,j))) +
     &            sig(i-1,j,2) *
     &             (dest(i-1,j+1) - dest(i-1,j) +
     &              2.d0 * (dest(i,j+1) - dest(i,j))) +
     &            sig(i,j-1,2) *
     &             (dest(i+1,j-1) - dest(i+1,j) +
     &              2.d0 * (dest(i,j-1) - dest(i,j))) +
     &            sig(i,j,2) *
     &             (dest(i+1,j+1) - dest(i+1,j) +
     &              2.d0 * (dest(i,j+1) - dest(i,j)))))
 20            res(i,j) = tmp + 0.5d0 *
     &           (sig(i-1,j-1,3) *
     &             (dest(i-1,j-1) - dest(i,j)) -
     &            sig(i-1,j,3) *
     &             (dest(i-1,j+1) - dest(i,j)) -
     &            sig(i,j-1,3) *
     &             (dest(i+1,j-1) - dest(i,j)) +
     &            sig(i,j,3) *
     &             (dest(i+1,j+1) - dest(i,j)))
      end if
      end
c-----------------------------------------------------------------------
c nine-point variable stencils
      subroutine hgcen_full(cen, cenl0, cenh0, cenl1, cenh1,
     & sigx, sigy, sbl0, sbh0, sbl1, sbh1,
     & regl0, regh0, regl1, regh1,
     & hx, hy, irz, imax)
      integer cenl0, cenh0, cenl1, cenh1
      integer sbl0, sbh0, sbl1, sbh1
      integer regl0, regh0, regl1, regh1
      double precision cen(cenl0:cenh0,cenl1:cenh1)
      double precision sigx(sbl0:sbh0,sbl1:sbh1)
      double precision sigy(sbl0:sbh0,sbl1:sbh1)
      double precision hx, hy
      integer irz, imax
      double precision hxm2, hym2, facrz, r0, r1
      integer i, j
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      if (irz .eq. 0) then
         do 10 j = regl1, regh1
            do 10 i = regl0, regh0
 10            cen(i,j) = 3.d0 / (hxm2 * (sigx(i-1,j-1) + sigx(i-1,j) +
     &                                    sigx(i,j-1) + sigx(i,j)) +
     &                            hym2 * (sigy(i-1,j-1) + sigy(i-1,j) +
     &                                    sigy(i,j-1) + sigy(i,j)))
      else
         facrz = 0.25d0 * hx * hym2
         do 20 i = regl0, regh0
            r1 = (i + 0.5d0) * hx
            r0 = r1 - hx
            if (i .eq. imax) then
               r1 = -r0
            end if
            do 20 j = regl1, regh1
 20            cen(i,j) = 3.d0 / (hxm2 * (sigx(i-1,j-1) + sigx(i-1,j) +
     &                                    sigx(i,j-1) + sigx(i,j)) +
     &                            hym2 * (sigy(i-1,j-1) + sigy(i-1,j) +
     &                                    sigy(i,j-1) + sigy(i,j)) +
     &                    facrz * ((sigy(i-1,j-1) + sigy(i-1,j)) / r0 -
     &                             (sigy(i,j-1)   + sigy(i,j))   / r1))
      end if
      end
c-----------------------------------------------------------------------
      subroutine hgrlx_full_stencil(cor, corl0, corh0, corl1, corh1,
     & res, resl0, resh0, resl1, resh1,
     & sigx, sigy, sfl0, sfh0, sfl1, sfh1,
     & cen, cenl0, cenh0, cenl1, cenh1,
     & regl0, regh0, regl1, regh1,
     & hx, hy, irz, imax)
      integer corl0, corh0, corl1, corh1
      integer resl0, resh0, resl1, resh1
      integer sfl0, sfh0, sfl1, sfh1
      integer cenl0, cenh0, cenl1, cenh1
      integer regl0, regh0, regl1, regh1
      double precision cor(corl0:corh0,corl1:corh1)
      double precision res(resl0:resh0,resl1:resh1)
      double precision sigx(sfl0:sfh0,sfl1:sfh1)
      double precision sigy(sfl0:sfh0,sfl1:sfh1)
      double precision cen(cenl0:cenh0,cenl1:cenh1)
      double precision hx, hy
      integer irz, imax
      double precision hxm2, hym2, fac, facrz, r0, r1
      integer i, j
      double precision AVG, AVGRZ
      AVG() = fac * (hxm2 *
     &          (sigx(i-1,j-1) *
     &            (cor(i-1,j-1) - cor(i,j-1) + 2.d0 * cor(i-1,j)) +
     &           sigx(i-1,j) *
     &            (cor(i-1,j+1) - cor(i,j+1) + 2.d0 * cor(i-1,j)) +
     &           sigx(i,j-1) *
     &            (cor(i+1,j-1) - cor(i,j-1) + 2.d0 * cor(i+1,j)) +
     &           sigx(i,j) *
     &            (cor(i+1,j+1) - cor(i,j+1) + 2.d0 * cor(i+1,j))) +
     &               hym2 *
     &          (sigy(i-1,j-1) *
     &            (cor(i-1,j-1) - cor(i-1,j) + 2.d0 * cor(i,j-1)) +
     &           sigy(i-1,j) *
     &            (cor(i-1,j+1) - cor(i-1,j) + 2.d0 * cor(i,j+1)) +
     &           sigy(i,j-1) *
     &            (cor(i+1,j-1) - cor(i+1,j) + 2.d0 * cor(i,j-1)) +
     &           sigy(i,j) *
     &            (cor(i+1,j+1) - cor(i+1,j) + 2.d0 * cor(i,j+1))))
      AVGRZ() = AVG() + facrz *
     &           ((sigy(i-1,j-1) / r0 - sigy(i,j-1) / r1) * cor(i,j-1) +
     &             (sigy(i-1,j)   / r0 - sigy(i,j)   / r1) * cor(i,j+1))
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      fac = 1.d0 / 6.d0
      if (irz .eq. 0) then
         if (regh1 - regl1 .gt. regh0 - regl0) then
            do 10 i = regl0, regh0
cdir$ ivdep
               do 10 j = regl1, regh1
                  cor(i,j) = (AVG() - res(i,j)) * cen(i,j)
 10               continue
         else
            do 20 j = regl1, regh1
cdir$ ivdep
               do 20 i = regl0, regh0
                  cor(i,j) = (AVG() - res(i,j)) * cen(i,j)
 20               continue
         end if
      else
         facrz = hx * hym2 / 12.d0
         if (regh1 - regl1 .gt. regh0 - regl0) then
            do 30 i = regl0, regh0
               r1 = (i + 0.5d0) * hx
               r0 = r1 - hx
               if (i .eq. imax) then
                  r1 = -r0
               end if
cdir$ ivdep
               do 30 j = regl1, regh1
                  cor(i,j) = (AVGRZ() - res(i,j)) * cen(i,j)
 30               continue
         else
            do 40 j = regl1, regh1
cdir$ ivdep
               do 40 i = regl0, regh0
                  r1 = (i + 0.5d0) * hx
                  r0 = r1 - hx
                  if (i .eq. imax) then
                     r1 = -r0
                  end if
                  cor(i,j) = (AVGRZ() - res(i,j)) * cen(i,j)
 40               continue
         end if
      end if
      end
c-----------------------------------------------------------------------
      subroutine hgrlxl_full_stencil(cor, corl0, corh0, corl1, corh1,
     & res, resl0, resh0, resl1, resh1,
     & sigx, sigy, sfl0, sfh0, sfl1, sfh1,
     & cen, cenl0, cenh0, cenl1, cenh1,
     & regl0, regh0, regl1, regh1,
     & doml0, domh0, doml1, domh1,
     & hx, hy, irz, imax, lsd)
      integer corl0, corh0, corl1, corh1
      integer resl0, resh0, resl1, resh1
      integer sfl0, sfh0, sfl1, sfh1
      integer cenl0, cenh0, cenl1, cenh1
      integer regl0, regh0, regl1, regh1
      integer doml0, domh0, doml1, domh1
      double precision cor(corl0:corh0,corl1:corh1)
      double precision res(resl0:resh0,resl1:resh1)
      double precision sigx(sfl0:sfh0,sfl1:sfh1)
      double precision sigy(sfl0:sfh0,sfl1:sfh1)
      double precision cen(cenl0:cenh0,cenl1:cenh1)
      double precision hx, hy
      double precision wrk(256)
      integer irz, imax, lsd
      double precision hxm2, hym2, fac, facrz, r0, r1, betm, aj
      integer i, j, jw, ipass
      double precision RHSL0, RHSL1, RHSRZL0, RHSRZL1
      double precision AVG, AVGRZ
      RHSL0() = (res(i,j) - fac * (hxm2 *
     &           (sigx(i-1,j-1) *
     &             (cor(i-1,j-1) - cor(i,j-1)) +
     &            sigx(i-1,j) *
     &             (cor(i-1,j+1) - cor(i,j+1)) +
     &            sigx(i,j-1) *
     &             (cor(i+1,j-1) - cor(i,j-1)) +
     &            sigx(i,j) * 
     &             (cor(i+1,j+1) - cor(i,j+1))) + 
     &                            hym2 * 
     &           (sigy(i-1,j-1) * 
     &             (cor(i-1,j-1) + 2.d0 * cor(i,j-1)) + 
     &            sigy(i-1,j) * 
     &             (cor(i-1,j+1) + 2.d0 * cor(i,j+1)) + 
     &            sigy(i,j-1) * 
     &             (cor(i+1,j-1) + 2.d0 * cor(i,j-1)) + 
     &            sigy(i,j) * 
     &             (cor(i+1,j+1) + 2.d0 * cor(i,j+1)))))
      RHSL1() = (res(i,j) - fac * (hxm2 *
     &          (sigx(i-1,j-1) *
     &             (cor(i-1,j-1) + 2.d0 * cor(i-1,j)) +
     &            sigx(i-1,j) *
     &             (cor(i-1,j+1) + 2.d0 * cor(i-1,j)) +
     &            sigx(i,j-1) *
     &             (cor(i+1,j-1) + 2.d0 * cor(i+1,j)) +
     &            sigx(i,j) *
     &             (cor(i+1,j+1) + 2.d0 * cor(i+1,j))) +
     &                            hym2 * 
     &           (sigy(i-1,j-1) * 
     &             (cor(i-1,j-1) - cor(i-1,j)) + 
     &            sigy(i-1,j) * 
     &             (cor(i-1,j+1) - cor(i-1,j)) + 
     &            sigy(i,j-1) * 
     &             (cor(i+1,j-1) - cor(i+1,j)) + 
     &            sigy(i,j) * 
     &             (cor(i+1,j+1) - cor(i+1,j)))))

      RHSRZL0() = (RHSL0() - facrz *
     &           ((sigy(i-1,j-1) / r0 - sigy(i,j-1) / r1) * cor(i,j-1) +
     &            (sigy(i-1,j)   / r0 - sigy(i,j)   / r1) * cor(i,j+1)))

      RHSRZL1()= RHSL1()
      AVG() = fac * (hxm2 *
     &          (sigx(i-1,j-1) *
     &            (cor(i-1,j-1) - cor(i,j-1) + 2.d0 * cor(i-1,j)) +
     &           sigx(i-1,j) *
     &            (cor(i-1,j+1) - cor(i,j+1) + 2.d0 * cor(i-1,j)) +
     &           sigx(i,j-1) *
     &            (cor(i+1,j-1) - cor(i,j-1) + 2.d0 * cor(i+1,j)) +
     &           sigx(i,j) *
     &            (cor(i+1,j+1) - cor(i,j+1) + 2.d0 * cor(i+1,j))) +
     &               hym2 *
     &          (sigy(i-1,j-1) *
     &            (cor(i-1,j-1) - cor(i-1,j) + 2.d0 * cor(i,j-1)) +
     &           sigy(i-1,j) *
     &            (cor(i-1,j+1) - cor(i-1,j) + 2.d0 * cor(i,j+1)) +
     &           sigy(i,j-1) *
     &            (cor(i+1,j-1) - cor(i+1,j) + 2.d0 * cor(i,j-1)) +
     &           sigy(i,j) *
     &            (cor(i+1,j+1) - cor(i+1,j) + 2.d0 * cor(i,j+1))))
      AVGRZ() = AVG() + facrz *
     &           ((sigy(i-1,j-1) / r0 - sigy(i,j-1) / r1) * cor(i,j-1) +
     &             (sigy(i-1,j)   / r0 - sigy(i,j)   / r1) * cor(i,j+1))
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
                  cor(i,j) = RHSL1() * betm
                  aj = fac *
     &              (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     &               hxm2 * (sigx(i-1,j) + sigx(i,j)))
                  wrk(1) = 2.d0 * aj * betm
               else
c interface to grid at same level:
                  aj = fac *
     &              (hym2 * 2.d0 * (sigy(i-1,j-1) + sigy(i,j-1)) -
     &               hxm2 * (sigx(i-1,j-1) + sigx(i,j-1)))
                  cor(i,j) = (RHSL1() - aj * cor(i,j-1)) * betm
                  aj = fac *
     &              (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     &               hxm2 * (sigx(i-1,j) + sigx(i,j)))
                  wrk(1) = aj * betm
               end if
c forward solve loop:
               do 10 j = regl1 + 1, regh1 - 1
                  jw = j - regl1
                  if (cen(i,j) .eq. 0.d0) then
                     betm = 0.d0
                  else
                     betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(jw))
                  end if
                  cor(i,j) = (RHSL1() - aj * cor(i,j-1)) * betm
                  aj = fac *
     &              (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     &               hxm2 * (sigx(i-1,j) + sigx(i,j)))
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
                  cor(i,j) = (RHSL1() - aj * cor(i,j-1)) * betm
               else if (jw .gt. 0) then
c interface to grid at same level:
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(jw))
                  cor(i,j) = RHSL1() - aj * cor(i,j-1)
                  aj = fac *
     &              (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     &               hxm2 * (sigx(i-1,j) + sigx(i,j)))
                  cor(i,j) = (cor(i,j) - aj * cor(i,j+1)) * betm
               end if
c back substitution loop:
               do 20 j = regh1 - 1, regl1, -1
                  jw = j - regl1
                  cor(i,j) = cor(i,j) - wrk(jw + 1) * cor(i,j+1)
 20               continue
c update neumann bdys:
               if (regl1 .eq. doml1) then
                  cor(i,regl1-1) = cor(i,regl1+1)
               end if
               if (regh1 .eq. domh1) then
                  cor(i,regh1+1) = cor(i,regh1-1)
               end if
 30            continue
         else
            print *, "Line solve not implemented in dimension", lsd
		stop
         end if
      else
         print *, "Line solve not implemented for rz:"
         print *, "   Using Gauss-Seidel instead."
	    stop
         facrz = hx * hym2 / 12.d0
         if (regh1 - regl1 .gt. regh0 - regl0) then
            do 40 i = regl0, regh0
               r1 = (i + 0.5d0) * hx
               r0 = r1 - hx
               if (i .eq. imax) then
                  r1 = -r0
               end if
cdir$ ivdep
               do 40 j = regl1, regh1
                  cor(i,j) = (AVGRZ() - res(i,j)) * cen(i,j)
 40               continue
         else
            do 50 j = regl1, regh1
cdir$ ivdep
               do 50 i = regl0, regh0
                  r1 = (i + 0.5d0) * hx
                  r0 = r1 - hx
                  if (i .eq. imax) then
                     r1 = -r0
                  end if
                  cor(i,j) = (AVGRZ() - res(i,j)) * cen(i,j)
 50               continue
         end if
      end if
      end
c-----------------------------------------------------------------------
      subroutine hgrlnf_full_stencil(cor, corl0, corh0, corl1, corh1,
     & res, resl0, resh0, resl1, resh1,
     & wrk, wrkl0, wrkh0, wrkl1, wrkh1,
     & sigx, sigy, sfl0, sfh0, sfl1, sfh1,
     & cen, cenl0, cenh0, cenl1, cenh1,
     & regl0, regh0, regl1, regh1,
     & doml0, domh0, doml1, domh1,
     & hx, hy, irz, imax, lsd, ipass)
      integer corl0, corh0, corl1, corh1
      integer resl0, resh0, resl1, resh1
      integer wrkl0, wrkh0, wrkl1, wrkh1
      integer sfl0, sfh0, sfl1, sfh1
      integer cenl0, cenh0, cenl1, cenh1
      integer regl0, regh0, regl1, regh1
      integer doml0, domh0, doml1, domh1
      double precision cor(corl0:corh0,corl1:corh1)
      double precision res(resl0:resh0,resl1:resh1)
      double precision wrk(wrkl0:wrkh0,wrkl1:wrkh1)
      double precision sigx(sfl0:sfh0,sfl1:sfh1)
      double precision sigy(sfl0:sfh0,sfl1:sfh1)
      double precision cen(cenl0:cenh0,cenl1:cenh1)
      double precision hx, hy
      integer irz, imax, lsd, ipass
      double precision hxm2, hym2, fac, facrz, r0, r1, betm, aj
      integer i, j, ioff
      double precision RHSL0, RHSL1, RHSRZL0, RHSRZL1
      double precision AVG, AVGRZ
      RHSL0() = (res(i,j) - fac * (hxm2 *
     &           (sigx(i-1,j-1) *
     &             (cor(i-1,j-1) - cor(i,j-1)) +
     &            sigx(i-1,j) *
     &             (cor(i-1,j+1) - cor(i,j+1)) +
     &            sigx(i,j-1) *
     &             (cor(i+1,j-1) - cor(i,j-1)) +
     &            sigx(i,j) * 
     &             (cor(i+1,j+1) - cor(i,j+1))) + 
     &                            hym2 * 
     &           (sigy(i-1,j-1) * 
     &             (cor(i-1,j-1) + 2.d0 * cor(i,j-1)) + 
     &            sigy(i-1,j) * 
     &             (cor(i-1,j+1) + 2.d0 * cor(i,j+1)) + 
     &            sigy(i,j-1) * 
     &             (cor(i+1,j-1) + 2.d0 * cor(i,j-1)) + 
     &            sigy(i,j) * 
     &             (cor(i+1,j+1) + 2.d0 * cor(i,j+1)))))
      RHSL1() = (res(i,j) - fac * (hxm2 *
     &          (sigx(i-1,j-1) *
     &             (cor(i-1,j-1) + 2.d0 * cor(i-1,j)) +
     &            sigx(i-1,j) *
     &             (cor(i-1,j+1) + 2.d0 * cor(i-1,j)) +
     &            sigx(i,j-1) *
     &             (cor(i+1,j-1) + 2.d0 * cor(i+1,j)) +
     &            sigx(i,j) *
     &             (cor(i+1,j+1) + 2.d0 * cor(i+1,j))) +
     &                            hym2 * 
     &           (sigy(i-1,j-1) * 
     &             (cor(i-1,j-1) - cor(i-1,j)) + 
     &            sigy(i-1,j) * 
     &             (cor(i-1,j+1) - cor(i-1,j)) + 
     &            sigy(i,j-1) * 
     &             (cor(i+1,j-1) - cor(i+1,j)) + 
     &            sigy(i,j) * 
     &             (cor(i+1,j+1) - cor(i+1,j)))))

      RHSRZL0() = (RHSL0() - facrz *
     &           ((sigy(i-1,j-1) / r0 - sigy(i,j-1) / r1) * cor(i,j-1) +
     &            (sigy(i-1,j)   / r0 - sigy(i,j)   / r1) * cor(i,j+1)))

      RHSRZL1()= RHSL1()
      AVG() = fac * (hxm2 *
     &          (sigx(i-1,j-1) *
     &            (cor(i-1,j-1) - cor(i,j-1) + 2.d0 * cor(i-1,j)) +
     &           sigx(i-1,j) *
     &            (cor(i-1,j+1) - cor(i,j+1) + 2.d0 * cor(i-1,j)) +
     &           sigx(i,j-1) *
     &            (cor(i+1,j-1) - cor(i,j-1) + 2.d0 * cor(i+1,j)) +
     &           sigx(i,j) *
     &            (cor(i+1,j+1) - cor(i,j+1) + 2.d0 * cor(i+1,j))) +
     &               hym2 *
     &          (sigy(i-1,j-1) *
     &            (cor(i-1,j-1) - cor(i-1,j) + 2.d0 * cor(i,j-1)) +
     &           sigy(i-1,j) *
     &            (cor(i-1,j+1) - cor(i-1,j) + 2.d0 * cor(i,j+1)) +
     &           sigy(i,j-1) *
     &            (cor(i+1,j-1) - cor(i+1,j) + 2.d0 * cor(i,j-1)) +
     &           sigy(i,j) *
     &            (cor(i+1,j+1) - cor(i+1,j) + 2.d0 * cor(i,j+1))))
      AVGRZ() = AVG() + facrz *
     &           ((sigy(i-1,j-1) / r0 - sigy(i,j-1) / r1) * cor(i,j-1) +
     &             (sigy(i-1,j)   / r0 - sigy(i,j)   / r1) * cor(i,j+1))
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      fac = 1.d0 / 6.d0
      if (irz .eq. 0) then
         if (lsd .eq. 0) then
            if (mod(regl1,2) .eq. 0) then
               ioff = ipass
            else
               ioff = 1 - ipass
            end if
            i = regl0
            do 10 j = regl1 + ioff, regh1, 2
               aj = fac *
     &           (hxm2 * 2.d0 * (sigx(i,j-1) + sigx(i,j)) -
     &            hym2 * (sigy(i,j-1) + sigy(i,j)))
               if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
                  wrk(i,j) = 0.d0
               else if (regl0 .eq. doml0) then
c neumann bdy:
                  betm = -cen(i,j)
                  cor(i,j) = RHSL0() * betm
                  wrk(i,j) = 2.d0 * aj * betm
               end if
               wrk(i+1,j) = aj
 10            continue
c forward solve loop:
            do 20 i = regl0 + 1, regh0 - 1
               do 20 j = regl1 + ioff, regh1, 2
                  aj = wrk(i,j)
                  if (cen(i,j) .eq. 0.d0) then
                     betm = 0.d0
                  else
                     betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i-1,j))
                  end if
                  cor(i,j) = (RHSL0() - aj * cor(i-1,j)) * betm
                  aj = fac *
     &              (hxm2 * 2.d0 * (sigx(i,j-1) + sigx(i,j)) -
     &               hym2 * (sigy(i,j-1) + sigy(i,j)))
                  wrk(i+1,j) = aj
 20               wrk(i,j) = aj * betm
            i = regh0
            do 30 j = regl1 + ioff, regh1, 2
               if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
               else if (regh0 .eq. domh0) then
c neumann bdy:
                  aj = 2.d0 * wrk(i,j)
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i-1,j))
                  cor(i,j) = (RHSL0() - aj * cor(i-1,j)) * betm
               else if (i .gt. regl0) then
c interface to grid at same level:
                  aj = wrk(i,j)
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i-1,j))
                  cor(i,j) = (RHSL0() - aj * cor(i-1,j)) * betm
                  aj = fac *
     &              (hxm2 * 2.d0 * (sigx(i,j-1) + sigx(i,j)) -
     &               hym2 * (sigy(i,j-1) + sigy(i,j)))
                  wrk(i,j) = aj * betm
               end if
 30            continue
         else
            if (mod(regl0,2) .eq. 0) then
               ioff = ipass
            else
               ioff = 1 - ipass
            end if
            j = regl1
            do 110 i = regl0 + ioff, regh0, 2
               aj = fac *
     &           (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     &            hxm2 * (sigx(i-1,j) + sigx(i,j)))
               if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
                  wrk(i,j) = 0.d0
               else if (regl1 .eq. doml1) then
c neumann bdy:
                  betm = -cen(i,j)
                  cor(i,j) = RHSL1() * betm
                  wrk(i,j) = 2.d0 * aj * betm
               end if
               wrk(i,j+1) = aj
 110           continue
c forward solve loop:
            do 120 j = regl1 + 1, regh1 - 1
               do 120 i = regl0 + ioff, regh0, 2
                  aj = wrk(i,j)
                  if (cen(i,j) .eq. 0.d0) then
                     betm = 0.d0
                  else
                     betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
                  end if
                  cor(i,j) = (RHSL1() - aj * cor(i,j-1)) * betm
                  aj = fac *
     &              (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     &               hxm2 * (sigx(i-1,j) + sigx(i,j)))
                  wrk(i,j+1) = aj
 120              wrk(i,j) = aj * betm
            j = regh1
            do 130 i = regl0 + ioff, regh0, 2
               if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
               else if (regh1 .eq. domh1) then
c neumann bdy:
                  aj = 2.d0 * wrk(i,j)
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
                  cor(i,j) = (RHSL1() - aj * cor(i,j-1)) * betm
               else if (j .gt. regl1) then
c interface to grid at same level:
                  aj = wrk(i,j)
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
                  cor(i,j) = (RHSL1() - aj * cor(i,j-1)) * betm
                  aj = fac *
     &              (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     &               hxm2 * (sigx(i-1,j) + sigx(i,j)))
                  wrk(i,j) = aj * betm
               end if
 130           continue
         end if
      else
         facrz = hx * hym2 / 12.d0
         if (lsd .eq. 0) then
            if (mod(regl1,2) .eq. 0) then
               ioff = ipass
            else
               ioff = 1 - ipass
            end if
            i = regl0
            r1 = (i + 0.5d0) * hx
            r0 = r1 - hx
            do 210 j = regl1 + ioff, regh1, 2
               aj = fac *
     &           (hxm2 * 2.d0 * (sigx(i,j-1) + sigx(i,j)) -
     &            hym2 * (sigy(i,j-1) + sigy(i,j)))
               if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
                  wrk(i,j) = 0.d0
               else if (regl0 .eq. doml0) then
c neumann bdy:
                  betm = -cen(i,j)
                  cor(i,j) = RHSRZL0() * betm
                  wrk(i,j) = 2.d0 * aj * betm
               end if
               wrk(i+1,j) = aj
 210           continue
c forward solve loop:
            do 220 i = regl0 + 1, regh0 - 1
               r1 = (i + 0.5d0) * hx
               r0 = r1 - hx
               do 220 j = regl1 + ioff, regh1, 2
                  aj = wrk(i,j)
                  if (cen(i,j) .eq. 0.d0) then
                     betm = 0.d0
                  else
                     betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i-1,j))
                  end if
                  cor(i,j) = (RHSRZL0() - aj * cor(i-1,j)) * betm
                  aj = fac *
     &              (hxm2 * 2.d0 * (sigx(i,j-1) + sigx(i,j)) -
     &               hym2 * (sigy(i,j-1) + sigy(i,j)))
                  wrk(i+1,j) = aj
 220              wrk(i,j) = aj * betm
            i = regh0
            r1 = (i + 0.5d0) * hx
            r0 = r1 - hx
            if (i .eq. imax) then
               r1 = -r0
            end if
            do 230 j = regl1 + ioff, regh1, 2
               if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
               else if (regh0 .eq. domh0) then
c neumann bdy:
                  aj = 2.d0 * wrk(i,j)
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i-1,j))
                  cor(i,j) = (RHSRZL0() - aj * cor(i-1,j)) * betm
               else if (i .gt. regl0) then
c interface to grid at same level:
                  aj = wrk(i,j)
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i-1,j))
                  cor(i,j) = (RHSRZL0() - aj * cor(i-1,j)) * betm
                  aj = fac *
     &              (hxm2 * 2.d0 * (sigx(i,j-1) + sigx(i,j)) -
     &               hym2 * (sigy(i,j-1) + sigy(i,j)))
                  wrk(i,j) = aj * betm
               end if
 230           continue
         else
            if (mod(regl0,2) .eq. 0) then
               ioff = ipass
            else
               ioff = 1 - ipass
            end if
            j = regl1
            do 310 i = regl0 + ioff, regh0, 2
               r1 = (i + 0.5d0) * hx
               r0 = r1 - hx
               if (i .eq. imax) then
                  r1 = -r0
               end if
               aj = fac *
     &                (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     &                 hxm2 * (sigx(i-1,j) + sigx(i,j))) +
     &              facrz * (sigy(i-1,j) / r0 - sigy(i,j) / r1)
               if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
                  wrk(i,j) = 0.d0
               else if (regl1 .eq. doml1) then
c neumann bdy:
                  betm = -cen(i,j)
                  cor(i,j) = RHSRZL1() * betm
                  wrk(i,j) = 2.d0 * aj * betm
               end if
               wrk(i,j+1) = aj
 310           continue
c forward solve loop:
            do 320 j = regl1 + 1, regh1 - 1
               do 320 i = regl0 + ioff, regh0, 2
                  r1 = (i + 0.5d0) * hx
                  r0 = r1 - hx
                  if (i .eq. imax) then
                     r1 = -r0
                  end if
                  aj = wrk(i,j)
                  if (cen(i,j) .eq. 0.d0) then
                     betm = 0.d0
                  else
                     betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
                  end if
                  cor(i,j) = (RHSRZL1() - aj * cor(i,j-1)) * betm
                  aj = fac *
     &                   (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     &                    hxm2 * (sigx(i-1,j) + sigx(i,j))) +
     &                 facrz * (sigy(i-1,j) / r0 - sigy(i,j) / r1)
                  wrk(i,j+1) = aj
 320              wrk(i,j) = aj * betm
            j = regh1
            do 330 i = regl0 + ioff, regh0, 2
               if (cen(i,j) .eq. 0.d0) then
c dirichlet bdy:
               else if (regh1 .eq. domh1) then
c neumann bdy:
                  aj = 2.d0 * wrk(i,j)
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
                  cor(i,j) = (RHSRZL1() - aj * cor(i,j-1)) * betm
               else if (j .gt. regl1) then
c interface to grid at same level:
                  r1 = (i + 0.5d0) * hx
                  r0 = r1 - hx
                  if (i .eq. imax) then
                     r1 = -r0
                  end if
                  aj = wrk(i,j)
                  betm = 1.d0 / (-1.d0 / cen(i,j) - aj * wrk(i,j-1))
                  cor(i,j) = (RHSRZL1() - aj * cor(i,j-1)) * betm
                  aj = fac *
     &                   (hym2 * 2.d0 * (sigy(i-1,j) + sigy(i,j)) -
     &                    hxm2 * (sigx(i-1,j) + sigx(i,j))) +
     &                 facrz * (sigy(i-1,j) / r0 - sigy(i,j) / r1)
                  wrk(i,j) = aj * betm
               end if
 330           continue
         end if
      end if
      end
c-----------------------------------------------------------------------
      subroutine hgres_full_stencil(res, resl0, resh0, resl1, resh1,
     & src,  srcl0, srch0, srcl1, srch1,
     & dest, destl0, desth0, destl1, desth1,
     & sigx, sigy, sfl0, sfh0, sfl1, sfh1,
     & regl0, regh0, regl1, regh1,
     & hx, hy, irz, imax)
      integer resl0, resh0, resl1, resh1
      integer srcl0, srch0, srcl1, srch1
      integer destl0, desth0, destl1, desth1
      integer sfl0, sfh0, sfl1, sfh1
      integer regl0, regh0, regl1, regh1
      double precision res(resl0:resh0,resl1:resh1)
      double precision src(srcl0:srch0,srcl1:srch1)
      double precision dest(destl0:desth0,destl1:desth1)
      double precision sigx(sfl0:sfh0,sfl1:sfh1)
      double precision sigy(sfl0:sfh0,sfl1:sfh1)
      double precision hx, hy
      integer irz, imax
      double precision hxm2, hym2, fac, r0, r1, tmp
      integer i, j
      hxm2 = 1.d0 / (hx*hx)
      hym2 = 1.d0 / (hy*hy)
      fac = 1.d0 / 6.d0
      if (regh1 - regl1 .gt. regh0 - regl0) then
         do 10 i = regl0, regh0
            do 10 j = regl1, regh1
               tmp = hxm2 *
     &           (sigx(i-1,j-1) *
     &             (dest(i-1,j-1) - dest(i,j-1) +
     &              2.d0 * (dest(i-1,j) - dest(i,j))) +
     &            sigx(i-1,j) *
     &             (dest(i-1,j+1) - dest(i,j+1) +
     &              2.d0 * (dest(i-1,j) - dest(i,j))) +
     &            sigx(i,j-1) *
     &             (dest(i+1,j-1) - dest(i,j-1) +
     &              2.d0 * (dest(i+1,j) - dest(i,j))) +
     &            sigx(i,j) *
     &             (dest(i+1,j+1) - dest(i,j+1) +
     &              2.d0 * (dest(i+1,j) - dest(i,j))))
 10            res(i,j) = src(i,j) - fac * (tmp + hym2 *
     &           (sigy(i-1,j-1) *
     &             (dest(i-1,j-1) - dest(i-1,j) +
     &              2.d0 * (dest(i,j-1) - dest(i,j))) +
     &            sigy(i-1,j) *
     &             (dest(i-1,j+1) - dest(i-1,j) +
     &              2.d0 * (dest(i,j+1) - dest(i,j))) +
     &            sigy(i,j-1) *
     &             (dest(i+1,j-1) - dest(i+1,j) +
     &              2.d0 * (dest(i,j-1) - dest(i,j))) +
     &            sigy(i,j) *
     &             (dest(i+1,j+1) - dest(i+1,j) +
     &              2.d0 * (dest(i,j+1) - dest(i,j)))))
      else
         do 20 j = regl1, regh1
            do 20 i = regl0, regh0
               tmp = hxm2 *
     &           (sigx(i-1,j-1) *
     &             (dest(i-1,j-1) - dest(i,j-1) +
     &              2.d0 * (dest(i-1,j) - dest(i,j))) +
     &            sigx(i-1,j) *
     &             (dest(i-1,j+1) - dest(i,j+1) +
     &              2.d0 * (dest(i-1,j) - dest(i,j))) +
     &            sigx(i,j-1) *
     &             (dest(i+1,j-1) - dest(i,j-1) +
     &              2.d0 * (dest(i+1,j) - dest(i,j))) +
     &            sigx(i,j) *
     &             (dest(i+1,j+1) - dest(i,j+1) +
     &              2.d0 * (dest(i+1,j) - dest(i,j))))
 20            res(i,j) = src(i,j) - fac * (tmp + hym2 *
     &           (sigy(i-1,j-1) *
     &             (dest(i-1,j-1) - dest(i-1,j) +
     &              2.d0 * (dest(i,j-1) - dest(i,j))) +
     &            sigy(i-1,j) *
     &             (dest(i-1,j+1) - dest(i-1,j) +
     &              2.d0 * (dest(i,j+1) - dest(i,j))) +
     &            sigy(i,j-1) *
     &             (dest(i+1,j-1) - dest(i+1,j) +
     &              2.d0 * (dest(i,j-1) - dest(i,j))) +
     &            sigy(i,j) *
     &             (dest(i+1,j+1) - dest(i+1,j) +
     &              2.d0 * (dest(i,j+1) - dest(i,j)))))
      end if
      if (irz .eq. 1) then
         fac = hx * hym2 / 12.d0
         do 30 i = regl0, regh0
            r1 = (i + 0.5d0) * hx
            r0 = r1 - hx
            if (i .eq. imax) then
               r1 = -r0
            end if
            do 30 j = regl1, regh1
 30            res(i,j) = res(i,j) - fac *
     &          ((sigy(i-1,j-1) * (dest(i,j-1) - dest(i,j)) +
     &            sigy(i-1,j)   * (dest(i,j+1) - dest(i,j))) / r0 -
     &           (sigy(i,j-1)   * (dest(i,j-1) - dest(i,j)) +
     &            sigy(i,j)     * (dest(i,j+1) - dest(i,j))) / r1)
      end if
      end

c-----------------------------------------------------------------------
      subroutine hgrlnb(cor, corl0, corh0, corl1, corh1,
     & wrk, wrkl0, wrkh0, wrkl1, wrkh1,
     & regl0, regh0, regl1, regh1,
     & lsd, ipass)
      integer corl0, corh0, corl1, corh1
      integer wrkl0, wrkh0, wrkl1, wrkh1
      integer regl0, regh0, regl1, regh1
      double precision cor(corl0:corh0,corl1:corh1)
      double precision wrk(wrkl0:wrkh0,wrkl1:wrkh1)
      integer lsd, ipass
      integer i, j, ioff
      if (lsd .eq. 0) then
         if (mod(regl1,2) .eq. 0) then
            ioff = ipass
         else
            ioff = 1 - ipass
         end if
c back substitution loop:
         do 10 i = regh0 - 1, regl0, -1
            do 10 j = regl1 + ioff, regh1, 2
               cor(i,j) = cor(i,j) - wrk(i,j) * cor(i+1,j)
 10            continue
      else
         if (mod(regl0,2) .eq. 0) then
            ioff = ipass
         else
            ioff = 1 - ipass
         end if
c back substitution loop:
         do 20 j = regh1 - 1, regl1, -1
            do 20 i = regl0 + ioff, regh0, 2
               cor(i,j) = cor(i,j) - wrk(i,j) * cor(i,j+1)
 20            continue
      end if
      end

c Unrolled indexing in these 3 routines uses the fact that each array
c has a border of width 1
c-----------------------------------------------------------------------
c Works for NODE-based data.
      subroutine hgip(v0, v1, mask,
     & regl0, regh0, regl1, regh1, sum)
      integer regl0, regh0, regl1, regh1
      double precision v0(*)
      double precision v1(*)
      double precision mask(*)
      double precision sum
      integer i, idiff
c      do 10 i = 1, (regh0 - regl0 + 1) * (regh1 - regl1 + 1)
      idiff = regh0 - regl0 + 1
      do 10 i = idiff + 2, idiff * (regh1 - regl1) - 1
 10      sum = sum + mask(i) * v0(i) * v1(i)
      end
c-----------------------------------------------------------------------
      subroutine hgcg1(r, p, z, x, w, c, mask,
     & regl0, regh0, regl1, regh1,
     & alpha, rho)
      integer regl0, regh0, regl1, regh1
      double precision r(*)
      double precision p(*)
      double precision z(*)
      double precision x(*)
      double precision w(*)
      double precision c(*)
      double precision mask(*)
      double precision alpha, rho
      integer i, idiff
c      do 10 i = 1, (regh0 - regl0 + 1) * (regh1 - regl1 + 1)
      idiff = regh0 - regl0 + 1
      do 10 i = idiff + 2, idiff * (regh1 - regl1) - 1
         r(i) = r(i) - alpha * w(i)
         x(i) = x(i) + alpha * p(i)
         z(i) = r(i) * c(i)
 10      rho = rho + mask(i) * z(i) * r(i)
      end

      subroutine hgcg2(p, z,
     & regl0, regh0, regl1, regh1,
     & alpha)
      integer regl0, regh0, regl1, regh1
      double precision p(*)
      double precision z(*)
      double precision alpha
      integer i, idiff
c      do 10 i = 1, (regh0 - regl0 + 1) * (regh1 - regl1 + 1)
      idiff = regh0 - regl0 + 1
      do 10 i = idiff + 2, idiff * (regh1 - regl1) - 1
 10      p(i) = alpha * p(i) + z(i)
      end
