c-----------------------------------------------------------------------
      subroutine hgfres_terrain(
     & res,   resl0, resh0, resl1, resh1,
     & src,   srcl0, srch0, srcl1, srch1,
     & fdst,  fdstl0, fdsth0, fdstl1, fdsth1,
     & cdst,  cdstl0, cdsth0, cdstl1, cdsth1,
     & sigmaf, sfl0, sfh0, sfl1, sfh1,
     & sigmac, scl0, sch0, scl1, sch1,
     &        regl0, regh0, regl1, regh1,
     & hx, hy, ir, jr, idim, idir, idd1)
      integer resl0, resh0, resl1, resh1
      integer srcl0, srch0, srcl1, srch1
      integer fdstl0, fdsth0, fdstl1, fdsth1
      integer cdstl0, cdsth0, cdstl1, cdsth1
      integer sfl0, sfh0, sfl1, sfh1
      integer scl0, sch0, scl1, sch1
      integer regl0, regh0, regl1, regh1
      double precision hx, hy
      double precision res(resl0:resh0,resl1:resh1)
      double precision src(srcl0:srch0,srcl1:srch1)
      double precision fdst(fdstl0:fdsth0,fdstl1:fdsth1)
      double precision cdst(cdstl0:cdsth0,cdstl1:cdsth1)
      double precision sigmaf(sfl0:sfh0,sfl1:sfh1, 3)
      double precision sigmac(scl0:sch0,scl1:sch1, 3)
      integer ir, jr, idim, idir
      double precision fac0, fac1, tmp
      integer i, j, is, js, m, n
      integer idd1
      if (idim .eq. 0) then
         i = regl0
         if (idir .eq. 1) then
            is = i - 1
         else
            is = i
         end if
         fac0 = 1.d0 / 6.d0
         do j = regl1, regh1
            res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     &          (sigmac(is,j-1,1) *
     &            (cdst(i-idir,j-1) - cdst(i,j-1) +
     &             2.d0 * (cdst(i-idir,j) - cdst(i,j))) +
     &           sigmac(is,j,1) *
     &            (cdst(i-idir,j+1) - cdst(i,j+1) +
     &             2.d0 * (cdst(i-idir,j) - cdst(i,j))) +
     &           sigmac(is,j-1,2) *
     &            (cdst(i-idir,j-1) - cdst(i-idir,j) +
     &             2.d0 * (cdst(i,j-1) - cdst(i,j))) +
     &           sigmac(is,j,2) *
     &            (cdst(i-idir,j+1) - cdst(i-idir,j) +
     &             2.d0 * (cdst(i,j+1) - cdst(i,j))) -
     &           3.d0 * (sigmac(is,j-(i-is),3) *
     &                     (cdst(i-idir,j-idir) - cdst(i,j)) -
     &                   sigmac(is,j-(1-(i-is)),3) *
     &                     (cdst(i-idir,j+idir) - cdst(i,j))))
         end do
         fac0 = fac0 / jr
         i = i * ir
         if (idir .eq. 1) then
            is = i
         else
            is = i - 1
         end if
         do n = 0, jr-1
            fac1 = (jr-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do j = jr*regl1, jr*regh1, jr
               tmp =
     &            sigmaf(is,j-n-1,1) *
     &             (fdst(i+idir,j-n-1) - fdst(i,j-n-1) +
     &              2.d0 * (fdst(i+idir,j-n) - fdst(i,j-n))) +
     &              sigmaf(is,j-n,1) *
     &             (fdst(i+idir,j-n+1) - fdst(i,j-n+1) +
     &              2.d0 * (fdst(i+idir,j-n) - fdst(i,j-n))) +
     &              sigmaf(is,j+n-1,1) *
     &             (fdst(i+idir,j+n-1) - fdst(i,j+n-1) +
     &              2.d0 * (fdst(i+idir,j+n) - fdst(i,j+n))) +
     &              sigmaf(is,j+n,1) *
     &             (fdst(i+idir,j+n+1) - fdst(i,j+n+1) +
     &              2.d0 * (fdst(i+idir,j+n) - fdst(i,j+n)))
               res(i,j) = res(i,j) - fac1 * (tmp +
     &              sigmaf(is,j-n-1,2) *
     &             (fdst(i+idir,j-n-1) - fdst(i+idir,j-n) +
     &              2.d0 * (fdst(i,j-n-1) - fdst(i,j-n))) +
     &              sigmaf(is,j-n,2) *
     &             (fdst(i+idir,j-n+1) - fdst(i+idir,j-n) +
     &              2.d0 * (fdst(i,j-n+1) - fdst(i,j-n))) +
     &              sigmaf(is,j+n-1,2) *
     &             (fdst(i+idir,j+n-1) - fdst(i+idir,j+n) +
     &              2.d0 * (fdst(i,j+n-1) - fdst(i,j+n))) +
     &              sigmaf(is,j+n,2) *
     &             (fdst(i+idir,j+n+1) - fdst(i+idir,j+n) +
     &              2.d0 * (fdst(i,j+n+1) - fdst(i,j+n))) -
     &              3.d0 * (sigmaf(is,j-n-(i-is),3) *
     &             (fdst(i+idir,j-n+idir) - fdst(i+idir,j-n)) -
     &                      sigmaf(is,j-n-(1-(i-is)),3) *
     &             (fdst(i+idir,j-n-idir) - fdst(i+idir,j-n)) +
     &                      sigmaf(is,j+n-(i-is),3) *
     &             (fdst(i+idir,j+n+idir) - fdst(i+idir,j+n)) -
     &                      sigmaf(is,j+n-(1-(i-is)),3) *
     &             (fdst(i+idir,j+n-idir) - fdst(i+idir,j+n))))
            end do
         end do
      else
         j = regl1
         if (idir .eq. 1) then
            js = j - 1
         else
            js = j
         end if
         fac0 = 1.d0 / 6.d0
         do i = regl0, regh0
            res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     &          (sigmac(i-1,js,1) *
     &            (cdst(i-1,j-idir) - cdst(i,j-idir) +
     &             2.d0 * (cdst(i-1,j) - cdst(i,j))) +
     &           sigmac(i,js,1) *
     &            (cdst(i+1,j-idir) - cdst(i,j-idir) +
     &             2.d0 * (cdst(i+1,j) - cdst(i,j))) +
     &           sigmac(i-1,js,2) *
     &            (cdst(i-1,j-idir) - cdst(i-1,j) +
     &             2.d0 * (cdst(i,j-idir) - cdst(i,j))) +
     &           sigmac(i,js,2) *
     &            (cdst(i+1,j-idir) - cdst(i+1,j) +
     &             2.d0 * (cdst(i,j-idir) - cdst(i,j))) -
     &           3.d0 * (sigmac(i-(j-js),js,3) *
     &                     (cdst(i-idir,j-idir) - cdst(i,j)) -
     &                   sigmac(i-(1-(j-js)),js,3) *
     &                     (cdst(i+idir,j-idir) - cdst(i,j))))
         end do
         fac0 = fac0 / ir
         j = j * jr
         if (idir .eq. 1) then
            js = j
         else
            js = j - 1
         end if
         do m = 0, ir-1
            fac1 = (ir-m) * fac0
            if (m .eq. 0) fac1 = 0.5d0 * fac1
            do i = ir*regl0, ir*regh0, ir
               tmp =
     &              sigmaf(i-m-1,js,1) *
     &             (fdst(i-m-1,j+idir) - fdst(i-m,j+idir) +
     &              2.d0 * (fdst(i-m-1,j) - fdst(i-m,j))) +
     &              sigmaf(i-m,js,1) *
     &             (fdst(i-m+1,j+idir) - fdst(i-m,j+idir) +
     &              2.d0 * (fdst(i-m+1,j) - fdst(i-m,j))) +
     &              sigmaf(i+m-1,js,1) *
     &             (fdst(i+m-1,j+idir) - fdst(i+m,j+idir) +
     &              2.d0 * (fdst(i+m-1,j) - fdst(i+m,j))) +
     &              sigmaf(i+m,js,1) *
     &             (fdst(i+m+1,j+idir) - fdst(i+m,j+idir) +
     &              2.d0 * (fdst(i+m+1,j) - fdst(i+m,j)))
               res(i,j) = res(i,j) - fac1 * (tmp +
     &            sigmaf(i-m-1,js,2) *
     &             (fdst(i-m-1,j+idir) - fdst(i-m-1,j) +
     &              2.d0 * (fdst(i-m,j+idir) - fdst(i-m,j))) +
     &              sigmaf(i-m,js,2) *
     &             (fdst(i-m+1,j+idir) - fdst(i-m+1,j) +
     &              2.d0 * (fdst(i-m,j+idir) - fdst(i-m,j))) +
     &              sigmaf(i+m-1,js,2) *
     &             (fdst(i+m-1,j+idir) - fdst(i+m-1,j) +
     &              2.d0 * (fdst(i+m,j+idir) - fdst(i+m,j))) +
     &              sigmaf(i+m,js,2) *
     &             (fdst(i+m+1,j+idir) - fdst(i+m+1,j) +
     &              2.d0 * (fdst(i+m,j+idir) - fdst(i+m,j))) -
     &              3.d0 * (sigmaf(i-m-(j-js),js,3) *
     &             (fdst(i-m+idir,j+idir) - fdst(i-m,j+idir)) -
     &                      sigmaf(i-m-(1-(j-js)),js,3) *
     &             (fdst(i-m-idir,j+idir) - fdst(i-m,j+idir)) +
     &                      sigmaf(i+m-(j-js),js,3) *
     &             (fdst(i+m+idir,j+idir) - fdst(i+m,j+idir)) -
     &                      sigmaf(i+m-(1-(j-js)),js,3) *
     &             (fdst(i+m-idir,j+idir) - fdst(i+m,j+idir))))
            end do
         end do
      end if
      end
c-----------------------------------------------------------------------
      subroutine hgcres_terrain(
     & res,   resl0, resh0, resl1, resh1,
     & src,   srcl0, srch0, srcl1, srch1,
     & fdst,  fdstl0, fdsth0, fdstl1, fdsth1,
     & cdst,  cdstl0, cdsth0, cdstl1, cdsth1,
     & sigmaf, sfl0, sfh0, sfl1, sfh1,
     & sigmac, scl0, sch0, scl1, sch1,
     &        regl0, regh0, regl1, regh1,
     & hx, hy, ir, jr, ga, idd)
      integer resl0, resh0, resl1, resh1
      integer srcl0, srch0, srcl1, srch1
      integer fdstl0, fdsth0, fdstl1, fdsth1
      integer cdstl0, cdsth0, cdstl1, cdsth1
      integer sfl0, sfh0, sfl1, sfh1
      integer scl0, sch0, scl1, sch1
      integer regl0, regh0, regl1, regh1
      double precision hx, hy
      double precision res(resl0:resh0,resl1:resh1)
      double precision src(srcl0:srch0,srcl1:srch1)
      double precision fdst(fdstl0:fdsth0,fdstl1:fdsth1)
      double precision cdst(cdstl0:cdsth0,cdstl1:cdsth1)
      double precision sigmaf(sfl0:sfh0,sfl1:sfh1, 3)
      double precision sigmac(scl0:sch0,scl1:sch1, 3)
      integer ir, jr, ga(0:1,0:1), idd
      double precision sum, fac, fac1
      integer ic, jc, if, jf, ii, ji, idir, jdir, m, n
      ic = regl0
      jc = regl1
      if = ic * ir
      jf = jc * jr
      sum = 0.d0
c quadrants
      do ji = 0, 1
         jdir = 2 * ji - 1
         do ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,ji) .eq. 1) then
               sum = sum + sigmaf(if+ii-1,jf+ji-1,1) *
     &            (fdst(if+idir,jf+jdir) - fdst(if,jf+jdir) +
     &             2.d0 * (fdst(if+idir,jf) - fdst(if,jf))) +
     &                     sigmaf(if+ii-1,jf+ji-1,2) *
     &            (fdst(if+idir,jf+jdir) - fdst(if+idir,jf) +
     &             2.d0 * (fdst(if,jf+jdir) - fdst(if,jf))) -
     &                     sigmaf(if+ii-1,jf+ji-1,3) * idir * jdir *
     &            3.d0 * (fdst(if+idir,jf+jdir) - fdst(if,jf))
            else
               sum = sum + sigmac(ic+ii-1,jc+ji-1,1) *
     &            (cdst(ic+idir,jc+jdir) - cdst(ic,jc+jdir) +
     &             2.d0 * (cdst(ic+idir,jc) - cdst(ic,jc))) +
     &                     sigmac(ic+ii-1,jc+ji-1,2) *
     &            (cdst(ic+idir,jc+jdir) - cdst(ic+idir,jc) +
     &             2.d0 * (cdst(ic,jc+jdir) - cdst(ic,jc))) -
     &                     sigmac(ic+ii-1,jc+ji-1,3) * idir * jdir *
     &            3.d0 * (cdst(ic+idir,jc+jdir) - cdst(ic,jc))
            end if
         end do
      end do
c edges
      do ji = 0, 1
         jdir = 2 * ji - 1
         do ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,ji) - ga(ii,1-ji) .eq. 1) then
               fac1 = 1.d0 / ir
               do m = idir, idir*(ir-1), idir
                  fac = (ir-abs(m)) * fac1
                  sum = sum + fac * (sigmaf(if+m-1,jf+ji-1,1) *
     &               (fdst(if+m-1,jf+jdir) - fdst(if+m,jf+jdir) +
     &                2.d0 * (fdst(if+m-1,jf) - fdst(if+m,jf))) +
     &                               sigmaf(if+m-1,jf+ji-1,2) *
     &               (fdst(if+m-1,jf+jdir) - fdst(if+m-1,jf) +
     &                2.d0 * (fdst(if+m,jf+jdir) - fdst(if+m,jf))) +
     &                               sigmaf(if+m-1,jf+ji-1,3) * jdir *
     &               3.d0 * (fdst(if+m-1,jf+jdir) - fdst(if+m,jf)) +
     &                               sigmaf(if+m,jf+ji-1,1) *
     &               (fdst(if+m+1,jf+jdir) - fdst(if+m,jf+jdir) +
     &                2.d0 * (fdst(if+m+1,jf) - fdst(if+m,jf))) +
     &                               sigmaf(if+m,jf+ji-1,2) *
     &               (fdst(if+m+1,jf+jdir) - fdst(if+m+1,jf) +
     &                2.d0 * (fdst(if+m,jf+jdir) - fdst(if+m,jf))) -
     &                               sigmaf(if+m,jf+ji-1,3) * jdir *
     &               3.d0 * (fdst(if+m+1,jf+jdir) - fdst(if+m,jf)))
               end do
            end if
            if (ga(ii,ji) - ga(1-ii,ji) .eq. 1) then
               fac1 = 1.d0 / jr
               do n = jdir, jdir*(jr-1), jdir
                  fac = (jr-abs(n)) * fac1
                  sum = sum + fac * (sigmaf(if+ii-1,jf+n-1,1) *
     &               (fdst(if+idir,jf+n-1) - fdst(if,jf+n-1) +
     &                2.d0 * (fdst(if+idir,jf+n) - fdst(if,jf+n))) +
     &                               sigmaf(if+ii-1,jf+n-1,2) *
     &               (fdst(if+idir,jf+n-1) - fdst(if+idir,jf+n) +
     &                2.d0 * (fdst(if,jf+n-1) - fdst(if,jf+n))) +
     &                               sigmaf(if+ii-1,jf+n-1,3) * idir *
     &               3.d0 * (fdst(if+idir,jf+n-1) - fdst(if,jf+n)) +
     &                               sigmaf(if+ii-1,jf+n,1) *
     &               (fdst(if+idir,jf+n+1) - fdst(if,jf+n+1) +
     &                2.d0 * (fdst(if+idir,jf+n) - fdst(if,jf+n))) +
     &                               sigmaf(if+ii-1,jf+n,2) *
     &               (fdst(if+idir,jf+n+1) - fdst(if+idir,jf+n) +
     &                2.d0 * (fdst(if,jf+n+1) - fdst(if,jf+n))) -
     &                               sigmaf(if+ii-1,jf+n,3) * idir *
     &               3.d0 * (fdst(if+idir,jf+n+1) - fdst(if,jf+n)))
               end do
            end if
         end do
      end do
c weighting
      res(if,jf) = src(if,jf) - sum / 6.d0
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
      do j = regl1, regh1
         do i = regl0, regh0
            cen(i,j) = 6.d0 / (2.d0 * (sig(i-1,j-1,1) + sig(i-1,j,1) +
     &                                 sig(i,j-1,1) + sig(i,j,1) +
     &                                 sig(i-1,j-1,2) + sig(i-1,j,2) +
     &                                 sig(i,j-1,2) + sig(i,j,2)) +
     &                         3.d0 * (sig(i-1,j,3) - sig(i-1,j-1,3) +
     &                                 sig(i,j-1,3) - sig(i,j,3)))
         end do
      end do
c$$$      write(unit = 10, fmt = *) 'cen - hgcen_full'
c$$$      write(unit = 10, fmt = *) cen
c$$$      write(unit = 10, fmt = *) 'sig'
c$$$      write(unit = 10, fmt = *) sig
      end
c-----------------------------------------------------------------------
      subroutine hgrlx_terrain(
     & cor,   corl0, corh0, corl1, corh1,
     & res,   resl0, resh0, resl1, resh1,
     & sig,   sfl0, sfh0, sfl1, sfh1,
     & cen,   cenl0, cenh0, cenl1, cenh1,
     &        regl0, regh0, regl1, regh1)
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
      do j = regl1, regh1
cdir$ ivdep
         do i = regl0, regh0
               cor(i,j) = (AVG() - res(i,j)) * cen(i,j)
         end do
      end do
c$$$      write(unit = 10, fmt = *) 'cor'
c$$$      write(unit = 10, fmt = *) cor
c$$$      write(unit = 10, fmt = *) 'sig'
c$$$      write(unit = 10, fmt = *) sig
c$$$      write(unit = 10, fmt = *) 'res'
c$$$      write(unit = 10, fmt = *) res
c$$$      write(unit = 10, fmt = *) 'cen'
c$$$      write(unit = 10, fmt = *) cen
c$$$      stop
      end
c-----------------------------------------------------------------------
      subroutine hgres_terrain(
     & res,   resl0, resh0, resl1, resh1,
     & src,   srcl0, srch0, srcl1, srch1,
     & dest,  destl0, desth0, destl1, desth1,
     & sig,   sfl0, sfh0, sfl1, sfh1,
     & cen,   cenl0, cenh0, cenl1, cenh1,
     &        regl0, regh0, regl1, regh1)
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
      do j = regl1, regh1
         do i = regl0, regh0
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
               res(i,j) = tmp + 0.5d0 *
     &           (sig(i-1,j-1,3) *
     &             (dest(i-1,j-1) - dest(i,j)) -
     &            sig(i-1,j,3) *
     &             (dest(i-1,j+1) - dest(i,j)) -
     &            sig(i,j-1,3) *
     &             (dest(i+1,j-1) - dest(i,j)) +
     &            sig(i,j,3) *
     &             (dest(i+1,j+1) - dest(i,j)))
         end do
      end do
      end
c-----------------------------------------------------------------------

      subroutine hgrlnf_terrain(
     & cor,   corl0, corh0, corl1, corh1,
     & res,   resl0, resh0, resl1, resh1,
     & wrk,   wrkl0, wrkh0, wrkl1, wrkh1,
     & sig,   sfl0, sfh0, sfl1, sfh1,
     & cen,   cenl0, cenh0, cenl1, cenh1,
     &        regl0, regh0, regl1, regh1,
     &        doml0, domh0, doml1, domh1,
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
         do j = regl1 + ioff, regh1, 2
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
         end do
c forward solve loop:
         do i = regl0 + 1, regh0 - 1
            do j = regl1 + ioff, regh1, 2
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
               wrk(i,j) = aj * betm
            end do
         end do
         i = regh0
         do j = regl1 + ioff, regh1, 2
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
         end do
      else
         if (mod(regl0,2) .eq. 0) then
            ioff = ipass
         else
            ioff = 1 - ipass
         end if
         j = regl1
         do i = regl0 + ioff, regh0, 2
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
         end do
c forward solve loop:
         do j = regl1 + 1, regh1 - 1
            do i = regl0 + ioff, regh0, 2
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
               wrk(i,j) = aj * betm 
            end do
         end do
         j = regh1
         do i = regl0 + ioff, regh0, 2
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
         end do
      end if
      end

      
