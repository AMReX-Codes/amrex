c variable density versions:

c Note---assumes fdst linearly interpolated from cdst along edge
      subroutine hgfres(
     & res,  resl0, resh0, resl1, resh1,
     & src,  srcl0, srch0, srcl1, srch1,
     & fdst, fdstl0, fdsth0, fdstl1, fdsth1,
     & cdst, cdstl0, cdsth0, cdstl1, cdsth1,
     & sigmaf, sfl0, sfh0, sfl1, sfh1,
     & sigmac, scl0, sch0, scl1, sch1,
     &       regl0, regh0, regl1, regh1,
     & hx, hy, ir, jr, idim, idir, irz)
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
      integer irz
      integer ir, jr, idim, idir
      double precision hxm2, hym2, fac0, fac1, tmp
      integer i, is, j, js, m, n

      if (irz .eq. 1 .and. regl0 .le. 0 .and. regh0 .ge. 0) then
         print *,'I DONT THINK WE SHOULD BE IN HGFRES AT I=0 '
         stop
      endif

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
         do j = regl1, regh1
            res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     &        (hxm2 *
     &          ((sigmac(is,j-1) + sigmac(is,j)) *
     &            (cdst(i-idir,j) - cdst(i,j))) +
     &         hym2 *
     &          (sigmac(is,j-1) *
     &            (cdst(i,j-1) - cdst(i,j)) +
     &           sigmac(is,j) *
     &            (cdst(i,j+1) - cdst(i,j))))
         end do
         fac0 = fac0 / (ir * jr * jr)
         hxm2 = ir * ir * hxm2
         hym2 = jr * jr * hym2
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
               tmp = hxm2 *
     &           ((sigmaf(is,j-n-1) + sigmaf(is,j-n)) *
     &             (fdst(i+idir,j-n) - fdst(i,j-n)) +
     &           (sigmaf(is,j+n-1) + sigmaf(is,j+n)) *
     &             (fdst(i+idir,j+n) - fdst(i,j+n)))
               res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     &            (sigmaf(is,j-n-1) *
     &             (fdst(i,j-n-1) - fdst(i,j-n)) +
     &             sigmaf(is,j-n) *
     &             (fdst(i,j-n+1) - fdst(i,j-n)) +
     &             sigmaf(is,j+n-1) *
     &             (fdst(i,j+n-1) - fdst(i,j+n)) +
     &             sigmaf(is,j+n) *
     &             (fdst(i,j+n+1) - fdst(i,j+n))))
            end do
         end do
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
         do i = regl0, regh0
            res(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     &        (hxm2 *
     &          (sigmac(i-1,js) *
     &            (cdst(i-1,j) - cdst(i,j)) +
     &           sigmac(i,js) *
     &            (cdst(i+1,j) - cdst(i,j))) +
     &         hym2 *
     &          ((sigmac(i-1,js) + sigmac(i,js)) *
     &            (cdst(i,j-idir) - cdst(i,j))))
         end do

c        This correction is *only* for the cross stencil
         if (irz .eq. 1 .and. regl0 .le. 0 .and. regh0 .ge. 0) then
            i = 0
            res(i*ir,j*jr) = res(i*ir,j*jr) + fac0 *
     &         hym2 * 0.5d0 *
     &          ((sigmac(i-1,js) + sigmac(i,js)) *
     &            (cdst(i,j-idir) - cdst(i,j)))
         endif

         fac0 = fac0 / (ir * ir * jr)
         hxm2 = ir * ir * hxm2
         hym2 = jr * jr * hym2
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
               tmp = hxm2 *
     &             (sigmaf(i-m-1,js) *
     &             (fdst(i-m-1,j) - fdst(i-m,j)) +
     &              sigmaf(i-m,js) *
     &             (fdst(i-m+1,j) - fdst(i-m,j)) +
     &              sigmaf(i+m-1,js) *
     &             (fdst(i+m-1,j) - fdst(i+m,j)) +
     &              sigmaf(i+m,js) *
     &             (fdst(i+m+1,j) - fdst(i+m,j)))
               res(i,j) = res(i,j) - fac1 * (tmp + hym2 *
     &           ((sigmaf(i-m-1,js) + sigmaf(i-m,js)) *
     &             (fdst(i-m,j+idir) - fdst(i-m,j)) +
     &              (sigmaf(i+m-1,js) + sigmaf(i+m,js)) *
     &             (fdst(i+m,j+idir) - fdst(i+m,j))))
            end do

            if (irz .eq. 1 .and. m .eq. 0 .and.
     &          regl0 .le. 0 .and. regh0 .ge. 0) then
               i = 0
               res(i,j) = res(i,j) + fac1 * hym2 * 0.5d0 *
     &           ((sigmaf(i-m-1,js) + sigmaf(i-m,js)) *
     &             (fdst(i-m,j+idir) - fdst(i-m,j)) +
     &              (sigmaf(i+m-1,js) + sigmaf(i+m,js)) *
     &             (fdst(i+m,j+idir) - fdst(i+m,j)))
            endif
         end do
      end if
      end

c Note---assumes fdst linearly interpolated from cdst along edges
      subroutine hgcres(
     & res,   resl0, resh0, resl1, resh1,
     & src,   srcl0, srch0, srcl1, srch1,
     & fdst,  fdstl0, fdsth0, fdstl1, fdsth1,
     & cdst,  cdstl0, cdsth0, cdstl1, cdsth1,
     & sigmaf, sfl0, sfh0, sfl1, sfh1,
     & sigmac, scl0, sch0, scl1, sch1,
     &        regl0, regh0, regl1, regh1,
     & hx, hy, ir, jr, ga, irz)
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
      integer ir, jr, ga(0:1,0:1), irz
      double precision hxm2, hym2, hxm2c, hym2c, sum, center,
     &   ffac, cfac, fac, fac1
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
      do ji = 0, 1
         jdir = 2 * ji - 1
         do ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,ji) .eq. 1) then
               center = center + ffac
               sum = sum + sigmaf(if+ii-1,jf+ji-1) *
     &           (hxm2 * (fdst(if+idir,jf) - fdst(if,jf)) +
     &            hym2 * (fdst(if,jf+jdir) - fdst(if,jf)))
               if (irz .eq. 1 .and. ic .eq. 0) then
                 sum = sum - sigmaf(if+ii-1,jf+ji-1) * 0.5d0 *
     &             (hym2 * (fdst(if,jf+jdir) - fdst(if,jf)))
               endif
            else
               center = center + cfac
               sum = sum + ir * jr * sigmac(ic+ii-1,jc+ji-1) *
     &           (hxm2c * (cdst(ic+idir,jc) - cdst(ic,jc)) +
     &            hym2c * (cdst(ic,jc+jdir) - cdst(ic,jc)))
               if (irz .eq. 1 .and. ic .eq. 0) then
                 sum = sum - ir * jr * sigmac(ic+ii-1,jc+ji-1) * 0.5d0 *
     &             (hym2c * (cdst(ic,jc+jdir) - cdst(ic,jc)))
               endif
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
               ffac = 0.5d0 * (ir-1)
               center = center + ffac
               do m = idir, idir*(ir-1), idir
                  fac = (ir-abs(m)) * fac1
                  sum = sum + fac *
     &              (hxm2 * (sigmaf(if+m-1,jf+ji-1) *
     &                        (fdst(if+m-1,jf) - fdst(if+m,jf)) +
     &                       sigmaf(if+m,jf+ji-1) *
     &                        (fdst(if+m+1,jf) - fdst(if+m,jf))) +
     &               hym2 *
     &                 (sigmaf(if+m-1,jf+ji-1) + sigmaf(if+m,jf+ji-1)) *
     &                 (fdst(if+m,jf+jdir) - fdst(if+m,jf)))
               end do
            end if
            if (ga(ii,ji) - ga(1-ii,ji) .eq. 1) then
               fac1 = 1.d0 / jr
               ffac = 0.5d0 * (jr-1)
               center = center + ffac
               do n = jdir, jdir*(jr-1), jdir
                  fac = (jr-abs(n)) * fac1
                  sum = sum + fac *
     &              (hxm2 *
     &                 (sigmaf(if+ii-1,jf+n-1) + sigmaf(if+ii-1,jf+n)) *
     &                 (fdst(if+idir,jf+n) - fdst(if,jf+n)) +
     &               hym2 * (sigmaf(if+ii-1,jf+n-1) *
     &                        (fdst(if,jf+n-1) - fdst(if,jf+n)) +
     &                       sigmaf(if+ii-1,jf+n) *
     &                        (fdst(if,jf+n+1) - fdst(if,jf+n))))
               end do
            end if
         end do
      end do
c weighting
      res(if,jf) = src(if,jf) - sum / center
      end
c-----------------------------------------------------------------------

      subroutine hgcen(
     & cen,   cenl0,cenh0,cenl1,cenh1,
     & signd, snl0,snh0,snl1,snh1,
     &        regl0,regh0,regl1,regh1,irz)
      integer cenl0,cenh0,cenl1,cenh1
      integer snl0,snh0,snl1,snh1
      integer regl0,regh0,regl1,regh1
      double precision cen(cenl0:cenh0,cenl1:cenh1)
      double precision signd(snl0:snh0,snl1:snh1, 2)
      double precision tmp
      integer irz
      integer i, j
      do j = regl1, regh1
         do i = regl0, regh0
            tmp = (signd(i-1,j,1) + signd(i,j,1) 
     &           + signd(i,j-1,2) + signd(i,j,2))
            if ( tmp .eq. 0.0D0 ) then
               cen(i,j) = 0.0D0
            else
               cen(i,j) = 1.0D0 / tmp
            end if
         end do
         if (irz .eq. 1 .and. regl0 .eq. 0) then
            i = 0
            tmp = (signd(i-1,j,1) + signd(i,j,1)
     &           + 0.5d0*(signd(i,j-1,2) + signd(i,j,2)))
            if ( tmp .eq. 0.0D0 ) then
               cen(i,j) = 0.0D0
            else
               cen(i,j) = 1.0D0 / tmp 
            end if
         end if
      end do
      end

c-----------------------------------------------------------------------

      subroutine hgres_cross(
     & res, resl0,resh0,resl1,resh1,
     & src, dest, signd,
     &      regl0,regh0,regl1,regh1,
     & irz)
      integer resl0,resh0,resl1,resh1
      integer regl0,regh0,regl1,regh1
      double precision res(resl0:resh0, resl1:resh1)
      double precision src(resl0:resh0, resl1:resh1)
      double precision dest(resl0:resh0, resl1:resh1)
      double precision signd(resl0:resh0, resl1:resh1,2)
      integer irz
      integer istart, iend
      integer i, j, jdiff, ly
      jdiff = resh0 - resl0 + 1
      ly = (resh1 - resl1 + 1) * jdiff
      istart = (regl1 - resl1) * jdiff + (regl0 - resl0)
      iend   = (regh1 - resl1) * jdiff + (regh0 - resl0)

      do j = regl1, regh1
          do i = regl0, regh0
          res(i,j) = (src(i,j) - (
     &    + signd(i-1,j,1)*(dest(i-1,j)-dest(i,j))
     &    + signd(i,j,1)  *(dest(i+1,j)-dest(i,j))
     &    + signd(i,j-1,2)*(dest(i,j-1)-dest(i,j))
     &    + signd(i,j,2)  *(dest(i,j+1)-dest(i,j))
     &    )
     &    )
        end do
        end do

      if (irz .eq. 1 .and. regl0 .eq. 0) then
        do j = regl1, regh1
        do i = regl0, regh0
           res(i,j) = (src(i,j) -
     &       (signd(i-1,j,1) * (dest(i-1,j) - dest(i,j)) +
     &        signd(i,j,1)   * (dest(i+1,j) - dest(i,j)) +
     &        signd(i,j-1,2) * (dest(i,j-1) - dest(i,j)) * 0.5d0 +
     &        signd(i,j,2)   * (dest(i,j+1) - dest(i,j)) * 0.5d0 ))
        end do
        end do
      endif
      end
c-----------------------------------------------------------------------
      subroutine hgscon(
     & signd, snl0,snh0,snl1,snh1,
     & sigx, sigy,
     &        scl0,sch0,scl1,sch1,
     &        regl0,regh0,regl1,regh1,
     & hx, hy)
      integer snl0,snh0,snl1,snh1
      integer scl0,sch0,scl1,sch1
      integer regl0,regh0,regl1,regh1
      double precision signd(snl0:snh0,snl1:snh1, 2)
      double precision sigx(scl0:sch0,scl1:sch1)
      double precision sigy(scl0:sch0,scl1:sch1)
      double precision hx, hy
      double precision facx, facy
      integer i, j
      facx = 0.5D0 / (hx*hx)
      facy = 0.5D0 / (hy*hy)
         do j = regl1, regh1
            do i = regl0-1, regh0
               signd(i,j,1) = facx *
     &               (sigx(i,j) + sigx(i,j-1))
            end do
         end do
         do j = regl1-1, regh1
            do i = regl0, regh0
               signd(i,j,2) = facy *
     &               (sigy(i-1,j) + sigy(i,j))
            end do
         end do
      end

c-----------------------------------------------------------------------

      subroutine hgrlx(
     & cor, res, sig, cen,
     &     resl0,resh0,resl1,resh1,
     &     regl0,regh0,regl1,regh1,irz)
      integer resl0,resh0,resl1,resh1
      integer regl0,regh0,regl1,regh1
      double precision cor(resl0:resh0,resl1:resh1)
      double precision res(resl0:resh0,resl1:resh1)
      double precision sig(resl0:resh0,resl1:resh1,2)
      double precision cen(resl0:resh0,resl1:resh1)
      double precision AVG
      double precision AVGREDGE
      integer irz
      integer istart, iend
      integer i, j, jdiff, ly, ipar
      AVGREDGE() = (sig(i-1,j,1) * cor(i-1,j) +
     &              sig(i,j,1)   * cor(i+1,j) +
     &              sig(i,j-1,2) * cor(i,j-1) * 0.5d0 +
     &              sig(i,j,2)   * cor(i,j+1) * 0.5d0 )
      AVG() = (sig(i-1,j,1)        * cor(i-1,j) +
     &         sig(i,j,1)          * cor(i+1,j) +
     &         sig(i,j-1,2)        * cor(i,j-1) +
     &         sig(i,j,2)          * cor(i,j+1))
      jdiff =  resh0 - resl0 + 1
      ly    = (resh1 - resl1 + 1) * jdiff
      istart = (regl1 - resl1) * jdiff + (regl0 - resl0)
      iend   = (regh1 - resl1) * jdiff + (regh0 - resl0)
      if (irz .eq. 0 .or. regl0 .gt. 0) then
        ipar = 1
        do j = regl1, regh1
          ipar = 1 - ipar
          do i = regl0 + ipar, regh0, 2
            cor(i,j) = (AVG()-res(i,j))*cen(i,j)
          end do
        end do
        ipar = 0
        do j = regl1, regh1
          ipar = 1 - ipar
          do i = regl0+ipar, regh0, 2
            cor(i,j) = (AVG()-res(i,j))*cen(i,j)
          end do
        end do
      else
c     Now irz = 1 and regl0 = 0, so we are touching the r=0 edge
        ipar = 1
        do j = regl1, regh1
          ipar = 1 - ipar
          do i = regl0 + ipar, regh0, 2
            if (i .eq. 0) then
              cor(i,j) = (AVGREDGE() - res(i,j)) * cen(i,j)
            else
              cor(i,j) = (AVG() - res(i,j)) * cen(i,j)
            endif
          end do
        end do
        ipar = 0
        do j = regl1, regh1
          ipar = 1 - ipar
          do i = regl0 + ipar, regh0, 2
            if ( i .eq. 0) then
              cor(i,j) = (AVGREDGE() - res(i,j)) * cen(i,j)
            else
              cor(i,j) = (AVG() - res(i,j)) * cen(i,j)
            endif
          end do
        end do
      endif
      end

c-----------------------------------------------------------------------
c NODE-based data, factor of 2 only.
      subroutine hgints(
     & dest,  destl0,desth0,destl1,desth1,
     &        regl0,regh0,regl1,regh1,
     & signd, snl0,snh0,snl1,snh1,
     & src,   srcl0,srch0,srcl1,srch1,
     &        bbl0,bbh0,bbl1,bbh1,
     & ir, jr)
      integer destl0,desth0,destl1,desth1
      integer regl0,regh0,regl1,regh1
      integer snl0,snh0,snl1,snh1
      integer srcl0,srch0,srcl1,srch1
      integer bbl0,bbh0,bbl1,bbh1
      integer ir, jr
      double precision dest(destl0:desth0,destl1:desth1)
      double precision signd(snl0:snh0,snl1:snh1, 2)
      double precision src(srcl0:srch0,srcl1:srch1)
      integer i, j, ic, jc
      do jc = bbl1, bbh1
         do ic = bbl0, bbh0
            dest(ir*ic,jr*jc) = src(ic,jc)
         end do
      end do
      if (ir .eq. 2) then
         do jc = bbl1, bbh1
            do ic = bbl0, bbh0-1
               i = ir * ic
               j = jr * jc
               dest(i+1,j) = (signd(i,j,1)  * src(ic,jc) +
     &                        signd(i+1,j,1) * src(ic+1,jc)) /
     &                       (signd(i,j,1) + signd(i+1,j,1))
            end do
         end do
      end if
      if (jr .eq. 2) then
         do jc = bbl1, bbh1-1
            do ic = bbl0, bbh0
               i = ir * ic
               j = jr * jc
               dest(i,j+1) = (signd(i,j,2)  * src(ic,jc) +
     &                        signd(i,j+1,2) * src(ic,jc+1)) /
     &                       (signd(i,j,2) + signd(i,j+1,2))
            end do
         end do
      end if
      if (ir .eq. 2 .and. jr .eq. 2) then
         do jc = bbl1, bbh1-1
            do ic = bbl0, bbh0-1
               i = ir * ic
               j = jr * jc
               dest(i+1,j+1) = (signd(i,j+1,1)   * dest(i,j+1) +
     &                          signd(i+1,j+1,1) * dest(i+2,j+1) +
     &                          signd(i+1,j,2)   * dest(i+1,j) +
     &                          signd(i+1,j+1,2) * dest(i+1,j+2)) /
     &                         (signd(i,j+1,1) + signd(i+1,j+1,1) +
     &                          signd(i+1,j,2) + signd(i+1,j+1,2))
            end do
         end do
      end if
      end

c-----------------------------------------------------------------------
c CELL-based data only.
      subroutine hgsrst(
     & destx, desty,
     &     destl0, desth0, destl1, desth1,
     &     regl0, regh0, regl1, regh1,
     & srcx, srcy,
     &     srcl0, srch0, srcl1, srch1,
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
         do j = regl1, regh1
            do i = regl0, regh0
               i2 = 2 * i
               j2 = 2 * j
               destx(i,j) = 1.d0 /
     &                      (1.d0 / (srcx(i2,j2)   + srcx(i2,j2+1)) +
     &                       1.d0 / (srcx(i2+1,j2) + srcx(i2+1,j2+1)))
               desty(i,j) = 1.d0 /
     &                      (1.d0 / (srcy(i2,j2)   + srcy(i2+1,j2)) +
     &                       1.d0 / (srcy(i2,j2+1) + srcy(i2+1,j2+1)))
            end do
         end do
      else if (ir .eq. 2) then
         do j = regl1, regh1
            do i = regl0, regh0
               i2 = 2 * i
               destx(i,j) = 2.d0 /
     &                      (1.d0 / srcx(i2,j) + 1.d0 / srcx(i2+1,j))
               desty(i,j) = 0.5d0 * (srcy(i2,j) + srcy(i2+1,j))
            end do
         end do
      else
         do j = regl1, regh1
            do i = regl0, regh0
               j2 = 2 * j
               destx(i,j) = 0.5d0 * (srcx(i,j2) + srcx(i,j2+1))
               desty(i,j) = 2.d0 /
     &                      (1.d0 / srcy(i,j2) + 1.d0 / srcy(i,j2+1))
            end do
         end do
      end if
      end
