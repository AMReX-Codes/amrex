c-----------------------------------------------------------------------
      subroutine hgavg(
     &     src,   srcl0, srch0, srcl1, srch1,
     &     rf,    fl0, fh0, fl1, fh1,
     &     fregl0, fregh0, fregl1, fregh1,
     &     hx, irz, imax, idense)
      integer srcl0, srch0, srcl1, srch1
      integer fl0, fh0, fl1, fh1
      integer fregl0, fregh0, fregl1, fregh1
      double precision src(srcl0:srch0,srcl1:srch1)
      double precision rf(fl0:fh0,fl1:fh1)
      double precision hx
      integer irz, imax
      double precision fac, r1, r0m, r1m
      integer i, j
      integer idense
      fac = 0.25d0
      do j = fregl1, fregh1
         do i = fregl0, fregh0
c     We *only* want this modification for the cross stencil
            if (irz .eq. 1 .and. i .eq. 0) then
               src(i,j) = src(i,j) + fac * (rf(i,j-1) + rf(i,j))
            else
               src(i,j) = src(i,j) + fac * (rf(i-1,j-1) + rf(i-1,j) +
     &              rf(i,j-1)   + rf(i,j))
            endif
         end do
      end do

      if ( idense .eq. 1 ) then
c     We dont want these extra terms for the cross stencil
         if (irz .eq. 1) then
            fac = hx / 24.d0
            r1  = (fregl0 - 0.5d0) * hx
            r1m = 1.d0 / r1
            do i = fregl0, fregh0
               r0m = r1m
               r1  = r1 + hx
               if (i .lt. imax) then
                  r1m = 1.d0 / r1
               else
                  r1m = -r0m
               end if
               do j = fregl1, fregh1
                  src(i,j) = src(i,j) + fac *
     &                 (r0m * (rf(i-1,j) + rf(i-1,j-1)) -
     &                 r1m * (rf(i,j-1) + rf(i,j)))
               end do
            end do
         end if
      end if
      end
c-----------------------------------------------------------------------
c Note---only generates values at coarse points along edge of fine grid
      subroutine hgfavg(
     & src,   srcl0, srch0, srcl1, srch1,
     & rc,    cl0, ch0, cl1, ch1,
     & rf,    fl0, fh0, fl1, fh1,
     &        cregl0, cregh0, cregl1, cregh1,
     & ir, jr, idim, idir,
     & hx, irz, imax, idense)
      integer srcl0, srch0, srcl1, srch1
      integer cl0, ch0, cl1, ch1
      integer fl0, fh0, fl1, fh1
      integer cregl0, cregh0, cregl1, cregh1
      double precision src(srcl0:srch0,srcl1:srch1)
      double precision rc(cl0:ch0,cl1:ch1)
      double precision rf(fl0:fh0,fl1:fh1)
      double precision hx
      integer ir, jr, idim, idir, irz, imax
      double precision fac0, fac, r, rfac, rfac0, rfac1
      double precision rfac0m, rfac1m, rfac0p, rfac1p
      integer i, j, irc, jrc, irf, jrf, m, n
      integer idense

      if (idim .eq. 0) then
         i = cregl0
         if (idir .eq. 1) then
            irc = i - 1
            irf = i * ir
         else
            irc = i
            irf = i * ir - 1
         end if
         fac0 = 0.5d0 * ir / (ir+1)
         do j = cregl1, cregh1
          src(i*ir,j*jr) = src(i*ir,j*jr) + fac0 *
     &         (rc(irc,j) + rc(irc,j-1))
         end do
         if ( idense .eq. 1 ) then
c     We dont want these extra terms for the cross stencil
            if (irz .eq. 1) then
               r = (irc + 0.5d0) * (hx * ir)
            rfac = (hx * ir) / (6.d0 * r)
            do j = cregl1, cregh1
               src(i*ir,j*jr) = src(i*ir,j*jr) + fac0 *
     &           (rfac * idir * (rc(irc,j) + rc(irc,j-1)))
            end do
            r = (irf + 0.5d0) * hx
            rfac = hx / (6.d0 * r)
         end if
         end if
         fac0 = fac0 / (ir * jr * jr)
         i = i * ir
         do n = 0, jr-1
            fac = (jr-n) * fac0
            if (n .eq. 0) fac = 0.5d0 * fac
            do j = jr*cregl1, jr*cregh1, jr
               src(i,j) = src(i,j) + fac *
     &           (rf(irf,j-n) + rf(irf,j-n-1) +
     &            rf(irf,j+n) + rf(irf,j+n-1))
            end do
            if ( idense .eq. 1 ) then
c     We dont want these extra terms for the cross stencil
            if (irz .eq. 1) then
               do j = jr*cregl1, jr*cregh1, jr
                  src(i,j) = src(i,j) - fac *
     &              (rfac * idir * (rf(irf,j-n) + rf(irf,j-n-1) +
     &                              rf(irf,j+n) + rf(irf,j+n-1)))
               end do
            end if
            end if
         end do
      else if (idim .eq. 1) then
         j = cregl1
         if (idir .eq. 1) then
            jrc = j - 1
            jrf = j * jr
         else
            jrc = j
            jrf = j * jr - 1
         end if
         fac0 = 0.5d0 * jr / (jr+1)
         do i = cregl0, cregh0
            src(i*ir,j*jr) = src(i*ir,j*jr) + fac0 *
     &        (rc(i,jrc) + rc(i-1,jrc))
         end do

         if (irz .eq. 1 .and. cregl0 .le. 0 .and. cregh0 .ge. 0) then
            i = 0
            src(i*ir,j*jr) = src(i*ir,j*jr) - fac0 * rc(i-1,jrc)
         endif

         if ( idense .eq. 1 ) then
c     We dont want these extra terms for the cross stencil
         if (irz .eq. 1 .and. cregh0 .lt. imax) then
            do i = cregl0, cregh0
               r = (i + 0.5d0) * (hx * ir)
               rfac0 = (hx * ir) / (6.d0 * (r - hx * ir))
               rfac1 = (hx * ir) / (6.d0 * r)
               src(i*ir,j*jr) = src(i*ir,j*jr) + fac0 *
     &           (rfac0 * rc(i-1,jrc) - rfac1 * rc(i,jrc))
            end do
         else if (irz .eq. 1) then
c This should only occur with a corner at the outer boundary, which
c should be handled by the cavg routine instead:
            i = cregh0
            r = (i - 0.5d0) * (hx * ir)
            rfac0 = (hx * ir) / (6.d0 * r)
            rfac1 = -rfac0
            src(i*ir,j*jr) = src(i*ir,j*jr) + fac0 *
     &           (rfac0 * rc(i-1,jrc) - rfac1 * rc(i,jrc))
         end if
         end if
         fac0 = fac0 / (ir * ir * jr)
         j = j * jr
         do m = 0, ir-1
            fac = (ir-m) * fac0
            if (m .eq. 0) fac = 0.5d0 * fac
            do i = ir*cregl0, ir*cregh0, ir
               src(i,j) = src(i,j) + fac *
     &           (rf(i-m,jrf) + rf(i-m-1,jrf) +
     &            rf(i+m,jrf) + rf(i+m-1,jrf))
            end do

            if (irz .eq. 1 .and. m .eq. 0 .and.
     @          cregl0 .le. 0 .and. cregh0 .ge. 0) then
               i = 0
               src(i*ir,j*jr) = src(i*ir,j*jr) - fac *
     @              (rf(i-m-1,jrf) + rf(i+m-1,jrf))
            endif
            if ( idense .eq. 1 ) then
c     We dont want these extra terms for the cross stencil
            if (irz .eq. 1 .and. cregh0 .lt. imax) then
               do i = ir*cregl0, ir*cregh0, ir
                  r = (i + 0.5d0) * hx
                  rfac0m = hx / (6.d0 * (r - (m + 1) * hx))
                  rfac1m = hx / (6.d0 * (r - m * hx))
                  rfac0p = hx / (6.d0 * (r + (m - 1) * hx))
                  rfac1p = hx / (6.d0 * (r + m * hx))
                  src(i,j) = src(i,j) + fac *
     &              (rfac0m * rf(i-m-1,jrf) - rfac1m * rf(i-m,jrf) +
     &               rfac0p * rf(i+m-1,jrf) - rfac1p * rf(i+m,jrf))
               end do
            else if (irz .eq. 1) then
c This should only occur with a corner at the outer boundary, which
c should be handled by the cavg routine instead:
               i = ir*cregh0
               r = (i + 0.5d0) * hx
               rfac0m = hx / (6.d0 * (r - (m + 1) * hx))
               if (m .eq. 0) then
                  rfac1m = -rfac0m
               else
                  rfac1m = hx / (6.d0 * (r - m * hx))
               end if
               rfac0p = -rfac1m
               rfac1p = -rfac0m
               src(i,j) = src(i,j) + fac *
     &              (rfac0m * rf(i-m-1,jrf) - rfac1m * rf(i-m,jrf) +
     &               rfac0p * rf(i+m-1,jrf) - rfac1p * rf(i+m,jrf))
            end if
            end if
         end do
      end if
      end
c-----------------------------------------------------------------------
c Note---only generates values at coarse points along edge of fine grid
      subroutine hgcavg(
     & src,   srcl0, srch0, srcl1, srch1,
     & rc,    cl0, ch0, cl1, ch1,
     & rf,    fl0, fh0, fl1, fh1,
     &        cregl0, cregh0, cregl1, cregh1,
     & ir, jr, ga, idd,
     & hx, irz, imax, idense)
      integer srcl0, srch0, srcl1, srch1
      integer cl0, ch0, cl1, ch1
      integer fl0, fh0, fl1, fh1
      integer cregl0, cregh0, cregl1, cregh1
      double precision src(srcl0:srch0,srcl1:srch1)
      double precision rc(cl0:ch0,cl1:ch1)
      double precision rf(fl0:fh0,fl1:fh1)
      double precision hx
      integer ir, jr, ga(0:1,0:1), irz, imax, idd
      double precision rm2, sum, center, fac, ffac, cfac, fac1
      double precision r, rfac, rfac0, rfac1
      integer ic, jc, if, jf, ii, ji, irc, jrc, irf, jrf
      integer m, n, idir, jdir
      integer idense
      rm2 = 1.d0 / (ir * jr)
      ic = cregl0
      jc = cregl1
      if = ic * ir
      jf = jc * jr

      sum = 0.d0
      center = 0.d0
c quadrants
      fac = rm2
      ffac = rm2
      cfac = 1.d0
      do ji = 0, 1
         jrf = jf + ji - 1
         jrc = jc + ji - 1
         do ii = 0, 1
            if (ga(ii,ji) .eq. 1) then
               irf = if + ii - 1
               center = center + ffac
               sum = sum + fac * rf(irf,jrf)
               if ( idense .eq. 0 ) then
c     We dont want these extra terms for the cross stencil
               if (irz .eq. 1) then
                  idir = 2 * ii - 1
                  r = (irf + 0.5d0) * hx
                  if (irf .lt. (ir * imax)) then
                     rfac =  hx / (6.d0 * r)
                  else
                     rfac = -hx / (6.d0 * (r - hx))
                  end if
                  sum = sum - fac * rfac * idir * rf(irf,jrf)
               end if
               end if
            else
               irc = ic + ii - 1
               center = center + cfac
               sum = sum + rc(irc,jrc)
               if ( idense .eq. 1 ) then
c     We dont want these extra terms for the cross stencil
               if (irz .eq. 1) then
                  idir = 2 * ii - 1
                  r = (irc + 0.5d0) * (hx * ir)
                  if (irc .lt. imax) then
                     rfac =  (hx * ir) / (6.d0 * r)
                  else
                     rfac = -(hx * ir) / (6.d0 * (r - (hx * ir)))
                  end if
                  sum = sum - rfac * idir * rc(irc,jrc)
               end if
               end if
            end if
         end do
      end do

c     We *only* want this modification for the 5-point stencil:
c      here we halve the average at the r=0 axis.
      if (irz .eq. 1 .and. ic .eq. 0) then
         sum = sum * 0.5d0
      endif

c edges
      do ji = 0, 1
         jdir = 2 * ji - 1
         jrf = jf + ji - 1
         do ii = 0, 1
            idir = 2 * ii - 1
            irf = if + ii - 1
            if (ga(ii,ji) - ga(ii,1-ji) .eq. 1) then
               fac1 = rm2 / ir
               ffac = (ir-1) * rm2
               center = center + ffac
               do m = idir, idir*(ir-1), idir
                  fac = (ir-abs(m)) * fac1
                  sum = sum + fac * (rf(if+m,jrf) + rf(if+m-1,jrf))
                  if ( idense .eq. 1 ) then
c     We dont want these extra terms for the cross stencil
                  if (irz .eq. 1) then
                     r = (if + 0.5d0) * hx
                     if (irf .lt. (ir * imax)) then
                        rfac0 =  hx / (6.d0 * (r + (m - 1) * hx))
                        rfac1 =  hx / (6.d0 * (r + m * hx))
                     else
                        rfac0 = -hx / (6.d0 * (r - m * hx))
                        rfac1 = -hx / (6.d0 * (r - (m + 1) * hx))
                     end if
                     sum = sum + fac *
     &                 (rfac0 * rf(if+m-1,jrf) - rfac1 * rf(if+m,jrf))
                  end if
                  end if
               end do
            end if
            if (ga(ii,ji) - ga(1-ii,ji) .eq. 1) then
               fac1 = rm2 / jr
               ffac = (jr-1) * rm2
               center = center + ffac
               do n = jdir, jdir*(jr-1), jdir
                  fac = (jr-abs(n)) * fac1
                  sum = sum + fac * (rf(irf,jf+n) + rf(irf,jf+n-1))
                  if ( idense .eq. 1 ) then
c     We dont want these extra terms for the cross stencil
                  if (irz .eq. 1) then
                     r = (irf + 0.5d0) * hx
                     rfac = hx / (6.d0 * r)
                     sum = sum - fac *
     &                 (rfac * idir * (rf(irf,jf+n) + rf(irf,jf+n-1)))
                  end if
                  end if
               end do
            end if
         end do
      end do
c weighting
      src(if,jf) = src(if,jf) + sum / center
      end
