c-----------------------------------------------------------------------
      subroutine hgavg(src, srcl0, srch0, srcl1, srch1, srcl2, srch2,
     & rf, fl0, fh0, fl1, fh1, fl2, fh2,
     & fregl0, fregh0, fregl1, fregh1, fregl2, fregh2)
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer fl0,fh0,fl1,fh1,fl2,fh2
      integer fregl0,fregh0,fregl1,fregh1,fregl2,fregh2
      real*8 src(srcl0:srch0,srcl1:srch1,srcl2:srch2)
      real*8 rf(fl0:fh0,fl1:fh1,fl2:fh2)
      real*8 fac
	parameter (fac = 0.125d0)
      integer i, j, k
      do k = fregl2, fregh2
         do j = fregl1, fregh1
            do i = fregl0, fregh0
               src(i,j,k) = src(i,j,k) + fac *
     &              (rf(i-1,j-1,k-1) + rf(i-1,j-1,k) +
     &               rf(i-1,j,k-1)   + rf(i-1,j,k) +
     &               rf(i,j-1,k-1)   + rf(i,j-1,k) +
     &               rf(i,j,k-1)     + rf(i,j,k))
            end do
         end do
      end do
      end
c-----------------------------------------------------------------------
c Note---only generates values at coarse points along edge of fine grid
      subroutine hgfavg(src,
     & srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & rc,
     & cl0,ch0,cl1,ch1,cl2,ch2,
     & rf,
     & fl0,fh0,fl1,fh1,fl2,fh2,
     & cregl0,cregh0,cregl1,cregh1,cregl2,cregh2,
     & ir, jr, kr, idim, idir)
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer cl0,ch0,cl1,ch1,cl2,ch2
      integer fl0,fh0,fl1,fh1,fl2,fh2
      integer cregl0,cregh0,cregl1,cregh1,cregl2,cregh2
      real*8 src(srcl0:srch0,srcl1:srch1,srcl2:srch2)
      real*8 rc(cl0:ch0,cl1:ch1,cl2:ch2)
      real*8 rf(fl0:fh0,fl1:fh1,fl2:fh2)
      integer ir, jr, kr, idim, idir
      real*8 fac0, fac1, fac
      integer i, j, k, irc, irf, jrc, jrf, krc, krf, l, m, n
      if (idim .eq. 0) then
         i = cregl0
         if (idir .eq. 1) then
            irc = i - 1
            irf = i * ir
         else
            irc = i
            irf = i * ir - 1
         end if
         fac0 = 0.25d0 * ir / (ir+1)
         do k = cregl2, cregh2
            do j = cregl1, cregh1
               src(i*ir,j*jr,k*kr) =
     &            src(i*ir,j*jr,k*kr) + fac0 *
     &              (rc(irc,j,k)   + rc(irc,j,k-1) +
     &               rc(irc,j-1,k) + rc(irc,j-1,k-1))
            end do
         end do
         fac0 = fac0 / (ir * jr * kr * jr * kr)
         i = i * ir
         do l = 0, kr-1
            fac1 = (kr-l) * fac0
            if (l .eq. 0) fac1 = 0.5d0 * fac1
            do n = 0, jr-1
               fac = (jr-n) * fac1
               if (n .eq. 0) fac = 0.5d0 * fac
               do k = kr*cregl2, kr*cregh2, kr
                  do j = jr*cregl1, jr*cregh1, jr
                     src(i,j,k) = src(i,j,k) + fac *
     &                 (rf(irf,j-n,k-l) + rf(irf,j-n,k-l-1) +
     &                  rf(irf,j-n-1,k-l) + rf(irf,j-n-1,k-l-1) +
     &                  rf(irf,j-n,k+l) + rf(irf,j-n,k+l-1) +
     &                  rf(irf,j-n-1,k+l) + rf(irf,j-n-1,k+l-1) +
     &                  rf(irf,j+n,k-l) + rf(irf,j+n,k-l-1) +
     &                  rf(irf,j+n-1,k-l) + rf(irf,j+n-1,k-l-1) +
     &                  rf(irf,j+n,k+l) + rf(irf,j+n,k+l-1) +
     &                  rf(irf,j+n-1,k+l) + rf(irf,j+n-1,k+l-1))
                  end do
               end do
            end do
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
         fac0 = 0.25d0 * jr / (jr+1)
         do k = cregl2, cregh2
            do i = cregl0, cregh0
               src(i*ir,j*jr,k*kr) =
     &            src(i*ir,j*jr,k*kr) + fac0 *
     &           (rc(i,jrc,k)   + rc(i-1,jrc,k) +
     &            rc(i,jrc,k-1) + rc(i-1,jrc,k-1))
            end do
         end do
         fac0 = fac0 / (ir * jr * kr * ir * kr)
         j = j * jr
         do l = 0, kr-1
            fac1 = (kr-l) * fac0
            if (l .eq. 0) fac1 = 0.5d0 * fac1
            do m = 0, ir-1
               fac = (ir-m) * fac1
               if (m .eq. 0) fac = 0.5d0 * fac
               do k = kr*cregl2, kr*cregh2, kr
                  do i = ir*cregl0, ir*cregh0, ir
                     src(i,j,k) = src(i,j,k) + fac *
     &                 (rf(i-m,jrf,k-l) + rf(i-m-1,jrf,k-l) +
     &                  rf(i-m,jrf,k-l-1) + rf(i-m-1,jrf,k-l-1) +
     &                  rf(i-m,jrf,k+l) + rf(i-m-1,jrf,k+l) +
     &                  rf(i-m,jrf,k+l-1) + rf(i-m-1,jrf,k+l-1) +
     &                  rf(i+m,jrf,k-l) + rf(i+m-1,jrf,k-l) +
     &                  rf(i+m,jrf,k-l-1) + rf(i+m-1,jrf,k-l-1) +
     &                  rf(i+m,jrf,k+l) + rf(i+m-1,jrf,k+l) +
     &                  rf(i+m,jrf,k+l-1) + rf(i+m-1,jrf,k+l-1))
                  end do
               end do
            end do
         end do
      else
         k = cregl2
         if (idir .eq. 1) then
            krc = k - 1
            krf = k * kr
         else
            krc = k
            krf = k * kr - 1
         end if
         fac0 = 0.25d0 * kr / (kr+1)
         do j = cregl1, cregh1
            do i = cregl0, cregh0
               src(i*ir,j*jr,k*kr) =
     &            src(i*ir,j*jr,k*kr) + fac0 *
     &           (rc(i,j,krc)   + rc(i-1,j,krc) +
     &            rc(i,j-1,krc) + rc(i-1,j-1,krc))
            end do
         end do
         fac0 = fac0 / (ir * jr * kr * ir * jr)
         k = k * kr
         do n = 0, jr-1
            fac1 = (jr-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do m = 0, ir-1
               fac = (ir-m) * fac1
               if (m .eq. 0) fac = 0.5d0 * fac
               do j = jr*cregl1, jr*cregh1, jr
                  do i = ir*cregl0, ir*cregh0, ir
                     src(i,j,k) = src(i,j,k) + fac *
     &                 (rf(i-m,j-n,krf) + rf(i-m-1,j-n,krf) +
     &                  rf(i-m,j-n-1,krf) + rf(i-m-1,j-n-1,krf) +
     &                  rf(i-m,j+n,krf) + rf(i-m-1,j+n,krf) +
     &                  rf(i-m,j+n-1,krf) + rf(i-m-1,j+n-1,krf) +
     &                  rf(i+m,j-n,krf) + rf(i+m-1,j-n,krf) +
     &                  rf(i+m,j-n-1,krf) + rf(i+m-1,j-n-1,krf) +
     &                  rf(i+m,j+n,krf) + rf(i+m-1,j+n,krf) +
     &                  rf(i+m,j+n-1,krf) + rf(i+m-1,j+n-1,krf))
                  end do
               end do
            end do
         end do
      end if
      end
c-----------------------------------------------------------------------
c Note---only generates values at coarse points along edge of fine grid
      subroutine hgeavg(src,
     & srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & rc,
     & cl0,ch0,cl1,ch1,cl2,ch2,
     & rf,
     & fl0,fh0,fl1,fh1,fl2,fh2,
     & cregl0,cregh0,cregl1,cregh1,cregl2,cregh2,
     & ir, jr, kr, ivect, ga)
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer cl0,ch0,cl1,ch1,cl2,ch2
      integer fl0,fh0,fl1,fh1,fl2,fh2
      integer cregl0,cregh0,cregl1,cregh1,cregl2,cregh2
      real*8 src(srcl0:srch0,srcl1:srch1,srcl2:srch2)
      real*8 rc(cl0:ch0,cl1:ch1,cl2:ch2)
      real*8 rf(fl0:fh0,fl1:fh1,fl2:fh2)
      integer ir, jr, kr, ivect(0:2), ga(0:1,0:1,0:1)
      real*8 rm3, center, cfac, ffac, fac0, fac1, fac
      integer ic, jc, kc, if, jf, kf, ii, ji, ki
      integer irc, irf, jrc, jrf, krc, krf, idir, jdir, kdir, l, m, n
      rm3 = 1.d0 / (ir * jr * kr)
      ic = cregl0
      jc = cregl1
      kc = cregl2
      if = ic * ir
      jf = jc * jr
      kf = kc * kr
      center = 0.d0
      if (ivect(0) .eq. 0) then
c quadrants
c each quadrant is two octants and their share of the two central edges
         ffac = 2.d0 * ir * rm3
         cfac = 2.d0
         do ki = 0, 1
            do ji = 0, 1
               if (ga(0,ji,ki) .eq. 1) then
                  center = center + ffac
               else
                  center = center + cfac
               end if
            end do
         end do
         do ki = 0, 1
            do ji = 0, 1
               if (ga(0,ji,ki) - ga(0,ji,1-ki) .eq. 1) then
                  ffac = 2.d0 * ir * (jr - 1) * rm3
                  center = center + ffac
               end if
               if (ga(0,ji,ki) - ga(0,1-ji,ki) .eq. 1) then
                  ffac = 2.d0 * ir * (kr - 1) * rm3
                  center = center + ffac
               end if
            end do
         end do
         center = 1.d0 / center
         fac1 = center * rm3 / ir
         do ki = 0, 1
            do ji = 0, 1
               if (ga(0,ji,ki) .eq. 1) then
                  krf = kf + ki - 1
                  jrf = jf + ji - 1
                  do m = 0, ir-1
                     fac = (ir-m) * fac1
                     if (m .eq. 0) fac = 0.5d0 * fac
                     do if = ir*cregl0, ir*cregh0, ir
                        src(if,jf,kf) = src(if,jf,kf) + fac *
     &                    (rf(if-m,jrf,krf) + rf(if-m-1,jrf,krf) +
     &                     rf(if+m,jrf,krf) + rf(if+m-1,jrf,krf))
                     end do
                  end do
               else
                  krc = kc + ki - 1
                  jrc = jc + ji - 1
                  do ic = cregl0, cregh0
                     if = ic * ir
                     src(if,jf,kf) = src(if,jf,kf) + center *
     &                 (rc(ic,jrc,krc) + rc(ic-1,jrc,krc))
                  end do
               end if
            end do
         end do
c faces
c each face is two faces and two sides of an edge
         do ki = 0, 1
            kdir = 2 * ki - 1
            do ji = 0, 1
               jdir = 2 * ji - 1
               if (ga(0,ji,ki) - ga(0,ji,1-ki) .eq. 1) then
                  krf = kf + ki - 1
                  fac0 = center * rm3 / (ir * jr)
                  do n = jdir, jdir*(jr-1), jdir
                     fac1 = (jr-abs(n)) * fac0
                     do m = 0, ir-1
                        fac = (ir-m) * fac1
                        if (m .eq. 0) fac = 0.5d0 * fac
                        do if = ir*cregl0, ir*cregh0, ir
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                    (rf(if-m,jf+n,krf) + rf(if-m-1,jf+n,krf) +
     &                     rf(if-m,jf+n-1,krf) + rf(if-m-1,jf+n-1,krf) +
     &                     rf(if+m,jf+n,krf) + rf(if+m-1,jf+n,krf) +
     &                     rf(if+m,jf+n-1,krf) + rf(if+m-1,jf+n-1,krf))
                        end do
                     end do
                  end do
               end if
               if (ga(0,ji,ki) - ga(0,1-ji,ki) .eq. 1) then
                  jrf = jf + ji - 1
                  fac0 = center * rm3 / (ir * kr)
                  do l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do m = 0, ir-1
                        fac = (ir-m) * fac1
                        if (m .eq. 0) fac = 0.5d0 * fac
                        do if = ir*cregl0, ir*cregh0, ir
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                    (rf(if-m,jrf,kf+l) + rf(if-m-1,jrf,kf+l) +
     &                     rf(if-m,jrf,kf+l-1) + rf(if-m-1,jrf,kf+l-1) +
     &                     rf(if+m,jrf,kf+l) + rf(if+m-1,jrf,kf+l) +
     &                     rf(if+m,jrf,kf+l-1) + rf(if+m-1,jrf,kf+l-1))
                        end do
                     end do
                  end do
               end if
            end do
         end do
      else if (ivect(1) .eq. 0) then
c quadrants
c each quadrant is two octants and their share of the two central edges
         ffac = 2.d0 * jr * rm3
         cfac = 2.d0
         do ki = 0, 1
            do ii = 0, 1
               if (ga(ii,0,ki) .eq. 1) then
                  center = center + ffac
               else
                  center = center + cfac
               end if
            end do
         end do
         do ki = 0, 1
            do ii = 0, 1
               if (ga(ii,0,ki) - ga(ii,0,1-ki) .eq. 1) then
                  ffac = 2.d0 * jr * (ir - 1) * rm3
                  center = center + ffac
               end if
               if (ga(ii,0,ki) - ga(1-ii,0,ki) .eq. 1) then
                  ffac = 2.d0 * jr * (kr - 1) * rm3
                  center = center + ffac
               end if
            end do
         end do
         center = 1.d0 / center
         fac1 = center * rm3 / jr
         do ki = 0, 1
            do ii = 0, 1
               if (ga(ii,0,ki) .eq. 1) then
                  krf = kf + ki - 1
                  irf = if + ii - 1
                  do n = 0, jr-1
                     fac = (jr-n) * fac1
                     if (n .eq. 0) fac = 0.5d0 * fac
                     do jf = jr*cregl1, jr*cregh1, jr
                        src(if,jf,kf) = src(if,jf,kf) + fac *
     &                    (rf(irf,jf-n,krf) + rf(irf,jf-n-1,krf) +
     &                     rf(irf,jf+n,krf) + rf(irf,jf+n-1,krf))
                     end do
	            end do
               else
                  krc = kc + ki - 1
                  irc = ic + ii - 1
                  do jc = cregl1, cregh1
                     jf = jc * jr
                     src(if,jf,kf) = src(if,jf,kf) + center *
     &                 (rc(irc,jc,krc) + rc(irc,jc-1,krc))
                  end do
               end if
            end do
         end do
c faces
c each face is two faces and two sides of an edge
         do ki = 0, 1
            kdir = 2 * ki - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,0,ki) - ga(ii,0,1-ki) .eq. 1) then
                  krf = kf + ki - 1
                  fac0 = center * rm3 / (ir * jr)
                  do  n = 0, jr-1
                     fac1 = (jr-n) * fac0
                     if (n .eq. 0) fac1 = 0.5d0 * fac1
                     do m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        do jf = jr*cregl1, jr*cregh1, jr
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                   (rf(if+m,jf-n,krf) + rf(if+m-1,jf-n,krf) +
     &                    rf(if+m,jf-n-1,krf) + rf(if+m-1,jf-n-1,krf) +
     &                    rf(if+m,jf+n,krf) + rf(if+m-1,jf+n,krf) +
     &                    rf(if+m,jf+n-1,krf) + rf(if+m-1,jf+n-1,krf))
                        end do
                     end do
                  end do
               end if
               if (ga(ii,0,ki) - ga(1-ii,0,ki) .eq. 1) then
                  irf = if + ii - 1
                  fac0 = center * rm3 / (jr * kr)
                  do l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do n = 0, jr-1
                        fac = (jr-n) * fac1
                        if (n .eq. 0) fac = 0.5d0 * fac
                        do jf = jr*cregl1, jr*cregh1, jr
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                  (rf(irf,jf-n,kf+l) + rf(irf,jf-n-1,kf+l) +
     &                   rf(irf,jf-n,kf+l-1) + rf(irf,jf-n-1,kf+l-1) +
     &                   rf(irf,jf+n,kf+l) + rf(irf,jf+n-1,kf+l) +
     &                   rf(irf,jf+n,kf+l-1) + rf(irf,jf+n-1,kf+l-1))
                        end do
                     end do
                  end do
               end if
            end do
         end do
      else
c quadrants
c each quadrant is two octants and their share of the two central edges
         ffac = 2.d0 * kr * rm3
         cfac = 2.d0
         do ji = 0, 1
            do ii = 0, 1
               if (ga(ii,ji,0) .eq. 1) then
                  center = center + ffac
               else
                  center = center + cfac
               end if
            end do
         end do
         do ji = 0, 1
            do ii = 0, 1
               if (ga(ii,ji,0) - ga(ii,1-ji,0) .eq. 1) then
                  ffac = 2.d0 * kr * (ir - 1) * rm3
                  center = center + ffac
               end if
               if (ga(ii,ji,0) - ga(1-ii,ji,0) .eq. 1) then
                  ffac = 2.d0 * kr * (jr - 1) * rm3
                  center = center + ffac
               end if
            end do
         end do
         center = 1.d0 / center
         fac1 = center * rm3 / kr
         do ji = 0, 1
            do ii = 0, 1
               if (ga(ii,ji,0) .eq. 1) then
                  jrf = jf + ji - 1
                  irf = if + ii - 1
                  do l = 0, kr-1
                     fac = (kr-l) * fac1
                     if (l .eq. 0) fac = 0.5d0 * fac
                     do kf = kr*cregl2, kr*cregh2, kr
                        src(if,jf,kf) = src(if,jf,kf) + fac *
     &                    (rf(irf,jrf,kf-l) + rf(irf,jrf,kf-l-1) +
     &                     rf(irf,jrf,kf+l) + rf(irf,jrf,kf+l-1))
                     end do
                  end do
               else
                  jrc = jc + ji - 1
                  irc = ic + ii - 1
                  do kc = cregl2, cregh2
                     kf = kc * kr
                     src(if,jf,kf) = src(if,jf,kf) + center *
     &                 (rc(irc,jrc,kc) + rc(irc,jrc,kc-1))
                  end do
               end if
            end do
         end do
c faces
c each face is two faces and two sides of an edge
         do ji = 0, 1
            jdir = 2 * ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,0) - ga(ii,1-ji,0) .eq. 1) then
                  jrf = jf + ji - 1
                  fac0 = center * rm3 / (ir * kr)
                  do l = 0, kr-1
                     fac1 = (kr-l) * fac0
                     if (l .eq. 0) fac1 = 0.5d0 * fac1
                     do m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        do kf = kr*cregl2, kr*cregh2, kr
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                  (rf(if+m,jrf,kf-l) + rf(if+m-1,jrf,kf-l) +
     &                   rf(if+m,jrf,kf-l-1) + rf(if+m-1,jrf,kf-l-1) +
     &                   rf(if+m,jrf,kf+l) + rf(if+m-1,jrf,kf+l) +
     &                   rf(if+m,jrf,kf+l-1) + rf(if+m-1,jrf,kf+l-1))
                        end do
                     end do
                  end do
               end if
               if (ga(ii,ji,0) - ga(1-ii,ji,0) .eq. 1) then
                  irf = if + ii - 1
                  fac0 = center * rm3 / (jr * kr)
                  do l = 0, kr-1
                     fac1 = (kr-l) * fac0
                     if (l .eq. 0) fac1 = 0.5d0 * fac1
                     do n = jdir, jdir*(jr-1), jdir
                        fac = (jr-abs(n)) * fac1
                        do kf = kr*cregl2, kr*cregh2, kr
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                  (rf(irf,jf+n,kf-l) + rf(irf,jf+n-1,kf-l) +
     &                   rf(irf,jf+n,kf-l-1) + rf(irf,jf+n-1,kf-l-1) +
     &                   rf(irf,jf+n,kf+l) + rf(irf,jf+n-1,kf+l) +
     &                   rf(irf,jf+n,kf+l-1) + rf(irf,jf+n-1,kf+l-1))
                        end do
                     end do
                  end do
               end if
            end do
         end do
      end if
      end
c-----------------------------------------------------------------------
c Note---only generates values at coarse points along edge of fine grid
      subroutine hgcavg(src,
     & srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & rc,
     & cl0,ch0,cl1,ch1,cl2,ch2,
     & rf,
     & fl0,fh0,fl1,fh1,fl2,fh2,
     & cregl0,cregh0,cregl1,cregh1,cregl2,cregh2,
     & ir, jr, kr, ga)
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer cl0,ch0,cl1,ch1,cl2,ch2
      integer fl0,fh0,fl1,fh1,fl2,fh2
      integer cregl0,cregh0,cregl1,cregh1,cregl2,cregh2
      real*8 src(srcl0:srch0,srcl1:srch1,srcl2:srch2)
      real*8 rc(cl0:ch0,cl1:ch1,cl2:ch2)
      real*8 rf(fl0:fh0,fl1:fh1,fl2:fh2)
      integer ir, jr, kr, ga(0:1,0:1,0:1)
      real*8 rm3, sum, center, cfac, ffac, fac0, fac1, fac
      integer ic, jc, kc, if, jf, kf, ii, ji, ki
      integer irc, irf, jrc, jrf, krc, krf, idir, jdir, kdir, l, m, n
      rm3 = 1.d0 / (ir * jr * kr)
      ic = cregl0
      jc = cregl1
      kc = cregl2
      if = ic * ir
      jf = jc * jr
      kf = kc * kr
      sum = 0.d0
      center = 0.d0
c octants
      fac = rm3
      ffac = rm3
      cfac = 1.d0
      do ki = 0, 1
         krf = kf + ki - 1
         krc = kc + ki - 1
         do ji = 0, 1
            jrf = jf + ji - 1
            jrc = jc + ji - 1
            do ii = 0, 1
               if (ga(ii,ji,ki) .eq. 1) then
                  irf = if + ii - 1
                  center = center + ffac
                  sum = sum + fac * rf(irf,jrf,krf)
               else
                  irc = ic + ii - 1
                  center = center + cfac
                  sum = sum + rc(irc,jrc,krc)
               end if
            end do
         end do
      end do
c faces
      do ki = 0, 1
         kdir = 2 * ki - 1
         krf = kf + ki - 1
         do ji = 0, 1
            jdir = 2 * ji - 1
            jrf = jf + ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               irf = if + ii - 1
               if (ga(ii,ji,ki) - ga(ii,ji,1-ki) .eq. 1) then
                  fac0 = rm3 / (ir * jr)
                  ffac = (ir-1) * (jr-1) * rm3
                  center = center + ffac
                  do n = jdir, jdir*(jr-1), jdir
                     fac1 = (jr-abs(n)) * fac0
                     do m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        sum = sum + fac *
     &                    (rf(if+m,jf+n,krf) + rf(if+m-1,jf+n,krf) +
     &                     rf(if+m,jf+n-1,krf) + rf(if+m-1,jf+n-1,krf))
                     end do
                  end do
               end if
               if (ga(ii,ji,ki) - ga(ii,1-ji,ki) .eq. 1) then
                  fac0 = rm3 / (ir * kr)
                  ffac = (ir-1) * (kr-1) * rm3
                  center = center + ffac
                  do l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        sum = sum + fac *
     &                    (rf(if+m,jrf,kf+l) + rf(if+m-1,jrf,kf+l) +
     &                     rf(if+m,jrf,kf+l-1) + rf(if+m-1,jrf,kf+l-1))
                     end do
                  end do
               end if
               if (ga(ii,ji,ki) - ga(1-ii,ji,ki) .eq. 1) then
                  fac0 = rm3 / (jr * kr)
                  ffac = (jr-1) * (kr-1) * rm3
                  center = center + ffac
                  do l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do n = jdir, jdir*(jr-1), jdir
                        fac = (jr-abs(n)) * fac1
                        sum = sum + fac *
     &                    (rf(irf,jf+n,kf+l) + rf(irf,jf+n,kf+l-1) +
     &                     rf(irf,jf+n-1,kf+l) + rf(irf,jf+n-1,kf+l-1))
                     end do
                  end do
               end if
            end do
         end do
      end do
c edges
      do ki = 0, 1
         kdir = 2 * ki - 1
         krf = kf + ki - 1
         do ji = 0, 1
            jdir = 2 * ji - 1
            jrf = jf + ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               irf = if + ii - 1
               if (ga(ii,ji,ki) -
     &             min(ga(ii,ji,1-ki), ga(ii,1-ji,ki), ga(ii,1-ji,1-ki))
     &             .eq. 1) then
                  fac1 = rm3 / ir
                  ffac = (ir-1) * rm3
                  center = center + ffac
                  do m = idir, idir*(ir-1), idir
                     fac = (ir-abs(m)) * fac1
                     sum = sum + fac *
     &                 (rf(if+m,jrf,krf) + rf(if+m-1,jrf,krf))
                  end do
               end if
               if (ga(ii,ji,ki) -
     &             min(ga(ii,ji,1-ki), ga(1-ii,ji,ki), ga(1-ii,ji,1-ki))
     &             .eq. 1) then
                  fac1 = rm3 / jr
                  ffac = (jr-1) * rm3
                  center = center + ffac
                  do n = jdir, jdir*(jr-1), jdir
                     fac = (jr-abs(n)) * fac1
                     sum = sum + fac *
     &                 (rf(irf,jf+n,krf) + rf(irf,jf+n-1,krf))
                  end do
               end if
               if (ga(ii,ji,ki) -
     &             min(ga(ii,1-ji,ki), ga(1-ii,ji,ki), ga(1-ii,1-ji,ki))
     &             .eq. 1) then
                  fac1 = rm3 / kr
                  ffac = (kr-1) * rm3
                  center = center + ffac
                  do l = kdir, kdir*(kr-1), kdir
                     fac = (kr-abs(l)) * fac1
                     sum = sum + fac *
     &                 (rf(irf,jrf,kf+l) + rf(irf,jrf,kf+l-1))
                  end do
               end if
            end do
         end do
      end do
c weighting
      src(if,jf,kf) = src(if,jf,kf) + sum / center
      end
