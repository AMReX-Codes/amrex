c-----------------------------------------------------------------------
      subroutine hggrad_terrain(
     & gpx, gpy, gpz,
     &       gpl0,gph0,gpl1,gph1,gpl2,gph2,
     & dest, destl0,desth0,destl1,desth1,destl2,desth2,
     &       fregl0,fregh0,fregl1,fregh1,fregl2,fregh2,
     & hx, hy, hz)
      integer gpl0,gph0,gpl1,gph1,gpl2,gph2
      integer destl0,desth0,destl1,desth1,destl2,desth2
      integer fregl0,fregh0,fregl1,fregh1,fregl2,fregh2
	double precision hx, hy, hz
      double precision gpx(gpl0:gph0,gpl1:gph1,gpl2:gph2)
      double precision gpy(gpl0:gph0,gpl1:gph1,gpl2:gph2)
      double precision gpz(gpl0:gph0,gpl1:gph1,gpl2:gph2)
      double precision dest(destl0:desth0,destl1:desth1,destl2:desth2)
      integer i, j, k
      do k = fregl2, fregh2
          do j = fregl1, fregh1
              do i = fregl0, fregh0
               gpx(i,j,k) = 0.25d0 * 
     &            (dest(i+1,j,k  ) + dest(i+1,j+1,k  ) +
     &             dest(i+1,j,k+1) + dest(i+1,j+1,k+1) -
     &             dest(i  ,j,k  ) - dest(i  ,j+1,k  ) -
     &             dest(i  ,j,k+1) - dest(i  ,j+1,k+1))
               gpy(i,j,k) = 0.25d0 * 
     &            (dest(i,j+1,k  ) + dest(i+1,j+1,k  ) +
     &             dest(i,j+1,k+1) + dest(i+1,j+1,k+1) -
     &             dest(i,j  ,k  ) - dest(i+1,j  ,k  ) -
     &             dest(i,j  ,k+1) - dest(i+1,j  ,k+1))
               gpz(i,j,k) = 0.25d0 * 
     &            (dest(i,j  ,k+1) + dest(i+1,j  ,k+1) +
     &             dest(i,j+1,k+1) + dest(i+1,j+1,k+1) -
     &             dest(i,j  ,k  ) - dest(i+1,j  ,k  ) -
     &             dest(i,j+1,k  ) - dest(i+1,j+1,k  ))
              end do
          end do
      end do
      end
c-----------------------------------------------------------------------
      subroutine hgdiv_terrain(
     & src, srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & uf, vf, wf,
     &      fl0,fh0,fl1,fh1,fl2,fh2,
     &      fregl0,fregh0,fregl1,fregh1,fregl2,fregh2,
     & hx, hy, hz)
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer fl0,fh0,fl1,fh1,fl2,fh2
      integer fregl0,fregh0,fregl1,fregh1,fregl2,fregh2
	double precision hx, hy, hz
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2)
      double precision uf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision vf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision wf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision fac
      integer i, j, k
      fac = 0.25d0
      do k = fregl2, fregh2
         do j = fregl1, fregh1
            do i = fregl0, fregh0
             src(i,j,k) = fac *
     &              (uf(i  ,j-1,k-1) - uf(i-1,j-1,k-1) +
     &               uf(i  ,j-1,k  ) - uf(i-1,j-1,k  ) +
     &               uf(i  ,j  ,k-1) - uf(i-1,j  ,k-1) +
     &               uf(i  ,j  ,k  ) - uf(i-1,j  ,k  ) +
     &               vf(i-1,j  ,k-1) - vf(i-1,j-1,k-1) +
     &               vf(i-1,j  ,k  ) - vf(i-1,j-1,k  ) +
     &               vf(i  ,j  ,k-1) - vf(i  ,j-1,k-1) +
     &               vf(i  ,j  ,k  ) - vf(i  ,j-1,k  ) +
     &               wf(i-1,j-1,k  ) - wf(i-1,j-1,k-1) +
     &               wf(i-1,j  ,k  ) - wf(i-1,j  ,k-1) +
     &               wf(i  ,j-1,k  ) - wf(i  ,j-1,k-1) +
     &               wf(i  ,j  ,k  ) - wf(i  ,j  ,k-1))
              end do
          end do
      end do
      end
c-----------------------------------------------------------------------
c Note---only generates values at coarse points along edge of fine grid
      subroutine hgfdiv_terrain(
     & src, srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & uc, vc, wc,
     &      cl0,ch0,cl1,ch1,cl2,ch2,
     & uf, vf, wf,
     &      fl0,fh0,fl1,fh1,fl2,fh2,
     &      cregl0,cregh0,cregl1,cregh1,cregl2,cregh2,
     & hx, hy, hz, ir, jr, kr, idim, idir)
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer cl0,ch0,cl1,ch1,cl2,ch2
      integer fl0,fh0,fl1,fh1,fl2,fh2
      integer cregl0,cregh0,cregl1,cregh1,cregl2,cregh2
      double precision hx, hy, hz
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2)
      double precision uc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision vc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision wc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision uf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision vf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision wf(fl0:fh0,fl1:fh1,fl2:fh2)
      integer ir, jr, kr, idim, idir
      double precision fac0, fac1, fac, tmp
      integer i, j, k, iuc, iuf, juc, juf, kuc, kuf, l, m, n
      if (idim .eq. 0) then
         i = cregl0
         if (idir .eq. 1) then
            iuc = i - 1
            iuf = i * ir
         else
            iuc = i
            iuf = i * ir - 1
         end if
         fac0 = 0.25d0
         do  k = cregl2, cregh2
            do  j = cregl1, cregh1
             src(i*ir,j*jr,k*kr) = fac0 *
     &           ((vc(iuc,j,k) - vc(iuc,j-1,k) +
     &             vc(iuc,j,k-1) - vc(iuc,j-1,k-1)) +
     &            (wc(iuc,j,k) - wc(iuc,j,k-1) +
     &             wc(iuc,j-1,k) - wc(iuc,j-1,k-1)) -
     &     idir * (uc(iuc,j,k) + uc(iuc,j,k-1) +
     &             uc(iuc,j-1,k) + uc(iuc,j-1,k-1)))
              end do
          end do
         fac0 = fac0 / (jr * kr)
         i = i * ir
         do l = 0, kr-1
            fac1 = (kr-l) * fac0
            if (l .eq. 0) fac1 = 0.5d0 * fac1
            do n = 0, jr-1
               fac = (jr-n) * fac1
               if (n .eq. 0) fac = 0.5d0 * fac
               do k = kr*cregl2, kr*cregh2, kr
                  do j = jr*cregl1, jr*cregh1, jr
                     tmp = idir *
     &                 (uf(iuf,j-n  ,k-l) + uf(iuf,j-n  ,k-l-1) +
     &                  uf(iuf,j-n-1,k-l) + uf(iuf,j-n-1,k-l-1) +
     &                  uf(iuf,j-n  ,k+l) + uf(iuf,j-n  ,k+l-1) +
     &                  uf(iuf,j-n-1,k+l) + uf(iuf,j-n-1,k+l-1) +
     &                  uf(iuf,j+n  ,k-l) + uf(iuf,j+n  ,k-l-1) +
     &                  uf(iuf,j+n-1,k-l) + uf(iuf,j+n-1,k-l-1) +
     &                  uf(iuf,j+n  ,k+l) + uf(iuf,j+n  ,k+l-1) +
     &                  uf(iuf,j+n-1,k+l) + uf(iuf,j+n-1,k+l-1))
                     tmp = tmp + 
     &                 (vf(iuf,j-n,k-l  ) - vf(iuf,j-n-1,k-l  ) +
     &                  vf(iuf,j-n,k-l-1) - vf(iuf,j-n-1,k-l-1) +
     &                  vf(iuf,j-n,k+l  ) - vf(iuf,j-n-1,k+l  ) +
     &                  vf(iuf,j-n,k+l-1) - vf(iuf,j-n-1,k+l-1) +
     &                  vf(iuf,j+n,k-l  ) - vf(iuf,j+n-1,k-l  ) +
     &                  vf(iuf,j+n,k-l-1) - vf(iuf,j+n-1,k-l-1) +
     &                  vf(iuf,j+n,k+l  ) - vf(iuf,j+n-1,k+l  ) +
     &                  vf(iuf,j+n,k+l-1) - vf(iuf,j+n-1,k+l-1))
                     src(i,j,k) = src(i,j,k) + fac * (tmp + 
     &                 (wf(iuf,j-n  ,k-l) - wf(iuf,j-n  ,k-l-1) +
     &                  wf(iuf,j-n-1,k-l) - wf(iuf,j-n-1,k-l-1) +
     &                  wf(iuf,j-n  ,k+l) - wf(iuf,j-n  ,k+l-1) +
     &                  wf(iuf,j-n-1,k+l) - wf(iuf,j-n-1,k+l-1) +
     &                  wf(iuf,j+n  ,k-l) - wf(iuf,j+n  ,k-l-1) +
     &                  wf(iuf,j+n-1,k-l) - wf(iuf,j+n-1,k-l-1) +
     &                  wf(iuf,j+n  ,k+l) - wf(iuf,j+n  ,k+l-1) +
     &                  wf(iuf,j+n-1,k+l) - wf(iuf,j+n-1,k+l-1)))
                  end do
               end do
            end do
         end do
      else if (idim .eq. 1) then
         j = cregl1
         if (idir .eq. 1) then
            juc = j - 1
            juf = j * jr
         else
            juc = j
            juf = j * jr - 1
         end if
         fac0 = 0.25d0
         do k = cregl2, cregh2
            do i = cregl0, cregh0
               src(i*ir,j*jr,k*kr) = fac0 *
     &           ((uc(i,juc,k) - uc(i-1,juc,k) +
     &             uc(i,juc,k-1) - uc(i-1,juc,k-1)) -
     &     idir * (vc(i,juc,k) + vc(i,juc,k-1) +
     &             vc(i-1,juc,k) + vc(i-1,juc,k-1)) +
     &            (wc(i,juc,k) - wc(i,juc,k-1) +
     &             wc(i-1,juc,k) - wc(i-1,juc,k-1)))
            end do
         end do
         fac0 = fac0 / (ir * kr)
         j = j * jr
         do l = 0, kr-1
            fac1 = (kr-l) * fac0
            if (l .eq. 0) fac1 = 0.5d0 * fac1
            do m = 0, ir-1
               fac = (ir-m) * fac1
               if (m .eq. 0) fac = 0.5d0 * fac
               do k = kr*cregl2, kr*cregh2, kr
                  do i = ir*cregl0, ir*cregh0, ir
                     tmp = 
     &                 (uf(i-m,juf,k-l  ) - uf(i-m-1,juf,k-l) +
     &                  uf(i-m,juf,k-l-1) - uf(i-m-1,juf,k-l-1) +
     &                  uf(i-m,juf,k+l  ) - uf(i-m-1,juf,k+l) +
     &                  uf(i-m,juf,k+l-1) - uf(i-m-1,juf,k+l-1) +
     &                  uf(i+m,juf,k-l  ) - uf(i+m-1,juf,k-l) +
     &                  uf(i+m,juf,k-l-1) - uf(i+m-1,juf,k-l-1) +
     &                  uf(i+m,juf,k+l  ) - uf(i+m-1,juf,k+l) +
     &                  uf(i+m,juf,k+l-1) - uf(i+m-1,juf,k+l-1))
                     tmp = tmp + idir *
     &                 (vf(i-m  ,juf,k-l) + vf(i-m  ,juf,k-l-1) +
     &                  vf(i-m-1,juf,k-l) + vf(i-m-1,juf,k-l-1) +
     &                  vf(i-m  ,juf,k+l) + vf(i-m  ,juf,k+l-1) +
     &                  vf(i-m-1,juf,k+l) + vf(i-m-1,juf,k+l-1) +
     &                  vf(i+m  ,juf,k-l) + vf(i+m  ,juf,k-l-1) +
     &                  vf(i+m-1,juf,k-l) + vf(i+m-1,juf,k-l-1) +
     &                  vf(i+m  ,juf,k+l) + vf(i+m  ,juf,k+l-1) +
     &                  vf(i+m-1,juf,k+l) + vf(i+m-1,juf,k+l-1))
                     src(i,j,k) = src(i,j,k) + fac * (tmp + 
     &                 (wf(i-m  ,juf,k-l) - wf(i-m  ,juf,k-l-1) +
     &                  wf(i-m-1,juf,k-l) - wf(i-m-1,juf,k-l-1) +
     &                  wf(i-m  ,juf,k+l) - wf(i-m  ,juf,k+l-1) +
     &                  wf(i-m-1,juf,k+l) - wf(i-m-1,juf,k+l-1) +
     &                  wf(i+m  ,juf,k-l) - wf(i+m  ,juf,k-l-1) +
     &                  wf(i+m-1,juf,k-l) - wf(i+m-1,juf,k-l-1) +
     &                  wf(i+m  ,juf,k+l) - wf(i+m  ,juf,k+l-1) +
     &                  wf(i+m-1,juf,k+l) - wf(i+m-1,juf,k+l-1)))
                  end do
               end do
            end do
         end do
      else
         k = cregl2
         if (idir .eq. 1) then
            kuc = k - 1
            kuf = k * kr
         else
            kuc = k
            kuf = k * kr - 1
         end if
         fac0 = 0.25d0
         do j = cregl1, cregh1
            do i = cregl0, cregh0
               src(i*ir,j*jr,k*kr) = fac0 *
     &           ((uc(i,j,kuc) - uc(i-1,j,kuc) +
     &             uc(i,j-1,kuc) - uc(i-1,j-1,kuc)) +
     &            (vc(i,j,kuc) - vc(i,j-1,kuc) +
     &             vc(i-1,j,kuc) - vc(i-1,j-1,kuc)) -
     &     idir * (wc(i,j,kuc) + wc(i,j-1,kuc) +
     &             wc(i-1,j,kuc) + wc(i-1,j-1,kuc)))
            end do
         end do
         fac0 = fac0 / (ir * jr)
         k = k * kr
         do n = 0, jr-1
            fac1 = (jr-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do m = 0, ir-1
               fac = (ir-m) * fac1
               if (m .eq. 0) fac = 0.5d0 * fac
               do j = jr*cregl1, jr*cregh1, jr
                  do i = ir*cregl0, ir*cregh0, ir
                     tmp = 
     &                 (uf(i-m,j-n  ,kuf) - uf(i-m-1,j-n  ,kuf) +
     &                  uf(i-m,j-n-1,kuf) - uf(i-m-1,j-n-1,kuf) +
     &                  uf(i-m,j+n  ,kuf) - uf(i-m-1,j+n  ,kuf) +
     &                  uf(i-m,j+n-1,kuf) - uf(i-m-1,j+n-1,kuf) +
     &                  uf(i+m,j-n  ,kuf) - uf(i+m-1,j-n  ,kuf) +
     &                  uf(i+m,j-n-1,kuf) - uf(i+m-1,j-n-1,kuf) +
     &                  uf(i+m,j+n  ,kuf) - uf(i+m-1,j+n  ,kuf) +
     &                  uf(i+m,j+n-1,kuf) - uf(i+m-1,j+n-1,kuf))
                     tmp = tmp + 
     &                 (vf(i-m  ,j-n,kuf) - vf(i-m  ,j-n-1,kuf) +
     &                  vf(i-m-1,j-n,kuf) - vf(i-m-1,j-n-1,kuf) +
     &                  vf(i-m  ,j+n,kuf) - vf(i-m  ,j+n-1,kuf) +
     &                  vf(i-m-1,j+n,kuf) - vf(i-m-1,j+n-1,kuf) +
     &                  vf(i+m  ,j-n,kuf) - vf(i+m  ,j-n-1,kuf) +
     &                  vf(i+m-1,j-n,kuf) - vf(i+m-1,j-n-1,kuf) +
     &                  vf(i+m  ,j+n,kuf) - vf(i+m  ,j+n-1,kuf) +
     &                  vf(i+m-1,j+n,kuf) - vf(i+m-1,j+n-1,kuf))
                     src(i,j,k) = src(i,j,k) + fac * (tmp + idir *
     &                 (wf(i-m  ,j-n,kuf) + wf(i-m  ,j-n-1,kuf) +
     &                  wf(i-m-1,j-n,kuf) + wf(i-m-1,j-n-1,kuf) +
     &                  wf(i-m  ,j+n,kuf) + wf(i-m  ,j+n-1,kuf) +
     &                  wf(i-m-1,j+n,kuf) + wf(i-m-1,j+n-1,kuf) +
     &                  wf(i+m  ,j-n,kuf) + wf(i+m  ,j-n-1,kuf) +
     &                  wf(i+m-1,j-n,kuf) + wf(i+m-1,j-n-1,kuf) +
     &                  wf(i+m  ,j+n,kuf) + wf(i+m  ,j+n-1,kuf) +
     &                  wf(i+m-1,j+n,kuf) + wf(i+m-1,j+n-1,kuf)))
                  end do
               end do
            end do
         end do
      end if
      end
c-----------------------------------------------------------------------
c Note---only generates values at coarse points along edge of fine grid
      subroutine hgediv_terrain(
     & src, srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & uc, vc, wc,
     &      cl0,ch0,cl1,ch1,cl2,ch2,
     & uf, vf, wf,
     &     fl0,fh0,fl1,fh1,fl2,fh2,
     &     cregl0,cregh0,cregl1,cregh1,cregl2,cregh2,
     & hx, hy, hz, ir, jr, kr, ivect, ga)
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer cl0,ch0,cl1,ch1,cl2,ch2
      integer fl0,fh0,fl1,fh1,fl2,fh2
      integer cregl0,cregh0,cregl1,cregh1,cregl2,cregh2
      double precision hx, hy, hz
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2)
      double precision uc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision vc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision wc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision uf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision vf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision wf(fl0:fh0,fl1:fh1,fl2:fh2)
      integer ir, jr, kr, ivect(0:2), ga(0:1,0:1,0:1)
      double precision fac0, fac1, fac
      integer ic, jc, kc, if, jf, kf, iuc, iuf, juc, juf, kuc, kuf
      integer ii, ji, ki, idir, jdir, kdir, l, m, n
      ic = cregl0
      jc = cregl1
      kc = cregl2
      if = ic * ir
      jf = jc * jr
      kf = kc * kr
      if (ivect(0) .eq. 0) then
         do if = ir*cregl0, ir*cregh0, ir
            src(if,jf,kf) = 0.d0
         end do
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac1 = 1.d0 / ir
         do ki = 0, 1
            kdir = 2 * ki - 1
            do ji = 0, 1
               jdir = 2 * ji - 1
               if (ga(0,ji,ki) .eq. 1) then
                  kuf = kf + ki - 1
                  juf = jf + ji - 1
                  do m = 0, ir-1
                     fac = (ir-m) * fac1
                     if (m .eq. 0) fac = 0.5d0 * fac
                     do if = ir*cregl0, ir*cregh0, ir
                        src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                 (
     &                    (uf(if-m,juf,kuf) - uf(if-m-1,juf,kuf) +
     &                     uf(if+m,juf,kuf) - uf(if+m-1,juf,kuf)) +
     &                                  jdir *
     &                    (vf(if-m,juf,kuf) + vf(if-m-1,juf,kuf) +
     &                     vf(if+m,juf,kuf) + vf(if+m-1,juf,kuf)) +
     &                                  kdir *
     &                    (wf(if-m,juf,kuf) + wf(if-m-1,juf,kuf) +
     &                     wf(if+m,juf,kuf) + wf(if+m-1,juf,kuf)))
                     end do
                  end do
               else
                  kuc = kc + ki - 1
                  juc = jc + ji - 1
                  do ic = cregl0, cregh0
                     if = ic * ir
                     src(if,jf,kf) = src(if,jf,kf) + 
     &                 (       (uc(ic,juc,kuc) - uc(ic-1,juc,kuc)) +
     &                  jdir * (vc(ic,juc,kuc) + vc(ic-1,juc,kuc)) +
     &                  kdir * (wc(ic,juc,kuc) + wc(ic-1,juc,kuc)))
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
                  kuf = kf + ki - 1
                  fac0 = 1.d0 / (ir * jr)
                  do n = jdir, jdir*(jr-1), jdir
                     fac1 = (jr-abs(n)) * fac0
                     do m = 0, ir-1
                        fac = (ir-m) * fac1
                        if (m .eq. 0) fac = 0.5d0 * fac
                        do if = ir*cregl0, ir*cregh0, ir
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                       (
     &            (uf(if-m  ,jf+n  ,kuf) - uf(if-m-1,jf+n  ,kuf) +
     &             uf(if-m  ,jf+n-1,kuf) - uf(if-m-1,jf+n-1,kuf) +
     &             uf(if+m  ,jf+n  ,kuf) - uf(if+m-1,jf+n  ,kuf) +
     &             uf(if+m  ,jf+n-1,kuf) - uf(if+m-1,jf+n-1,kuf)) +
c
     &            (vf(if-m  ,jf+n  ,kuf) - vf(if-m  ,jf+n-1,kuf) +
     &             vf(if-m-1,jf+n  ,kuf) - vf(if-m-1,jf+n-1,kuf) +
     &             vf(if+m  ,jf+n  ,kuf) - vf(if+m  ,jf+n-1,kuf) +
     &             vf(if+m-1,jf+n  ,kuf) - vf(if+m-1,jf+n-1,kuf)) +
     &                                         kdir *
     &            (wf(if-m  ,jf+n,kuf) + wf(if-m  ,jf+n-1,kuf) +
     &             wf(if-m-1,jf+n,kuf) + wf(if-m-1,jf+n-1,kuf) +
     &             wf(if+m  ,jf+n,kuf) + wf(if+m  ,jf+n-1,kuf) +
     &             wf(if+m-1,jf+n,kuf) + wf(if+m-1,jf+n-1,kuf)))
                        end do
                     end do
                  end do
               end if
               if (ga(0,ji,ki) - ga(0,1-ji,ki) .eq. 1) then
                  juf = jf + ji - 1
                  fac0 = 1.d0 / (ir * kr)
                  do l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do m = 0, ir-1
                        fac = (ir-m) * fac1
                        if (m .eq. 0) fac = 0.5d0 * fac
                        do if = ir*cregl0, ir*cregh0, ir
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                       (
     &            (uf(if-m,juf,kf+l  ) - uf(if-m-1,juf,kf+l  ) +
     &             uf(if-m,juf,kf+l-1) - uf(if-m-1,juf,kf+l-1) +
     &             uf(if+m,juf,kf+l  ) - uf(if+m-1,juf,kf+l  ) +
     &             uf(if+m,juf,kf+l-1) - uf(if+m-1,juf,kf+l-1)) +
     &                                        jdir *
     &            (vf(if-m  ,juf,kf+l) + vf(if-m  ,juf,kf+l-1) +
     &             vf(if-m-1,juf,kf+l) + vf(if-m-1,juf,kf+l-1) +
     &             vf(if+m  ,juf,kf+l) + vf(if+m  ,juf,kf+l-1) +
     &             vf(if+m-1,juf,kf+l) + vf(if+m-1,juf,kf+l-1)) +
c
     &            (wf(if-m  ,juf,kf+l) - wf(if-m  ,juf,kf+l-1) +
     &             wf(if-m-1,juf,kf+l) - wf(if-m-1,juf,kf+l-1) +
     &             wf(if+m  ,juf,kf+l) - wf(if+m  ,juf,kf+l-1) +
     &             wf(if+m-1,juf,kf+l) - wf(if+m-1,juf,kf+l-1)))
                        end do
                     end do
                  end do
               end if
            end do
         end do
c weighting
         do if = ir*cregl0, ir*cregh0, ir
            src(if,jf,kf) = .25d0 * src(if,jf,kf)
         end do

      else if (ivect(1) .eq. 0) then
         do jf = jr*cregl1, jr*cregh1, jr
            src(if,jf,kf) = 0.d0
         end do
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac1 = 1.d0 / jr
         do ki = 0, 1
            kdir = 2 * ki - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,0,ki) .eq. 1) then
                  kuf = kf + ki - 1
                  iuf = if + ii - 1
                  do  n = 0, jr-1
                     fac = (jr-n) * fac1
                     if (n .eq. 0) fac = 0.5d0 * fac
                     do jf = jr*cregl1, jr*cregh1, jr
                        src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                 (idir *
     &                    (uf(iuf,jf-n,kuf) + uf(iuf,jf-n-1,kuf) +
     &                     uf(iuf,jf+n,kuf) + uf(iuf,jf+n-1,kuf)) +
c
     &                    (vf(iuf,jf-n,kuf) - vf(iuf,jf-n-1,kuf) +
     &                     vf(iuf,jf+n,kuf) - vf(iuf,jf+n-1,kuf)) +
     &                                  kdir *
     &                    (wf(iuf,jf-n,kuf) + wf(iuf,jf-n-1,kuf) +
     &                     wf(iuf,jf+n,kuf) + wf(iuf,jf+n-1,kuf)))
                     end do
                  end do
               else
                  kuc = kc + ki - 1
                  iuc = ic + ii - 1
                  do jc = cregl1, cregh1
                     jf = jc * jr
                     src(if,jf,kf) = src(if,jf,kf) + 
     &                 (idir * (uc(iuc,jc,kuc) + uc(iuc,jc-1,kuc)) +
     &                         (vc(iuc,jc,kuc) - vc(iuc,jc-1,kuc)) +
     &                  kdir * (wc(iuc,jc,kuc) + wc(iuc,jc-1,kuc)))
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
                  kuf = kf + ki - 1
                  fac0 = 1.d0 / (ir * jr)
                  do n = 0, jr-1
                     fac1 = (jr-n) * fac0
                     if (n .eq. 0) fac1 = 0.5d0 * fac1
                     do m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        do jf = jr*cregl1, jr*cregh1, jr
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                       (
     &            (uf(if+m  ,jf-n  ,kuf) - uf(if+m-1,jf-n  ,kuf) +
     &             uf(if+m  ,jf-n-1,kuf) - uf(if+m-1,jf-n-1,kuf) +
     &             uf(if+m  ,jf+n  ,kuf) - uf(if+m-1,jf+n  ,kuf) +
     &             uf(if+m  ,jf+n-1,kuf) - uf(if+m-1,jf+n-1,kuf)) +
     &            (vf(if+m  ,jf-n  ,kuf) - vf(if+m  ,jf-n-1,kuf) +
     &             vf(if+m-1,jf-n  ,kuf) - vf(if+m-1,jf-n-1,kuf) +
     &             vf(if+m  ,jf+n  ,kuf) - vf(if+m  ,jf+n-1,kuf) +
     &             vf(if+m-1,jf+n  ,kuf) - vf(if+m-1,jf+n-1,kuf)) +
     &                                        kdir *
     &            (wf(if+m  ,jf-n,kuf) + wf(if+m  ,jf-n-1,kuf) +
     &             wf(if+m-1,jf-n,kuf) + wf(if+m-1,jf-n-1,kuf) +
     &             wf(if+m  ,jf+n,kuf) + wf(if+m  ,jf+n-1,kuf) +
     &             wf(if+m-1,jf+n,kuf) + wf(if+m-1,jf+n-1,kuf)))
                        end do
                     end do
                  end do
               end if
               if (ga(ii,0,ki) - ga(1-ii,0,ki) .eq. 1) then
                  iuf = if + ii - 1
                  fac0 = 1.d0 / (jr * kr)
                  do l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do n = 0, jr-1
                        fac = (jr-n) * fac1
                        if (n .eq. 0) fac = 0.5d0 * fac
                        do jf = jr*cregl1, jr*cregh1, jr
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                       (idir *
     &            (uf(iuf,jf-n  ,kf+l  ) + uf(iuf,jf-n-1,kf+l  ) +
     &             uf(iuf,jf-n  ,kf+l-1) + uf(iuf,jf-n-1,kf+l-1) +
     &             uf(iuf,jf+n  ,kf+l  ) + uf(iuf,jf+n-1,kf+l  ) +
     &             uf(iuf,jf+n  ,kf+l-1) + uf(iuf,jf+n-1,kf+l-1)) +
     &            (vf(iuf,jf-n  ,kf+l  ) - vf(iuf,jf-n-1,kf+l  ) +
     &             vf(iuf,jf-n  ,kf+l-1) - vf(iuf,jf-n-1,kf+l-1) +
     &             vf(iuf,jf+n  ,kf+l  ) - vf(iuf,jf+n-1,kf+l  ) +
     &             vf(iuf,jf+n  ,kf+l-1) - vf(iuf,jf+n-1,kf+l-1)) +
     &            (wf(iuf,jf-n  ,kf+l  ) - wf(iuf,jf-n  ,kf+l-1) +
     &             wf(iuf,jf-n-1,kf+l  ) - wf(iuf,jf-n-1,kf+l-1) +
     &             wf(iuf,jf+n  ,kf+l  ) - wf(iuf,jf+n  ,kf+l-1) +
     &             wf(iuf,jf+n-1,kf+l  ) - wf(iuf,jf+n-1,kf+l-1)))
                        end do
                     end do
                  end do
               end if
            end do
         end do
c weighting
         do jf = jr*cregl1, jr*cregh1, jr
            src(if,jf,kf) = .25d0 * src(if,jf,kf)
         end do

      else
         do kf = kr*cregl2, kr*cregh2, kr
            src(if,jf,kf) = 0.d0
         end do
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac1 = 1.d0 / kr
         do ji = 0, 1
            jdir = 2 * ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,0) .eq. 1) then
                  juf = jf + ji - 1
                  iuf = if + ii - 1
                  do l = 0, kr-1
                     fac = (kr-l) * fac1
                     if (l .eq. 0) fac = 0.5d0 * fac
                     do kf = kr*cregl2, kr*cregh2, kr
                        src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                 (idir *
     &                    (uf(iuf,juf,kf-l) + uf(iuf,juf,kf-l-1) +
     &                     uf(iuf,juf,kf+l) + uf(iuf,juf,kf+l-1)) +
     &                                  jdir *
     &                    (vf(iuf,juf,kf-l) + vf(iuf,juf,kf-l-1) +
     &                     vf(iuf,juf,kf+l) + vf(iuf,juf,kf+l-1)) +
c
     &                    (wf(iuf,juf,kf-l) - wf(iuf,juf,kf-l-1) +
     &                     wf(iuf,juf,kf+l) - wf(iuf,juf,kf+l-1)))
                     end do
                  end do
               else
                  juc = jc + ji - 1
                  iuc = ic + ii - 1
                  do kc = cregl2, cregh2
                     kf = kc * kr
                     src(if,jf,kf) = src(if,jf,kf) + 
     &                 (idir * (uc(iuc,juc,kc) + uc(iuc,juc,kc-1)) +
     &                  jdir * (vc(iuc,juc,kc) + vc(iuc,juc,kc-1)) +
     &                         (wc(iuc,juc,kc) - wc(iuc,juc,kc-1)))
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
                  juf = jf + ji - 1
                  fac0 = 1.d0 / (ir * kr)
                  do l = 0, kr-1
                     fac1 = (kr-l) * fac0
                     if (l .eq. 0) fac1 = 0.5d0 * fac1
                     do m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        do kf = kr*cregl2, kr*cregh2, kr
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                       (
     &            (uf(if+m,juf,kf-l  ) - uf(if+m-1,juf,kf-l) +
     &             uf(if+m,juf,kf-l-1) - uf(if+m-1,juf,kf-l-1) +
     &             uf(if+m,juf,kf+l  ) - uf(if+m-1,juf,kf+l) +
     &             uf(if+m,juf,kf+l-1) - uf(if+m-1,juf,kf+l-1)) +
     &                                        jdir *
     &            (vf(if+m  ,juf,kf-l) + vf(if+m  ,juf,kf-l-1) +
     &             vf(if+m-1,juf,kf-l) + vf(if+m-1,juf,kf-l-1) +
     &             vf(if+m  ,juf,kf+l) + vf(if+m  ,juf,kf+l-1) +
     &             vf(if+m-1,juf,kf+l) + vf(if+m-1,juf,kf+l-1)) +
c                                             
     &            (wf(if+m  ,juf,kf-l) - wf(if+m  ,juf,kf-l-1) +
     &             wf(if+m-1,juf,kf-l) - wf(if+m-1,juf,kf-l-1) +
     &             wf(if+m  ,juf,kf+l) - wf(if+m  ,juf,kf+l-1) +
     &             wf(if+m-1,juf,kf+l) - wf(if+m-1,juf,kf+l-1)))
                        end do
                     end do
                  end do
               end if
               if (ga(ii,ji,0) - ga(1-ii,ji,0) .eq. 1) then
                  iuf = if + ii - 1
                  fac0 = 1.d0 / (jr * kr)
                  do l = 0, kr-1
                     fac1 = (kr-l) * fac0
                     if (l .eq. 0) fac1 = 0.5d0 * fac1
                     do n = jdir, jdir*(jr-1), jdir
                        fac = (jr-abs(n)) * fac1
                        do kf = kr*cregl2, kr*cregh2, kr
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                       (idir *
     &            (uf(iuf,jf+n,kf-l  ) + uf(iuf,jf+n-1,kf-l  ) +
     &             uf(iuf,jf+n,kf-l-1) + uf(iuf,jf+n-1,kf-l-1) +
     &             uf(iuf,jf+n,kf+l  ) + uf(iuf,jf+n-1,kf+l  ) +
     &             uf(iuf,jf+n,kf+l-1) + uf(iuf,jf+n-1,kf+l-1)) +
c
     &            (vf(iuf,jf+n,kf-l  ) - vf(iuf,jf+n-1,kf-l  ) +
     &             vf(iuf,jf+n,kf-l-1) - vf(iuf,jf+n-1,kf-l-1) +
     &             vf(iuf,jf+n,kf+l  ) - vf(iuf,jf+n-1,kf+l  ) +
     &             vf(iuf,jf+n,kf+l-1) - vf(iuf,jf+n-1,kf+l-1)) +
c                                             
     &            (wf(iuf,jf+n,kf-l  ) - wf(iuf,jf+n  ,kf-l-1) +
     &             wf(iuf,jf+n-1,kf-l) - wf(iuf,jf+n-1,kf-l-1) +
     &             wf(iuf,jf+n,kf+l  ) - wf(iuf,jf+n  ,kf+l-1) +
     &             wf(iuf,jf+n-1,kf+l) - wf(iuf,jf+n-1,kf+l-1)))
                        end do
                     end do
                  end do
               end if
            end do
         end do
c weighting
         do kf = kr*cregl2, kr*cregh2, kr
            src(if,jf,kf) = .25d0 * src(if,jf,kf)
         end do
      end if
      end
c-----------------------------------------------------------------------
c Note---only generates values at coarse points along edge of fine grid
      subroutine hgcdiv_terrain(
     & src, srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & uc, vc, wc,
     &      cl0,ch0,cl1,ch1,cl2,ch2,
     & uf, vf, wf,
     &      fl0,fh0,fl1,fh1,fl2,fh2,
     &      cregl0,cregh0,cregl1,cregh1,cregl2,cregh2,
     & hx, hy, hz, ir, jr, kr, ga)
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer cl0,ch0,cl1,ch1,cl2,ch2
      integer fl0,fh0,fl1,fh1,fl2,fh2
      integer cregl0,cregh0,cregl1,cregh1,cregl2,cregh2
      double precision hx, hy, hz
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2)
      double precision uc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision vc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision wc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision uf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision vf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision wf(fl0:fh0,fl1:fh1,fl2:fh2)
      integer ir, jr, kr, ga(0:1,0:1,0:1)
      double precision sum, fac0, fac1, fac
      integer ic, jc, kc, if, jf, kf, iuc, iuf, juc, juf, kuc, kuf
      integer ii, ji, ki, idir, jdir, kdir, l, m, n
      ic = cregl0
      jc = cregl1
      kc = cregl2
      if = ic * ir
      jf = jc * jr
      kf = kc * kr
      sum = 0.d0
c octants
      do ki = 0, 1
         kdir = 2 * ki - 1
         kuf = kf + ki - 1
         kuc = kc + ki - 1
         do ji = 0, 1
            jdir = 2 * ji - 1
            juf = jf + ji - 1
            juc = jc + ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,ki) .eq. 1) then
                  iuf = if + ii - 1
                  sum = sum + 
     &              (idir * uf(iuf,juf,kuf) +
     &               jdir * vf(iuf,juf,kuf) +
     &               kdir * wf(iuf,juf,kuf))
               else
                  iuc = ic + ii - 1
                  sum = sum + 
     &              (idir * uc(iuc,juc,kuc) +
     &               jdir * vc(iuc,juc,kuc) +
     &               kdir * wc(iuc,juc,kuc))
               end if
            end do
         end do
      end do
c faces
      do ki = 0, 1
         kdir = 2 * ki - 1
         kuf = kf + ki - 1
         do ji = 0, 1
            jdir = 2 * ji - 1
            juf = jf + ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               iuf = if + ii - 1
               if (ga(ii,ji,ki) - ga(ii,ji,1-ki) .eq. 1) then
                  fac0 = 1.d0 / (ir * jr)
                  do n = jdir, jdir*(jr-1), jdir
                     fac1 = (jr-abs(n)) * fac0
                     do m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        sum = sum + fac * (
     &               (uf(if+m,jf+n  ,kuf) - uf(if+m-1,jf+n  ,kuf) +
     &                uf(if+m,jf+n-1,kuf) - uf(if+m-1,jf+n-1,kuf)) +
c
     &               (vf(if+m  ,jf+n,kuf) - vf(if+m  ,jf+n-1,kuf) +
     &                vf(if+m-1,jf+n,kuf) - vf(if+m-1,jf+n-1,kuf)) +
     &                                     kdir *
     &               (wf(if+m  ,jf+n,kuf) + wf(if+m  ,jf+n-1,kuf) +
     &                wf(if+m-1,jf+n,kuf) + wf(if+m-1,jf+n-1,kuf)))
                     end do
                  end do
               end if
               if (ga(ii,ji,ki) - ga(ii,1-ji,ki) .eq. 1) then
                  fac0 = 1.d0 / (ir * kr)
                  do l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        sum = sum + fac * (
     &                (uf(if+m  ,juf,kf+l  ) - uf(if+m-1,juf,kf+l  ) +
     &                 uf(if+m  ,juf,kf+l-1) - uf(if+m-1,juf,kf+l-1)) +
     &                                     jdir *
     &                (vf(if+m  ,juf,kf+l  ) + vf(if+m  ,juf,kf+l-1) +
     &                 vf(if+m-1,juf,kf+l  ) + vf(if+m-1,juf,kf+l-1)) +
c
     &                (wf(if+m  ,juf,kf+l  ) - wf(if+m  ,juf,kf+l-1) +
     &                 wf(if+m-1,juf,kf+l  ) - wf(if+m-1,juf,kf+l-1)))
                     end do
                  end do
               end if
               if (ga(ii,ji,ki) - ga(1-ii,ji,ki) .eq. 1) then
                  fac0 = 1.d0 / (jr * kr)
                  do l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do n = jdir, jdir*(jr-1), jdir
                        fac = (jr-abs(n)) * fac1
                        sum = sum + fac * (idir *
     &                (uf(iuf,jf+n  ,kf+l  ) + uf(iuf,jf+n  ,kf+l-1) +
     &                 uf(iuf,jf+n-1,kf+l  ) + uf(iuf,jf+n-1,kf+l-1)) +
c
     &                (vf(iuf,jf+n  ,kf+l  ) - vf(iuf,jf+n-1,kf+l  ) +
     &                 vf(iuf,jf+n  ,kf+l-1) - vf(iuf,jf+n-1,kf+l-1)) +
c
     &                (wf(iuf,jf+n  ,kf+l  ) - wf(iuf,jf+n  ,kf+l-1) +
     &                 wf(iuf,jf+n-1,kf+l  ) - wf(iuf,jf+n-1,kf+l-1)))
                     end do
                  end do
               end if
            end do
         end do
      end do
c edges
      do ki = 0, 1
         kdir = 2 * ki - 1
         kuf = kf + ki - 1
         do ji = 0, 1
            jdir = 2 * ji - 1
            juf = jf + ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               iuf = if + ii - 1
               if (ga(ii,ji,ki) -
     &             min(ga(ii,ji,1-ki), ga(ii,1-ji,ki), ga(ii,1-ji,1-ki))
     &             .eq. 1) then
                  fac1 = 1.d0 / ir
                  do m = idir, idir*(ir-1), idir
                     fac = (ir-abs(m)) * fac1
                     sum = sum + fac * (
     &                 (uf(if+m,juf,kuf) - uf(if+m-1,juf,kuf)) +
     &                                  jdir *
     &                 (vf(if+m,juf,kuf) + vf(if+m-1,juf,kuf)) +
     &                                  kdir *
     &                 (wf(if+m,juf,kuf) + wf(if+m-1,juf,kuf)))
                  end do
               end if
               if (ga(ii,ji,ki) -
     &             min(ga(ii,ji,1-ki), ga(1-ii,ji,ki), ga(1-ii,ji,1-ki))
     &             .eq. 1) then
                  fac1 = 1.d0 / jr
                  do n = jdir, jdir*(jr-1), jdir
                     fac = (jr-abs(n)) * fac1
                     sum = sum + fac * (idir *
     &                 (uf(iuf,jf+n,kuf) + uf(iuf,jf+n-1,kuf)) +
c
     &                 (vf(iuf,jf+n,kuf) - vf(iuf,jf+n-1,kuf)) +
     &                                  kdir *
     &                 (wf(iuf,jf+n,kuf) + wf(iuf,jf+n-1,kuf)))
                  end do
               end if
               if (ga(ii,ji,ki) -
     &             min(ga(ii,1-ji,ki), ga(1-ii,ji,ki), ga(1-ii,1-ji,ki))
     &             .eq. 1) then
                  fac1 = 1.d0 / kr
                  do l = kdir, kdir*(kr-1), kdir
                     fac = (kr-abs(l)) * fac1
                     sum = sum + fac * (idir *
     &                 (uf(iuf,juf,kf+l) + uf(iuf,juf,kf+l-1)) +
     &                                  jdir *
     &                 (vf(iuf,juf,kf+l) + vf(iuf,juf,kf+l-1)) +
c
     &                 (wf(iuf,juf,kf+l) - wf(iuf,juf,kf+l-1)))
                  end do
               end if
            end do
         end do
      end do
c weighting
      src(if,jf,kf) = .25d0 * sum 
      end
c-----------------------------------------------------------------------
      subroutine hggrad(
     & gpx, gpy, gpz,
     &       gpl0,gph0,gpl1,gph1,gpl2,gph2,
     & dest, destl0,desth0,destl1,desth1,destl2,desth2,
     &       fregl0,fregh0,fregl1,fregh1,fregl2,fregh2,
     & hx, hy, hz)
      integer gpl0,gph0,gpl1,gph1,gpl2,gph2
      integer destl0,desth0,destl1,desth1,destl2,desth2
      integer fregl0,fregh0,fregl1,fregh1,fregl2,fregh2
      double precision gpx(gpl0:gph0,gpl1:gph1,gpl2:gph2)
      double precision gpy(gpl0:gph0,gpl1:gph1,gpl2:gph2)
      double precision gpz(gpl0:gph0,gpl1:gph1,gpl2:gph2)
      double precision dest(destl0:desth0,destl1:desth1,destl2:desth2)
      double precision hx, hy, hz
      double precision hxm1h, hym1h, hzm1h
      integer i, j, k
      hxm1h = 0.25d0 / hx
      hym1h = 0.25d0 / hy
      hzm1h = 0.25d0 / hz
      do k = fregl2, fregh2
         do j = fregl1, fregh1
            do i = fregl0, fregh0
               gpx(i,j,k) = hxm1h * (dest(i+1,j,k) + dest(i+1,j+1,k) +
     &                           dest(i+1,j,k+1) + dest(i+1,j+1,k+1) -
     &                           dest(i,j,k) - dest(i,j+1,k) -
     &                           dest(i,j,k+1) - dest(i,j+1,k+1))
               gpy(i,j,k) = hym1h * (dest(i,j+1,k) + dest(i+1,j+1,k) +
     &                           dest(i,j+1,k+1) + dest(i+1,j+1,k+1) -
     &                           dest(i,j,k) - dest(i+1,j,k) -
     &                           dest(i,j,k+1) - dest(i+1,j,k+1))
               gpz(i,j,k) = hzm1h * (dest(i,j,k+1) + dest(i+1,j,k+1) +
     &                           dest(i,j+1,k+1) + dest(i+1,j+1,k+1) -
     &                           dest(i,j,k) - dest(i+1,j,k) -
     &                           dest(i,j+1,k) - dest(i+1,j+1,k))
            end do
         end do
      end do
      end
c-----------------------------------------------------------------------
      subroutine hgdiv(
     & src, srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & uf, vf, wf,
     &      fl0,fh0,fl1,fh1,fl2,fh2,
     &      fregl0,fregh0,fregl1,fregh1,fregl2,fregh2,
     & hx, hy, hz)
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer fl0,fh0,fl1,fh1,fl2,fh2
      integer fregl0,fregh0,fregl1,fregh1,fregl2,fregh2
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2)
      double precision uf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision vf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision wf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision hx, hy, hz
      double precision hxm1, hym1, hzm1, fac
      integer i, j, k
      hxm1 = 1.d0 / hx
      hym1 = 1.d0 / hy
      hzm1 = 1.d0 / hz
      fac = 0.25d0
      do k = fregl2, fregh2
         do j = fregl1, fregh1
            do i = fregl0, fregh0
               src(i,j,k) = fac *
     &              (hxm1 * (uf(i  ,j-1,k-1) - uf(i-1,j-1,k-1) +
     &                       uf(i  ,j-1,k  ) - uf(i-1,j-1,k) +
     &                       uf(i  ,j  ,k-1) - uf(i-1,j  ,k-1) +
     &                       uf(i  ,j  ,k  ) - uf(i-1,j  ,k  )) +
     &               hym1 * (vf(i-1,j  ,k-1) - vf(i-1,j-1,k-1) +
     &                       vf(i-1,j  ,k  ) - vf(i-1,j-1,k  ) +
     &                       vf(i  ,j  ,k-1) - vf(i  ,j-1,k-1) +
     &                       vf(i  ,j  ,k  ) - vf(i  ,j-1,k  )) +
     &               hzm1 * (wf(i-1,j-1,k  ) - wf(i-1,j-1,k-1) +
     &                       wf(i-1,j  ,k  ) - wf(i-1,j  ,k-1) +
     &                       wf(i  ,j-1,k  ) - wf(i  ,j-1,k-1) +
     &                       wf(i  ,j  ,k  ) - wf(i  ,j  ,k-1)))
            end do
         end do
      end do
      end
c-----------------------------------------------------------------------
c Note---only generates values at coarse points along edge of fine grid
      subroutine hgfdiv(
     & src, srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & uc, vc, wc,
     &      cl0,ch0,cl1,ch1,cl2,ch2,
     & uf, vf, wf,
     &      fl0,fh0,fl1,fh1,fl2,fh2,
     & cregl0,cregh0,cregl1,cregh1,cregl2,cregh2,
     & hx, hy, hz, ir, jr, kr, idim, idir, i1, i2)
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer cl0,ch0,cl1,ch1,cl2,ch2
      integer fl0,fh0,fl1,fh1,fl2,fh2
      integer cregl0,cregh0,cregl1,cregh1,cregl2,cregh2
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2)
      double precision uc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision vc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision wc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision uf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision vf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision wf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision hx, hy, hz
      integer ir, jr, kr, idim, idir, i1, i2
      double precision hxm1, hym1, hzm1, fac0, fac1, fac, tmp
      integer i, j, k, iuc, iuf, juc, juf, kuc, kuf, l, m, n
      if (idim .eq. 0) then
         i = cregl0
         if (idir .eq. 1) then
            iuc = i - 1
            iuf = i * ir
         else
            iuc = i
            iuf = i * ir - 1
         end if
         fac0 = 0.5d0 * ir / (ir + 1.d0)
         hxm1 = 1.d0 / (ir * hx)
         hym1 = 1.d0 / (jr * hy)
         hzm1 = 1.d0 / (kr * hz)
         do k = cregl2, cregh2
            do j = cregl1, cregh1
               src(i*ir,j*jr,k*kr) = fac0 *
     &           (hym1 * (vc(iuc,j  ,k  ) - vc(iuc,j-1,k  ) +
     &                    vc(iuc,j  ,k-1) - vc(iuc,j-1,k-1)) +
     &            hzm1 * (wc(iuc,j  ,k  ) - wc(iuc,j  ,k-1) +
     &                    wc(iuc,j-1,k  ) - wc(iuc,j-1,k-1)) -
     &            hxm1 * idir * (uc(iuc,j  ,k) + uc(iuc,j  ,k-1) +
     &                           uc(iuc,j-1,k) + uc(iuc,j-1,k-1)))
            end do
         end do
         fac0 = fac0 / (ir * jr * kr * jr * kr)
         hxm1 = ir * hxm1
         hym1 = jr * hym1
         hzm1 = kr * hzm1
         i = i * ir
         do l = 0, kr-1
            fac1 = (kr-l) * fac0
            if (l .eq. 0) fac1 = 0.5d0 * fac1
            do n = 0, jr-1
               fac = (jr-n) * fac1
               if (n .eq. 0) fac = 0.5d0 * fac
               do k = kr*cregl2, kr*cregh2, kr
                  do j = jr*cregl1, jr*cregh1, jr
                     tmp = hxm1 * idir *
     &                 (uf(iuf,j-n  ,k-l) + uf(iuf,j-n  ,k-l-1) +
     &                  uf(iuf,j-n-1,k-l) + uf(iuf,j-n-1,k-l-1) +
     &                  uf(iuf,j-n  ,k+l) + uf(iuf,j-n  ,k+l-1) +
     &                  uf(iuf,j-n-1,k+l) + uf(iuf,j-n-1,k+l-1) +
     &                  uf(iuf,j+n  ,k-l) + uf(iuf,j+n  ,k-l-1) +
     &                  uf(iuf,j+n-1,k-l) + uf(iuf,j+n-1,k-l-1) +
     &                  uf(iuf,j+n  ,k+l) + uf(iuf,j+n  ,k+l-1) +
     &                  uf(iuf,j+n-1,k+l) + uf(iuf,j+n-1,k+l-1))
                     tmp = tmp + hym1 *
     &                 (vf(iuf,j-n,k-l  ) - vf(iuf,j-n-1,k-l  ) +
     &                  vf(iuf,j-n,k-l-1) - vf(iuf,j-n-1,k-l-1) +
     &                  vf(iuf,j-n,k+l  ) - vf(iuf,j-n-1,k+l  ) +
     &                  vf(iuf,j-n,k+l-1) - vf(iuf,j-n-1,k+l-1) +
     &                  vf(iuf,j+n,k-l  ) - vf(iuf,j+n-1,k-l  ) +
     &                  vf(iuf,j+n,k-l-1) - vf(iuf,j+n-1,k-l-1) +
     &                  vf(iuf,j+n,k+l  ) - vf(iuf,j+n-1,k+l  ) +
     &                  vf(iuf,j+n,k+l-1) - vf(iuf,j+n-1,k+l-1))
                     src(i,j,k) = src(i,j,k) + fac * (tmp + hzm1 *
     &                 (wf(iuf,j-n,k-l  ) - wf(iuf,j-n  ,k-l-1) +
     &                  wf(iuf,j-n-1,k-l) - wf(iuf,j-n-1,k-l-1) +
     &                  wf(iuf,j-n,k+l  ) - wf(iuf,j-n  ,k+l-1) +
     &                  wf(iuf,j-n-1,k+l) - wf(iuf,j-n-1,k+l-1) +
     &                  wf(iuf,j+n,k-l  ) - wf(iuf,j+n  ,k-l-1) +
     &                  wf(iuf,j+n-1,k-l) - wf(iuf,j+n-1,k-l-1) +
     &                  wf(iuf,j+n,k+l  ) - wf(iuf,j+n  ,k+l-1) +
     &                  wf(iuf,j+n-1,k+l) - wf(iuf,j+n-1,k+l-1)))
                  end do
               end do
            end do
         end do
      else if (idim .eq. 1) then
         j = cregl1
         if (idir .eq. 1) then
            juc = j - 1
            juf = j * jr
         else
            juc = j
            juf = j * jr - 1
         end if
         fac0 = 0.5d0 * jr / (jr + 1.d0)
         hxm1 = 1.d0 / (ir * hx)
         hym1 = 1.d0 / (jr * hy)
         hzm1 = 1.d0 / (kr * hz)
         do k = cregl2, cregh2
            do i = cregl0, cregh0
               src(i*ir,j*jr,k*kr) = fac0 *
     &           (hxm1 * (uc(i,juc,k  ) - uc(i-1,juc,k  ) +
     &                    uc(i,juc,k-1) - uc(i-1,juc,k-1)) -
     &            hym1 * idir * (vc(i  ,juc,k) + vc(i  ,juc,k-1) +
     &                           vc(i-1,juc,k) + vc(i-1,juc,k-1)) +
     &            hzm1 * (wc(i  ,juc,k) - wc(i  ,juc,k-1) +
     &                    wc(i-1,juc,k) - wc(i-1,juc,k-1)))
            end do
         end do
         fac0 = fac0 / (ir * jr * kr * ir * kr)
         hxm1 = ir * hxm1
         hym1 = jr * hym1
         hzm1 = kr * hzm1
         j = j * jr
         do l = 0, kr-1
            fac1 = (kr-l) * fac0
            if (l .eq. 0) fac1 = 0.5d0 * fac1
            do m = 0, ir-1
               fac = (ir-m) * fac1
               if (m .eq. 0) fac = 0.5d0 * fac
               do k = kr*cregl2, kr*cregh2, kr
                  do i = ir*cregl0, ir*cregh0, ir
                     tmp = hxm1 *
     &                 (uf(i-m,juf,k-l  ) - uf(i-m-1,juf,k-l  ) +
     &                  uf(i-m,juf,k-l-1) - uf(i-m-1,juf,k-l-1) +
     &                  uf(i-m,juf,k+l  ) - uf(i-m-1,juf,k+l  ) +
     &                  uf(i-m,juf,k+l-1) - uf(i-m-1,juf,k+l-1) +
     &                  uf(i+m,juf,k-l  ) - uf(i+m-1,juf,k-l  ) +
     &                  uf(i+m,juf,k-l-1) - uf(i+m-1,juf,k-l-1) +
     &                  uf(i+m,juf,k+l  ) - uf(i+m-1,juf,k+l  ) +
     &                  uf(i+m,juf,k+l-1) - uf(i+m-1,juf,k+l-1))
                     tmp = tmp + hym1 * idir *
     &                 (vf(i-m  ,juf,k-l) + vf(i-m  ,juf,k-l-1) +
     &                  vf(i-m-1,juf,k-l) + vf(i-m-1,juf,k-l-1) +
     &                  vf(i-m  ,juf,k+l) + vf(i-m  ,juf,k+l-1) +
     &                  vf(i-m-1,juf,k+l) + vf(i-m-1,juf,k+l-1) +
     &                  vf(i+m  ,juf,k-l) + vf(i+m  ,juf,k-l-1) +
     &                  vf(i+m-1,juf,k-l) + vf(i+m-1,juf,k-l-1) +
     &                  vf(i+m  ,juf,k+l) + vf(i+m  ,juf,k+l-1) +
     &                  vf(i+m-1,juf,k+l) + vf(i+m-1,juf,k+l-1))
                     src(i,j,k) = src(i,j,k) + fac * (tmp + hzm1 *
     &                 (wf(i-m  ,juf,k-l) - wf(i-m  ,juf,k-l-1) +
     &                  wf(i-m-1,juf,k-l) - wf(i-m-1,juf,k-l-1) +
     &                  wf(i-m  ,juf,k+l) - wf(i-m  ,juf,k+l-1) +
     &                  wf(i-m-1,juf,k+l) - wf(i-m-1,juf,k+l-1) +
     &                  wf(i+m  ,juf,k-l) - wf(i+m  ,juf,k-l-1) +
     &                  wf(i+m-1,juf,k-l) - wf(i+m-1,juf,k-l-1) +
     &                  wf(i+m  ,juf,k+l) - wf(i+m  ,juf,k+l-1) +
     &                  wf(i+m-1,juf,k+l) - wf(i+m-1,juf,k+l-1)))
                  end do
               end do
            end do
         end do
      else
         k = cregl2
         if (idir .eq. 1) then
            kuc = k - 1
            kuf = k * kr
         else
            kuc = k
            kuf = k * kr - 1
         end if
         fac0 = 0.5d0 * kr / (kr + 1.d0)
         hxm1 = 1.d0 / (ir * hx)
         hym1 = 1.d0 / (jr * hy)
         hzm1 = 1.d0 / (kr * hz)
         do j = cregl1, cregh1
            do i = cregl0, cregh0
               src(i*ir,j*jr,k*kr) = fac0 *
     &           (hxm1 * (uc(i  ,j  ,kuc) - uc(i-1,j  ,kuc) +
     &                    uc(i  ,j-1,kuc) - uc(i-1,j-1,kuc)) +
     &            hym1 * (vc(i  ,j  ,kuc) - vc(i  ,j-1,kuc) +
     &                    vc(i-1,j  ,kuc) - vc(i-1,j-1,kuc)) -
     &            hzm1 * idir * (wc(i  ,j,kuc) + wc(i  ,j-1,kuc) +
     &                           wc(i-1,j,kuc) + wc(i-1,j-1,kuc)))
            end do
         end do
         fac0 = fac0 / (ir * jr * kr * ir * jr)
         hxm1 = ir * hxm1
         hym1 = jr * hym1
         hzm1 = kr * hzm1
         k = k * kr
         do n = 0, jr-1
            fac1 = (jr-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do m = 0, ir-1
               fac = (ir-m) * fac1
               if (m .eq. 0) fac = 0.5d0 * fac
               do j = jr*cregl1, jr*cregh1, jr
                  do i = ir*cregl0, ir*cregh0, ir
                     tmp = hxm1 *
     &                 (uf(i-m,j-n  ,kuf) - uf(i-m-1,j-n  ,kuf) +
     &                  uf(i-m,j-n-1,kuf) - uf(i-m-1,j-n-1,kuf) +
     &                  uf(i-m,j+n  ,kuf) - uf(i-m-1,j+n  ,kuf) +
     &                  uf(i-m,j+n-1,kuf) - uf(i-m-1,j+n-1,kuf) +
     &                  uf(i+m,j-n  ,kuf) - uf(i+m-1,j-n  ,kuf) +
     &                  uf(i+m,j-n-1,kuf) - uf(i+m-1,j-n-1,kuf) +
     &                  uf(i+m,j+n  ,kuf) - uf(i+m-1,j+n  ,kuf) +
     &                  uf(i+m,j+n-1,kuf) - uf(i+m-1,j+n-1,kuf))
                     tmp = tmp + hym1 *
     &                 (vf(i-m  ,j-n,kuf) - vf(i-m  ,j-n-1,kuf) +
     &                  vf(i-m-1,j-n,kuf) - vf(i-m-1,j-n-1,kuf) +
     &                  vf(i-m  ,j+n,kuf) - vf(i-m  ,j+n-1,kuf) +
     &                  vf(i-m-1,j+n,kuf) - vf(i-m-1,j+n-1,kuf) +
     &                  vf(i+m  ,j-n,kuf) - vf(i+m  ,j-n-1,kuf) +
     &                  vf(i+m-1,j-n,kuf) - vf(i+m-1,j-n-1,kuf) +
     &                  vf(i+m  ,j+n,kuf) - vf(i+m  ,j+n-1,kuf) +
     &                  vf(i+m-1,j+n,kuf) - vf(i+m-1,j+n-1,kuf))
                    src(i,j,k) = src(i,j,k) + fac * (tmp + hzm1 * idir *
     &                 (wf(i-m  ,j-n,kuf) + wf(i-m  ,j-n-1,kuf) +
     &                  wf(i-m-1,j-n,kuf) + wf(i-m-1,j-n-1,kuf) +
     &                  wf(i-m  ,j+n,kuf) + wf(i-m  ,j+n-1,kuf) +
     &                  wf(i-m-1,j+n,kuf) + wf(i-m-1,j+n-1,kuf) +
     &                  wf(i+m  ,j-n,kuf) + wf(i+m  ,j-n-1,kuf) +
     &                  wf(i+m-1,j-n,kuf) + wf(i+m-1,j-n-1,kuf) +
     &                  wf(i+m  ,j+n,kuf) + wf(i+m  ,j+n-1,kuf) +
     &                  wf(i+m-1,j+n,kuf) + wf(i+m-1,j+n-1,kuf)))
                  end do
               end do
            end do
         end do
      end if
      end
c-----------------------------------------------------------------------
c Note---only generates values at coarse points along edge of fine grid
      subroutine hgediv(
     & src, srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & uc, vc, wc,
     &      cl0,ch0,cl1,ch1,cl2,ch2,
     & uf, vf, wf,
     &      fl0,fh0,fl1,fh1,fl2,fh2,
     &      cregl0,cregh0,cregl1,cregh1,cregl2,cregh2,
     & hx, hy, hz, ir, jr, kr, ivect, ga)
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer cl0,ch0,cl1,ch1,cl2,ch2
      integer fl0,fh0,fl1,fh1,fl2,fh2
      integer cregl0,cregh0,cregl1,cregh1,cregl2,cregh2
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2)
      double precision uc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision vc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision wc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision uf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision vf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision wf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision hx, hy, hz
      integer ir, jr, kr, ivect(0:2), ga(0:1,0:1,0:1)
      double precision r3, hxm1, hym1, hzm1, hxm1c, hym1c, hzm1c
      double precision center, cfac, ffac, fac0, fac1, fac
      integer ic, jc, kc, if, jf, kf, iuc, iuf, juc, juf, kuc, kuf
      integer ii, ji, ki, idir, jdir, kdir, l, m, n
      r3 = ir * jr * kr
      hxm1c = 1.d0 / (ir * hx)
      hym1c = 1.d0 / (jr * hy)
      hzm1c = 1.d0 / (kr * hz)
      hxm1 = ir * hxm1c
      hym1 = jr * hym1c
      hzm1 = kr * hzm1c
      ic = cregl0
      jc = cregl1
      kc = cregl2
      if = ic * ir
      jf = jc * jr
      kf = kc * kr
      center = 0.d0
      if (ivect(0) .eq. 0) then
         do if = ir*cregl0, ir*cregh0, ir
            src(if,jf,kf) = 0.d0
         end do
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac1 = 1.d0 / ir
         ffac = ir
         cfac = r3
         do ki = 0, 1
            kdir = 2 * ki - 1
            do ji = 0, 1
               jdir = 2 * ji - 1
               if (ga(0,ji,ki) .eq. 1) then
                  kuf = kf + ki - 1
                  juf = jf + ji - 1
                  center = center + ffac
                  do m = 0, ir-1
                     fac = (ir-m) * fac1
                     if (m .eq. 0) fac = 0.5d0 * fac
                     do if = ir*cregl0, ir*cregh0, ir
                        src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                 (hxm1 *
     &                    (uf(if-m,juf,kuf) - uf(if-m-1,juf,kuf) +
     &                     uf(if+m,juf,kuf) - uf(if+m-1,juf,kuf)) +
     &                                  hym1 * jdir *
     &                    (vf(if-m,juf,kuf) + vf(if-m-1,juf,kuf) +
     &                     vf(if+m,juf,kuf) + vf(if+m-1,juf,kuf)) +
     &                                  hzm1 * kdir *
     &                    (wf(if-m,juf,kuf) + wf(if-m-1,juf,kuf) +
     &                     wf(if+m,juf,kuf) + wf(if+m-1,juf,kuf)))
                     end do
                  end do
               else
                  kuc = kc + ki - 1
                  juc = jc + ji - 1
                  center = center + cfac
                  do ic = cregl0, cregh0
                     if = ic * ir
                     src(if,jf,kf) = src(if,jf,kf) + r3 *
     &            (hxm1c *        (uc(ic,juc,kuc) - uc(ic-1,juc,kuc)) +
     &             hym1c * jdir * (vc(ic,juc,kuc) + vc(ic-1,juc,kuc)) +
     &             hzm1c * kdir * (wc(ic,juc,kuc) + wc(ic-1,juc,kuc)))
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
                  kuf = kf + ki - 1
                  fac0 = 1.d0 / (ir * jr)
                  ffac = ir * (jr - 1)
                  center = center + ffac
                  do n = jdir, jdir*(jr-1), jdir
                     fac1 = (jr-abs(n)) * fac0
                     do m = 0, ir-1
                        fac = (ir-m) * fac1
                        if (m .eq. 0) fac = 0.5d0 * fac
                        do if = ir*cregl0, ir*cregh0, ir
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                       (hxm1 *
     &                (uf(if-m,jf+n,kuf) - uf(if-m-1,jf+n,kuf) +
     &                 uf(if-m,jf+n-1,kuf) - uf(if-m-1,jf+n-1,kuf) +
     &                 uf(if+m,jf+n,kuf) - uf(if+m-1,jf+n,kuf) +
     &                 uf(if+m,jf+n-1,kuf) - uf(if+m-1,jf+n-1,kuf)) +
     &                                        hym1 *
     &                (vf(if-m,jf+n,kuf) - vf(if-m,jf+n-1,kuf) +
     &                 vf(if-m-1,jf+n,kuf) - vf(if-m-1,jf+n-1,kuf) +
     &                 vf(if+m,jf+n,kuf) - vf(if+m,jf+n-1,kuf) +
     &                 vf(if+m-1,jf+n,kuf) - vf(if+m-1,jf+n-1,kuf)) +
     &                                        hzm1 * kdir *
     &                (wf(if-m,jf+n,kuf) + wf(if-m,jf+n-1,kuf) +
     &                 wf(if-m-1,jf+n,kuf) + wf(if-m-1,jf+n-1,kuf) +
     &                 wf(if+m,jf+n,kuf) + wf(if+m,jf+n-1,kuf) +
     &                 wf(if+m-1,jf+n,kuf) + wf(if+m-1,jf+n-1,kuf)))
                        end do
                     end do
                  end do
               end if
               if (ga(0,ji,ki) - ga(0,1-ji,ki) .eq. 1) then
                  juf = jf + ji - 1
                  fac0 = 1.d0 / (ir * kr)
                  ffac = ir * (kr - 1)
                  center = center + ffac
                  do l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do m = 0, ir-1
                        fac = (ir-m) * fac1
                        if (m .eq. 0) fac = 0.5d0 * fac
                        do if = ir*cregl0, ir*cregh0, ir
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                       (hxm1 *
     &                (uf(if-m,juf,kf+l) - uf(if-m-1,juf,kf+l) +
     &                 uf(if-m,juf,kf+l-1) - uf(if-m-1,juf,kf+l-1) +
     &                 uf(if+m,juf,kf+l) - uf(if+m-1,juf,kf+l) +
     &                 uf(if+m,juf,kf+l-1) - uf(if+m-1,juf,kf+l-1)) +
     &                                        hym1 * jdir *
     &                (vf(if-m,juf,kf+l) + vf(if-m,juf,kf+l-1) +
     &                 vf(if-m-1,juf,kf+l) + vf(if-m-1,juf,kf+l-1) +
     &                 vf(if+m,juf,kf+l) + vf(if+m,juf,kf+l-1) +
     &                 vf(if+m-1,juf,kf+l) + vf(if+m-1,juf,kf+l-1)) +
     &                                        hzm1 *
     &                (wf(if-m,juf,kf+l) - wf(if-m,juf,kf+l-1) +
     &                 wf(if-m-1,juf,kf+l) - wf(if-m-1,juf,kf+l-1) +
     &                 wf(if+m,juf,kf+l) - wf(if+m,juf,kf+l-1) +
     &                 wf(if+m-1,juf,kf+l) - wf(if+m-1,juf,kf+l-1)))
                        end do
                     end do
                  end do
               end if
            end do
         end do
c weighting 
         do if = ir*cregl0, ir*cregh0, ir
            src(if,jf,kf) = src(if,jf,kf) / center
         end do
      else if (ivect(1) .eq. 0) then
         do jf = jr*cregl1, jr*cregh1, jr
            src(if,jf,kf) = 0.d0
         end do
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac1 = 1.d0 / jr
         ffac = jr
         cfac = r3
         do ki = 0, 1
            kdir = 2 * ki - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,0,ki) .eq. 1) then
                  kuf = kf + ki - 1
                  iuf = if + ii - 1
                  center = center + ffac
                  do n = 0, jr-1
                     fac = (jr-n) * fac1
                     if (n .eq. 0) fac = 0.5d0 * fac
                     do jf = jr*cregl1, jr*cregh1, jr
                        src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                 (hxm1 * idir *
     &                    (uf(iuf,jf-n,kuf) + uf(iuf,jf-n-1,kuf) +
     &                     uf(iuf,jf+n,kuf) + uf(iuf,jf+n-1,kuf)) +
     &                                  hym1 *
     &                    (vf(iuf,jf-n,kuf) - vf(iuf,jf-n-1,kuf) +
     &                     vf(iuf,jf+n,kuf) - vf(iuf,jf+n-1,kuf)) +
     &                                  hzm1 * kdir *
     &                    (wf(iuf,jf-n,kuf) + wf(iuf,jf-n-1,kuf) +
     &                     wf(iuf,jf+n,kuf) + wf(iuf,jf+n-1,kuf)))
                     end do
                  end do
               else
                  kuc = kc + ki - 1
                  iuc = ic + ii - 1
                  center = center + cfac
                  do jc = cregl1, cregh1
                     jf = jc * jr
                     src(if,jf,kf) = src(if,jf,kf) + r3 *
     &            (hxm1c * idir * (uc(iuc,jc,kuc) + uc(iuc,jc-1,kuc)) +
     &             hym1c *        (vc(iuc,jc,kuc) - vc(iuc,jc-1,kuc)) +
     &             hzm1c * kdir * (wc(iuc,jc,kuc) + wc(iuc,jc-1,kuc)))
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
                  kuf = kf + ki - 1
                  fac0 = 1.d0 / (ir * jr)
                  ffac = jr * (ir - 1)
                  center = center + ffac
                  do n = 0, jr-1
                     fac1 = (jr-n) * fac0
                     if (n .eq. 0) fac1 = 0.5d0 * fac1
                     do m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        do jf = jr*cregl1, jr*cregh1, jr
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                       (hxm1 *
     &                (uf(if+m,jf-n,kuf) - uf(if+m-1,jf-n,kuf) +
     &                 uf(if+m,jf-n-1,kuf) - uf(if+m-1,jf-n-1,kuf) +
     &                 uf(if+m,jf+n,kuf) - uf(if+m-1,jf+n,kuf) +
     &                 uf(if+m,jf+n-1,kuf) - uf(if+m-1,jf+n-1,kuf)) +
     &                                        hym1 *
     &                (vf(if+m,jf-n,kuf) - vf(if+m,jf-n-1,kuf) +
     &                 vf(if+m-1,jf-n,kuf) - vf(if+m-1,jf-n-1,kuf) +
     &                 vf(if+m,jf+n,kuf) - vf(if+m,jf+n-1,kuf) +
     &                 vf(if+m-1,jf+n,kuf) - vf(if+m-1,jf+n-1,kuf)) +
     &                                        hzm1 * kdir *
     &                (wf(if+m,jf-n,kuf) + wf(if+m,jf-n-1,kuf) +
     &                 wf(if+m-1,jf-n,kuf) + wf(if+m-1,jf-n-1,kuf) +
     &                 wf(if+m,jf+n,kuf) + wf(if+m,jf+n-1,kuf) +
     &                 wf(if+m-1,jf+n,kuf) + wf(if+m-1,jf+n-1,kuf)))
                        end do
                     end do
                  end do
               end if
               if (ga(ii,0,ki) - ga(1-ii,0,ki) .eq. 1) then
                  iuf = if + ii - 1
                  fac0 = 1.d0 / (jr * kr)
                  ffac = jr * (kr - 1)
                  center = center + ffac
                  do l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do n = 0, jr-1
                        fac = (jr-n) * fac1
                        if (n .eq. 0) fac = 0.5d0 * fac
                        do jf = jr*cregl1, jr*cregh1, jr
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                       (hxm1 * idir *
     &                (uf(iuf,jf-n,kf+l) + uf(iuf,jf-n-1,kf+l) +
     &                 uf(iuf,jf-n,kf+l-1) + uf(iuf,jf-n-1,kf+l-1) +
     &                 uf(iuf,jf+n,kf+l) + uf(iuf,jf+n-1,kf+l) +
     &                 uf(iuf,jf+n,kf+l-1) + uf(iuf,jf+n-1,kf+l-1)) +
     &                                        hym1 *
     &                (vf(iuf,jf-n,kf+l) - vf(iuf,jf-n-1,kf+l) +
     &                 vf(iuf,jf-n,kf+l-1) - vf(iuf,jf-n-1,kf+l-1) +
     &                 vf(iuf,jf+n,kf+l) - vf(iuf,jf+n-1,kf+l) +
     &                 vf(iuf,jf+n,kf+l-1) - vf(iuf,jf+n-1,kf+l-1)) +
     &                                        hzm1 *
     &                (wf(iuf,jf-n,kf+l) - wf(iuf,jf-n,kf+l-1) +
     &                 wf(iuf,jf-n-1,kf+l) - wf(iuf,jf-n-1,kf+l-1) +
     &                 wf(iuf,jf+n,kf+l) - wf(iuf,jf+n,kf+l-1) +
     &                 wf(iuf,jf+n-1,kf+l) - wf(iuf,jf+n-1,kf+l-1)))
                        end do
                     end do
                  end do
               end if
            end do
         end do
c weighting
         do jf = jr*cregl1, jr*cregh1, jr
            src(if,jf,kf) = src(if,jf,kf) / center
         end do
      else
         do kf = kr*cregl2, kr*cregh2, kr
            src(if,jf,kf) = 0.d0
         end do
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac1 = 1.d0 / kr
         ffac = kr
         cfac = r3
         do ji = 0, 1
            jdir = 2 * ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,0) .eq. 1) then
                  juf = jf + ji - 1
                  iuf = if + ii - 1
                  center = center + ffac
                  do l = 0, kr-1
                     fac = (kr-l) * fac1
                     if (l .eq. 0) fac = 0.5d0 * fac
                     do kf = kr*cregl2, kr*cregh2, kr
                        src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                 (hxm1 * idir *
     &                    (uf(iuf,juf,kf-l) + uf(iuf,juf,kf-l-1) +
     &                     uf(iuf,juf,kf+l) + uf(iuf,juf,kf+l-1)) +
     &                                  hym1 * jdir *
     &                    (vf(iuf,juf,kf-l) + vf(iuf,juf,kf-l-1) +
     &                     vf(iuf,juf,kf+l) + vf(iuf,juf,kf+l-1)) +
     &                                  hzm1 *
     &                    (wf(iuf,juf,kf-l) - wf(iuf,juf,kf-l-1) +
     &                     wf(iuf,juf,kf+l) - wf(iuf,juf,kf+l-1)))
                     end do
                  end do
               else
                  juc = jc + ji - 1
                  iuc = ic + ii - 1
                  center = center + cfac
                  do kc = cregl2, cregh2
                     kf = kc * kr
                     src(if,jf,kf) = src(if,jf,kf) + r3 *
     &             (hxm1c * idir * (uc(iuc,juc,kc) + uc(iuc,juc,kc-1)) +
     &              hym1c * jdir * (vc(iuc,juc,kc) + vc(iuc,juc,kc-1)) +
     &              hzm1c *        (wc(iuc,juc,kc) - wc(iuc,juc,kc-1)))
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
                  juf = jf + ji - 1
                  fac0 = 1.d0 / (ir * kr)
                  ffac = kr * (ir - 1)
                  center = center + ffac
                  do l = 0, kr-1
                     fac1 = (kr-l) * fac0
                     if (l .eq. 0) fac1 = 0.5d0 * fac1
                     do m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        do kf = kr*cregl2, kr*cregh2, kr
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                       (hxm1 *
     &                (uf(if+m,juf,kf-l) - uf(if+m-1,juf,kf-l) +
     &                 uf(if+m,juf,kf-l-1) - uf(if+m-1,juf,kf-l-1) +
     &                 uf(if+m,juf,kf+l) - uf(if+m-1,juf,kf+l) +
     &                 uf(if+m,juf,kf+l-1) - uf(if+m-1,juf,kf+l-1)) +
     &                                        hym1 * jdir *
     &                (vf(if+m,juf,kf-l) + vf(if+m,juf,kf-l-1) +
     &                 vf(if+m-1,juf,kf-l) + vf(if+m-1,juf,kf-l-1) +
     &                 vf(if+m,juf,kf+l) + vf(if+m,juf,kf+l-1) +
     &                 vf(if+m-1,juf,kf+l) + vf(if+m-1,juf,kf+l-1)) +
     &                                        hzm1 *
     &                (wf(if+m,juf,kf-l) - wf(if+m,juf,kf-l-1) +
     &                 wf(if+m-1,juf,kf-l) - wf(if+m-1,juf,kf-l-1) +
     &                 wf(if+m,juf,kf+l) - wf(if+m,juf,kf+l-1) +
     &                 wf(if+m-1,juf,kf+l) - wf(if+m-1,juf,kf+l-1)))
                        end do
                     end do
                  end do
               end if
               if (ga(ii,ji,0) - ga(1-ii,ji,0) .eq. 1) then
                  iuf = if + ii - 1
                  fac0 = 1.d0 / (jr * kr)
                  ffac = kr * (jr - 1)
                  center = center + ffac
                  do l = 0, kr-1
                     fac1 = (kr-l) * fac0
                     if (l .eq. 0) fac1 = 0.5d0 * fac1
                     do n = jdir, jdir*(jr-1), jdir
                        fac = (jr-abs(n)) * fac1
                        do kf = kr*cregl2, kr*cregh2, kr
                           src(if,jf,kf) = src(if,jf,kf) + fac *
     &                                       (hxm1 * idir *
     &                (uf(iuf,jf+n,kf-l) + uf(iuf,jf+n-1,kf-l) +
     &                 uf(iuf,jf+n,kf-l-1) + uf(iuf,jf+n-1,kf-l-1) +
     &                 uf(iuf,jf+n,kf+l) + uf(iuf,jf+n-1,kf+l) +
     &                 uf(iuf,jf+n,kf+l-1) + uf(iuf,jf+n-1,kf+l-1)) +
     &                                        hym1 *
     &                (vf(iuf,jf+n,kf-l) - vf(iuf,jf+n-1,kf-l) +
     &                 vf(iuf,jf+n,kf-l-1) - vf(iuf,jf+n-1,kf-l-1) +
     &                 vf(iuf,jf+n,kf+l) - vf(iuf,jf+n-1,kf+l) +
     &                 vf(iuf,jf+n,kf+l-1) - vf(iuf,jf+n-1,kf+l-1)) +
     &                                        hzm1 *
     &                (wf(iuf,jf+n,kf-l) - wf(iuf,jf+n,kf-l-1) +
     &                 wf(iuf,jf+n-1,kf-l) - wf(iuf,jf+n-1,kf-l-1) +
     &                 wf(iuf,jf+n,kf+l) - wf(iuf,jf+n,kf+l-1) +
     &                 wf(iuf,jf+n-1,kf+l) - wf(iuf,jf+n-1,kf+l-1)))
                        end do
                     end do
                  end do
               end if
            end do
         end do
c weighting
         do kf = kr*cregl2, kr*cregh2, kr
            src(if,jf,kf) = src(if,jf,kf) / center
         end do
      end if
      end
c-----------------------------------------------------------------------
c Note---only generates values at coarse points along edge of fine grid
      subroutine hgcdiv(
     & src, srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & uc, vc, wc,
     &      cl0,ch0,cl1,ch1,cl2,ch2,
     & uf, vf, wf,
     &      fl0,fh0,fl1,fh1,fl2,fh2,
     &      cregl0,cregh0,cregl1,cregh1,cregl2,cregh2,
     & hx, hy, hz, ir, jr, kr, ga)
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer cl0,ch0,cl1,ch1,cl2,ch2
      integer fl0,fh0,fl1,fh1,fl2,fh2
      integer cregl0,cregh0,cregl1,cregh1,cregl2,cregh2
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2)
      double precision uc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision vc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision wc(cl0:ch0,cl1:ch1,cl2:ch2)
      double precision uf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision vf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision wf(fl0:fh0,fl1:fh1,fl2:fh2)
      double precision hx, hy, hz
      integer ir, jr, kr, ga(0:1,0:1,0:1)
      double precision r3, hxm1, hym1, hzm1, hxm1c, hym1c, hzm1c
      double precision sum, center, cfac, ffac, fac0, fac1, fac
      integer ic, jc, kc, if, jf, kf, iuc, iuf, juc, juf, kuc, kuf
      integer ii, ji, ki, idir, jdir, kdir, l, m, n
      r3 = ir * jr * kr
      hxm1c = 1.d0 / (ir * hx)
      hym1c = 1.d0 / (jr * hy)
      hzm1c = 1.d0 / (kr * hz)
      hxm1 = ir * hxm1c
      hym1 = jr * hym1c
      hzm1 = kr * hzm1c
      ic = cregl0
      jc = cregl1
      kc = cregl2
      if = ic * ir
      jf = jc * jr
      kf = kc * kr
      sum = 0.d0
      center = 0.d0
c octants
      fac = 1.d0
      ffac = 0.5d0
      cfac = 0.5d0 * r3
      do ki = 0, 1
         kdir = 2 * ki - 1
         kuf = kf + ki - 1
         kuc = kc + ki - 1
         do ji = 0, 1
            jdir = 2 * ji - 1
            juf = jf + ji - 1
            juc = jc + ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,ki) .eq. 1) then
                  iuf = if + ii - 1
                  center = center + ffac
                  sum = sum + fac *
     &              (hxm1 * idir * uf(iuf,juf,kuf) +
     &               hym1 * jdir * vf(iuf,juf,kuf) +
     &               hzm1 * kdir * wf(iuf,juf,kuf))
               else
                  iuc = ic + ii - 1
                  center = center + cfac
                  sum = sum + r3 *
     &              (hxm1c * idir * uc(iuc,juc,kuc) +
     &               hym1c * jdir * vc(iuc,juc,kuc) +
     &               hzm1c * kdir * wc(iuc,juc,kuc))
               end if
            end do
         end do
      end do
c faces
      do ki = 0, 1
         kdir = 2 * ki - 1
         kuf = kf + ki - 1
         do ji = 0, 1
            jdir = 2 * ji - 1
            juf = jf + ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               iuf = if + ii - 1
               if (ga(ii,ji,ki) - ga(ii,ji,1-ki) .eq. 1) then
                  fac0 = 1.d0 / (ir * jr)
                  ffac = 0.5d0 * (ir-1) * (jr-1)
                  center = center + ffac
                  do n = jdir, jdir*(jr-1), jdir
                     fac1 = (jr-abs(n)) * fac0
                     do m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        sum = sum + fac * (hxm1 *
     &               (uf(if+m,jf+n,kuf) - uf(if+m-1,jf+n,kuf) +
     &                uf(if+m,jf+n-1,kuf) - uf(if+m-1,jf+n-1,kuf)) +
     &                                     hym1 *
     &               (vf(if+m,jf+n,kuf) - vf(if+m,jf+n-1,kuf) +
     &                vf(if+m-1,jf+n,kuf) - vf(if+m-1,jf+n-1,kuf)) +
     &                                     hzm1 * kdir *
     &               (wf(if+m,jf+n,kuf) + wf(if+m,jf+n-1,kuf) +
     &                wf(if+m-1,jf+n,kuf) + wf(if+m-1,jf+n-1,kuf)))
                     end do
                  end do
               end if
               if (ga(ii,ji,ki) - ga(ii,1-ji,ki) .eq. 1) then
                  fac0 = 1.d0 / (ir * kr)
                  ffac = 0.5d0 * (ir-1) * (kr-1)
                  center = center + ffac
                  do l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        sum = sum + fac * (hxm1 *
     &               (uf(if+m,juf,kf+l) - uf(if+m-1,juf,kf+l) +
     &                uf(if+m,juf,kf+l-1) - uf(if+m-1,juf,kf+l-1)) +
     &                                     hym1 * jdir *
     &               (vf(if+m,juf,kf+l) + vf(if+m,juf,kf+l-1) +
     &                vf(if+m-1,juf,kf+l) + vf(if+m-1,juf,kf+l-1)) +
     &                                     hzm1 *
     &               (wf(if+m,juf,kf+l) - wf(if+m,juf,kf+l-1) +
     &                wf(if+m-1,juf,kf+l) - wf(if+m-1,juf,kf+l-1)))
                     end do
                  end do
               end if
               if (ga(ii,ji,ki) - ga(1-ii,ji,ki) .eq. 1) then
                  fac0 = 1.d0 / (jr * kr)
                  ffac = 0.5d0 * (jr-1) * (kr-1)
                  center = center + ffac
                  do l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do n = jdir, jdir*(jr-1), jdir
                        fac = (jr-abs(n)) * fac1
                        sum = sum + fac * (hxm1 * idir *
     &               (uf(iuf,jf+n,kf+l) + uf(iuf,jf+n,kf+l-1) +
     &                uf(iuf,jf+n-1,kf+l) + uf(iuf,jf+n-1,kf+l-1)) +
     &                                     hym1 *
     &               (vf(iuf,jf+n,kf+l) - vf(iuf,jf+n-1,kf+l) +
     &                vf(iuf,jf+n,kf+l-1) - vf(iuf,jf+n-1,kf+l-1)) +
     &                                     hzm1 *
     &               (wf(iuf,jf+n,kf+l) - wf(iuf,jf+n,kf+l-1) +
     &                wf(iuf,jf+n-1,kf+l) - wf(iuf,jf+n-1,kf+l-1)))
                     end do
                  end do
               end if
            end do
         end do
      end do
c edges
      do ki = 0, 1
         kdir = 2 * ki - 1
         kuf = kf + ki - 1
         do ji = 0, 1
            jdir = 2 * ji - 1
            juf = jf + ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               iuf = if + ii - 1
               if (ga(ii,ji,ki) -
     &             min(ga(ii,ji,1-ki), ga(ii,1-ji,ki), ga(ii,1-ji,1-ki))
     &             .eq. 1) then
                  fac1 = 1.d0 / ir
                  ffac = 0.5d0 * (ir-1)
                  center = center + ffac
                  do m = idir, idir*(ir-1), idir
                     fac = (ir-abs(m)) * fac1
                     sum = sum + fac * (hxm1 *
     &                 (uf(if+m,juf,kuf) - uf(if+m-1,juf,kuf)) +
     &                                  hym1 * jdir *
     &                 (vf(if+m,juf,kuf) + vf(if+m-1,juf,kuf)) +
     &                                  hzm1 * kdir *
     &                 (wf(if+m,juf,kuf) + wf(if+m-1,juf,kuf)))
                  end do
               end if
               if (ga(ii,ji,ki) -
     &             min(ga(ii,ji,1-ki), ga(1-ii,ji,ki), ga(1-ii,ji,1-ki))
     &             .eq. 1) then
                  fac1 = 1.d0 / jr
                  ffac = 0.5d0 * (jr-1)
                  center = center + ffac
                  do n = jdir, jdir*(jr-1), jdir
                     fac = (jr-abs(n)) * fac1
                     sum = sum + fac * (hxm1 * idir *
     &                 (uf(iuf,jf+n,kuf) + uf(iuf,jf+n-1,kuf)) +
     &                                  hym1 *
     &                 (vf(iuf,jf+n,kuf) - vf(iuf,jf+n-1,kuf)) +
     &                                  hzm1 * kdir *
     &                 (wf(iuf,jf+n,kuf) + wf(iuf,jf+n-1,kuf)))
                  end do
               end if
               if (ga(ii,ji,ki) -
     &             min(ga(ii,1-ji,ki), ga(1-ii,ji,ki), ga(1-ii,1-ji,ki))
     &             .eq. 1) then
                  fac1 = 1.d0 / kr
                  ffac = 0.5d0 * (kr-1)
                  center = center + ffac
                  do l = kdir, kdir*(kr-1), kdir
                     fac = (kr-abs(l)) * fac1
                     sum = sum + fac * (hxm1 * idir *
     &                 (uf(iuf,juf,kf+l) + uf(iuf,juf,kf+l-1)) +
     &                                  hym1 * jdir *
     &                 (vf(iuf,juf,kf+l) + vf(iuf,juf,kf+l-1)) +
     &                                  hzm1 *
     &                 (wf(iuf,juf,kf+l) - wf(iuf,juf,kf+l-1)))
                  end do
               end if
            end do
         end do
      end do
c weighting
      src(if,jf,kf) = sum / center
      end
