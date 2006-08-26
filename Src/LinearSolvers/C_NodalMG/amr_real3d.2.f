
c-----------------------------------------------------------------------
c NODE-based data only.
c Fills coarse region defined by reg.
c Handles any corner geometry except all-coarse.
      subroutine ancr2(
     & dest, destl0,desth0,destl1,desth1,destl2,desth2,
     &       regl0,regh0,regl1,regh1,regl2,regh2,
     & src,  srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & ir, jr, kr, ncomp, integ, ga, i2)
      integer ncomp
      integer destl0,desth0,destl1,desth1,destl2,desth2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer ir, jr, kr, ga(0:1,0:1,0:1), integ
      double precision dest(destl0:desth0,destl1:desth1,destl2:desth2,
     &   ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2, ncomp)
      double precision cube, center, cfac, fac0, fac1, fac2, fac
      integer i, j, k, ii, ji, ki, idir, jdir, kdir, l, m, n,nc, i2
      i = regl0
      j = regl1
      k = regl2
      do nc = 1, ncomp
      dest(i,j,k,nc) = 0.0D0
      cube = ir * jr * kr
      if (integ .eq. 0) then
         center = 1.0D0 / cube
         fac0 = 1.0D0 / (cube**2)
         cfac = 0.1250D0 * cube * fac0 * (ir-1) * (jr-1) * (kr-1)
      else
         center = 1.0D0
         fac0 = 1.0D0 / cube
      end if
c octants
      do ki = 0, 1
         kdir = 2 * ki - 1
         do ji = 0, 1
            jdir = 2 * ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,ki) .eq. 1) then
                  do l = kdir, kdir*(kr-1), kdir
                     fac2 = (kr-abs(l)) * fac0
                     do n = jdir, jdir*(jr-1), jdir
                        fac1 = (jr-abs(n)) * fac2
                        do m = idir, idir*(ir-1), idir
                           fac = (ir-abs(m)) * fac1
                           dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                       fac * src(i*ir+m,j*jr+n,k*kr+l,nc)
                        end do
                     end do
                  end do
               else if (integ .eq. 0) then
                  center = center + cfac
               end if
            end do
         end do
      end do
c faces
      fac2 = kr * fac0
      cfac = 0.25D0 * cube * fac0 * (ir-1) * (jr-1)
      do ji = 0, 1
         jdir = 2 * ji - 1
         do ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,ji,0) + ga(ii,ji,1) .eq. 2) then
               do n = jdir, jdir*(jr-1), jdir
                  fac1 = (jr-abs(n)) * fac2
                  do m = idir, idir*(ir-1), idir
                     fac = (ir-abs(m)) * fac1
                     dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                 fac * src(i*ir+m,j*jr+n,k*kr,nc)
                  end do
               end do
            else if (integ .eq. 0) then
               center = center + cfac
            end if
         end do
      end do
      fac2 = jr * fac0
      cfac = 0.25D0 * cube * fac0 * (ir-1) * (kr-1)
      do ki = 0, 1
         kdir = 2 * ki - 1
         do ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,0,ki) + ga(ii,1,ki) .eq. 2) then
               do l = kdir, kdir*(kr-1), kdir
                  fac1 = (kr-abs(l)) * fac2
                  do m = idir, idir*(ir-1), idir
                     fac = (ir-abs(m)) * fac1
                     dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                 fac * src(i*ir+m,j*jr,k*kr+l,nc)
                  end do
               end do
            else if (integ .eq. 0) then
               center = center + cfac
            end if
         end do
      end do
      fac2 = ir * fac0
      cfac = 0.25D0 * cube * fac0 * (jr-1) * (kr-1)
      do ki = 0, 1
         kdir = 2 * ki - 1
         do ji = 0, 1
            jdir = 2 * ji - 1
            if (ga(0,ji,ki) + ga(1,ji,ki) .eq. 2) then
               do l = kdir, kdir*(kr-1), kdir
                  fac1 = (kr-abs(l)) * fac2
                  do n = jdir, jdir*(jr-1), jdir
                     fac = (jr-abs(n)) * fac1
                     dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                 fac * src(i*ir,j*jr+n,k*kr+l,nc)
                  end do
               end do
            else if (integ .eq. 0) then
               center = center + cfac
            end if
         end do
      end do
c edges
      fac1 = jr * kr * fac0
      cfac = 0.5D0 * cube * fac0 * (ir-1)
      do ii = 0, 1
         idir = 2 * ii - 1
         if (ga(ii,0,0) + ga(ii,0,1) +
     &       ga(ii,1,0) + ga(ii,1,1) .eq. 4) then
            do m = idir, idir*(ir-1), idir
               fac = (ir-abs(m)) * fac1
               dest(i,j,k,nc) = dest(i,j,k,nc) +
     &           fac * src(i*ir+m,j*jr,k*kr,nc)
            end do
         else if (integ .eq. 0) then
            center = center + cfac
         end if
      end do
      fac1 = ir * kr * fac0
      cfac = 0.5D0 * cube * fac0 * (jr-1)
      do ji = 0, 1
         jdir = 2 * ji - 1
         if (ga(0,ji,0) + ga(0,ji,1) +
     &       ga(1,ji,0) + ga(1,ji,1) .eq. 4) then
            do n = jdir, jdir*(jr-1), jdir
               fac = (jr-abs(n)) * fac1
               dest(i,j,k,nc) = dest(i,j,k,nc) +
     &            fac * src(i*ir,j*jr+n,k*kr,nc)
            end do
         else if (integ .eq. 0) then
            center = center + cfac
         end if
      end do
      fac1 = ir * jr * fac0
      cfac = 0.5D0 * cube * fac0 * (kr-1)
      do ki = 0, 1
         kdir = 2 * ki - 1
         if (ga(0,0,ki) + ga(0,1,ki) +
     &       ga(1,0,ki) + ga(1,1,ki) .eq. 4) then
            do l = kdir, kdir*(kr-1), kdir
               fac = (kr-abs(l)) * fac1
               dest(i,j,k,nc) = dest(i,j,k,nc) +
     &           fac * src(i*ir,j*jr,k*kr+l,nc)
            end do
         else if (integ .eq. 0) then
            center = center + cfac
         end if
      end do
c center
      dest(i,j,k,nc) = dest(i,j,k,nc) +
     &  center * src(i*ir,j*jr,k*kr,nc)
      end do
      end
c-----------------------------------------------------------------------
c NODE-based data only.
c Fills coarse region defined by reg.
c Handles any edge geometry except all-coarse or all-fine.
      subroutine aner2(
     & dest, destl0,desth0,destl1,desth1,destl2,desth2,
     &       regl0,regh0,regl1,regh1,regl2,regh2,
     & src,  srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & ir, jr, kr, ncomp, integ, ivect, ga)
      integer ncomp
      integer destl0,desth0,destl1,desth1,destl2,desth2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer ir, jr, kr, ivect(0:2), ga(0:1,0:1,0:1), integ
      double precision dest(destl0:desth0,destl1:desth1,destl2:desth2,
     &   ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2,ncomp)
      double precision cube, center, cfac, fac0, fac1, fac2, fac
      integer i, j, k, ii, ji, ki, idir, jdir, kdir, l, m, n, nc
      cube = ir * jr * kr
      if (ivect(0) .eq. 0) then
         j = regl1
         k = regl2
         do nc = 1, ncomp
         do i = regl0, regh0
            dest(i,j,k,nc) = 0.0D0
         end do
c center gets center plus two edges
         if (integ .eq. 0) then
            center = ir / cube
            fac0 = 1.0D0 / (cube**2)
            cfac = 0.25D0 * cube * fac0 * ir * (jr-1) * (kr-1)
         else
            center = 1.0D0
            fac0 = 1.0D0 / cube
         end if
c quadrants
c each quadrant is two octants and a face
         do ki = 0, 1
            kdir = 2 * ki - 1
            do ji = 0, 1
               jdir = 2 * ji - 1
               if (ga(0,ji,ki) .eq. 1) then
                  do l = kdir, kdir*(kr-1), kdir
                     fac2 = (kr-abs(l)) * fac0
                     do n = jdir, jdir*(jr-1), jdir
                        fac1 = (jr-abs(n)) * fac2
                        do m = 0, ir-1
                           fac = (ir-m) * fac1
                           if (m .eq. 0) fac = 0.5D0 * fac
                           do i = regl0, regh0
                              dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                          fac * (src(i*ir-m,j*jr+n,k*kr+l,nc) +
     &                                 src(i*ir+m,j*jr+n,k*kr+l,nc))
                           end do
                        end do
                     end do
                  end do
               else if (integ .eq. 0) then
                  center = center + cfac
               end if
            end do
         end do
c faces
c each face is two faces and an edge
         fac2 = kr * fac0
         cfac = 0.5D0 * cube * fac0 * ir * (jr-1)
         do ji = 0, 1
            jdir = 2 * ji - 1
            if (ga(0,ji,0) + ga(0,ji,1) .eq. 2) then
               do n = jdir, jdir*(jr-1), jdir
                  fac1 = (jr-abs(n)) * fac2
                  do m = 0, ir-1
                     fac = (ir-m) * fac1
                     if (m .eq. 0) fac = 0.5D0 * fac
                     do i = regl0, regh0
                        dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                    fac * (src(i*ir-m,j*jr+n,k*kr,nc) +
     &                           src(i*ir+m,j*jr+n,k*kr,nc))
                     end do
                  end do
               end do
            else if (integ .eq. 0) then
               center = center + cfac
            end if
         end do
         fac2 = jr * fac0
         cfac = 0.5D0 * cube * fac0 * ir * (kr-1)
         do ki = 0, 1
            kdir = 2 * ki - 1
            if (ga(0,0,ki) + ga(0,1,ki) .eq. 2) then
               do l = kdir, kdir*(kr-1), kdir
                  fac1 = (kr-abs(l)) * fac2
                  do m = 0, ir-1
                     fac = (ir-m) * fac1
                     if (m .eq. 0) fac = 0.5D0 * fac
                     do i = regl0, regh0
                        dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                    fac * (src(i*ir-m,j*jr,k*kr+l,nc) +
     &                           src(i*ir+m,j*jr,k*kr+l,nc))
                     end do
                  end do
               end do
            else if (integ .eq. 0) then
               center = center + cfac
            end if
         end do
c center
         do i = regl0, regh0
            dest(i,j,k,nc) = dest(i,j,k,nc) +
     &        center * src(i*ir,j*jr,k*kr,nc)
         end do
         end do
      else if (ivect(1) .eq. 0) then
         i = regl0
         k = regl2
         do nc = 1, ncomp
         do  j = regl1, regh1
            dest(i,j,k,nc) = 0.0D0
         end do
c center gets center plus two edges
         if (integ .eq. 0) then
            center = jr / cube
            fac0 = 1.0D0 / (cube**2)
            cfac = 0.25D0 * cube * fac0 * jr * (ir-1) * (kr-1)
         else
            center = 1.0D0
            fac0 = 1.0D0 / cube
         end if
c quadrants
c each quadrant is two octants and a face
         do  ki = 0, 1
            kdir = 2 * ki - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,0,ki) .eq. 1) then
                  do l = kdir, kdir*(kr-1), kdir
                     fac2 = (kr-abs(l)) * fac0
                     do n = 0, jr-1
                        fac1 = (jr-n) * fac2
                        if (n .eq. 0) fac1 = 0.5D0 * fac1
                        do m = idir, idir*(ir-1), idir
                           fac = (ir-abs(m)) * fac1
                           do j = regl1, regh1
                              dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                          fac * (src(i*ir+m,j*jr-n,k*kr+l,nc) +
     &                                 src(i*ir+m,j*jr+n,k*kr+l,nc))
                           end do
                        end do
                     end do
                  end do
               else if (integ .eq. 0) then
                  center = center + cfac
               end if
            end do
         end do
c faces
c each face is two faces and an edge
         fac2 = kr * fac0
         cfac = 0.5D0 * cube * fac0 * jr * (ir-1)
         do ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,0,0) + ga(ii,0,1) .eq. 2) then
               do n = 0, jr-1
                  fac1 = (jr-n) * fac2
                  if (n .eq. 0) fac1 = 0.5D0 * fac1
                  do m = idir, idir*(ir-1), idir
                     fac = (ir-abs(m)) * fac1
                     do j = regl1, regh1
                        dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                    fac * (src(i*ir+m,j*jr-n,k*kr,nc) +
     &                           src(i*ir+m,j*jr+n,k*kr,nc))
                     end do
                  end do
               end do
            else if (integ .eq. 0) then
               center = center + cfac
            end if
         end do
         fac2 = ir * fac0
         cfac = 0.5D0 * cube * fac0 * jr * (kr-1)
         do ki = 0, 1
            kdir = 2 * ki - 1
            if (ga(0,0,ki) + ga(1,0,ki) .eq. 2) then
               do  l = kdir, kdir*(kr-1), kdir
                  fac1 = (kr-abs(l)) * fac2
                  do n = 0, jr-1
                     fac = (jr-n) * fac1
                     if (n .eq. 0) fac = 0.5D0 * fac
                     do j = regl1, regh1
                        dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                    fac * (src(i*ir,j*jr-n,k*kr+l,nc) +
     &                           src(i*ir,j*jr+n,k*kr+l,nc))
                     end do
                  end do
               end do
            else if (integ .eq. 0) then
               center = center + cfac
            end if
         end do
c center
         do j = regl1, regh1
            dest(i,j,k,nc) = dest(i,j,k,nc) +
     &        center * src(i*ir,j*jr,k*kr,nc)
         end do
         end do
      else
         i = regl0
         j = regl1
         do nc = 1, ncomp
         do k = regl2, regh2
            dest(i,j,k,nc) = 0.0D0
         end do
c center gets center plus two edges
         if (integ .eq. 0) then
            center = kr / cube
            fac0 = 1.0D0 / (cube**2)
            cfac = 0.25D0 * cube * fac0 * kr * (ir-1) * (jr-1)
         else
            center = 1.0D0
            fac0 = 1.0D0 / cube
         end if
c quadrants
c each quadrant is two octants and a face
         do ji = 0, 1
            jdir = 2 * ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,0) .eq. 1) then
                  do l = 0, kr-1
                     fac2 = (kr-l) * fac0
                     if (l .eq. 0) fac2 = 0.5D0 * fac2
                     do n = jdir, jdir*(jr-1), jdir
                        fac1 = (jr-abs(n)) * fac2
                        do m = idir, idir*(ir-1), idir
                           fac = (ir-abs(m)) * fac1
                           do k = regl2, regh2
                              dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                          fac * (src(i*ir+m,j*jr+n,k*kr-l,nc) +
     &                                 src(i*ir+m,j*jr+n,k*kr+l,nc))
                           end do
                        end do
                     end do
                  end do
               else if (integ .eq. 0) then
                  center = center + cfac
               end if
            end do
         end do
c faces
c each face is two faces and an edge
         fac2 = jr * fac0
         cfac = 0.5D0 * cube * fac0 * kr * (ir-1)
         do ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,0,0) + ga(ii,1,0) .eq. 2) then
               do l = 0, kr-1
                  fac1 = (kr-l) * fac2
                  if (l .eq. 0) fac1 = 0.5D0 * fac1
                  do m = idir, idir*(ir-1), idir
                     fac = (ir-abs(m)) * fac1
                     do k = regl2, regh2
                        dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                    fac * (src(i*ir+m,j*jr,k*kr-l,nc) +
     &                           src(i*ir+m,j*jr,k*kr+l,nc))
                     end do
                  end do
               end do
            else if (integ .eq. 0) then
               center = center + cfac
            end if
         end do
         fac2 = ir * fac0
         cfac = 0.5D0 * cube * fac0 * kr * (jr-1)
         do ji = 0, 1
            jdir = 2 * ji - 1
            if (ga(0,ji,0) + ga(1,ji,0) .eq. 2) then
               do l = 0, kr-1
                  fac1 = (kr-l) * fac2
                  if (l .eq. 0) fac1 = 0.5D0 * fac1
                  do n = jdir, jdir*(jr-1), jdir
                     fac = (jr-abs(n)) * fac1
                     do k = regl2, regh2
                        dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                    fac * (src(i*ir,j*jr+n,k*kr-l,nc) +
     &                           src(i*ir,j*jr+n,k*kr+l,nc))
                     end do
                  end do
               end do
            else if (integ .eq. 0) then
               center = center + cfac
            end if
            end do
c center
         do k = regl2, regh2
            dest(i,j,k,nc) = dest(i,j,k,nc) +
     &        center * src(i*ir,j*jr,k*kr,nc)
         end do
         end do
      end if
      end
