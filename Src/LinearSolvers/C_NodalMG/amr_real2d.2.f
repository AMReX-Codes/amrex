c-----------------------------------------------------------------------
c NODE-based data only.
c Fills coarse region defined by reg.
c Handles any corner geometry except all-coarse.
      subroutine ancr2(
     & dest, destl0,desth0,destl1,desth1,
     &       regl0,regh0,regl1,regh1,
     & src,  srcl0,srch0,srcl1,srch1,
     & ir, jr, ncomp, integ, ga, i2)
      integer ncomp
      integer destl0,desth0,destl1,desth1
      integer regl0,regh0,regl1,regh1
      integer srcl0,srch0,srcl1,srch1
      integer ir, jr, ga(0:1,0:1), integ
      double precision dest(destl0:desth0,destl1:desth1,
     &   ncomp)
      double precision src(srcl0:srch0,srcl1:srch1, ncomp)
      double precision cube, center, cfac, fac0, fac1, fac
      integer i, j, ii, ji, idir, jdir, m, n, nc, i2

      i = regl0
      j = regl1
      do nc = 1, ncomp
      dest(i,j,nc) = 0.0D0
      cube = ir * jr
      if (integ .eq. 0) then
         center = 1.0D0 / cube
         fac0 = 1.0D0 / (cube**2)
         cfac = 0.25D0 * cube * fac0 * (ir-1) * (jr-1)
      else
         center = 1.0D0
         fac0 = 1.0D0 / cube
      end if

c octants
         do ji = 0, 1
            jdir = 2 * ji - 1
            do ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji) .eq. 1) then
                     do n = jdir, jdir*(jr-1), jdir
                        fac1 = (jr-abs(n)) * fac0
                        do m = idir, idir*(ir-1), idir
                           fac = (ir-abs(m)) * fac1
                           dest(i,j,nc) = dest(i,j,nc) +
     &                       fac * src(i*ir+m,j*jr+n,nc)
                        end do
                     end do
               else if (integ .eq. 0) then
                  center = center + cfac
               end if
            end do
          end do
c faces
      fac1 = jr * fac0
      cfac = 0.50D0 * cube * fac0 * (ir-1)
         do ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,0) + ga(ii,1) .eq. 2) then
                  do m = idir, idir*(ir-1), idir
                     fac = (ir-abs(m)) * fac1
                     dest(i,j,nc) = dest(i,j,nc) +
     &                 fac * src(i*ir+m,j*jr,nc)
                  end do
            else if (integ .eq. 0) then
               center = center + cfac
            end if
         end do
      fac1 = ir * fac0
      cfac = 0.50D0 * cube * fac0 * (jr-1)
         do ji = 0, 1
            jdir = 2 * ji - 1
            if (ga(0,ji) + ga(1,ji) .eq. 2) then
                  do n = jdir, jdir*(jr-1), jdir
                     fac = (jr-abs(n)) * fac1
                     dest(i,j,nc) = dest(i,j,nc) +
     &                 fac * src(i*ir,j*jr+n,nc)
                  end do
            else if (integ .eq. 0) then
               center = center + cfac
            end if
         end do
c center
      dest(i,j,nc) = dest(i,j,nc) + center * src(i*ir,j*jr,nc)
      end do

      end
