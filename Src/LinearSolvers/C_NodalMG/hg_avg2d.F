
c#include "hg_version.H"

#define dims(a) a l0,a h0,a l1,a h1
#define dimdec(a) a l0:a h0,a l1:a h1
c#define dimdecu(a) (a h0-a l0+1)*(a h1-a l1+1)
#define dimdecu(a) *

      subroutine hgavg(src,
     @ dims(src),
     @ rf,
     @ dims(f),
     @ dims(freg), hx, irz, imax)
      integer dims(src)
      integer dims(f)
      integer dims(freg)
      real*8 src(dimdec(src))
      real*8 rf(dimdec(f))
      real*8 hx
      integer irz, imax
      real*8 fac, r1, r0m, r1m
      integer i, j
      fac = 0.25d0
      do 10 j = fregl1, fregh1
         do 10 i = fregl0, fregh0
 10         src(i,j) = src(i,j) + fac * (rf(i-1,j-1) + rf(i-1,j) +
     @                                   rf(i,j-1)   + rf(i,j))
      if (irz .eq. 1) then
         fac = hx / 24.d0
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
     @            (r0m * (rf(i-1,j) + rf(i-1,j-1)) -
     @             r1m * (rf(i,j-1) + rf(i,j)))
      endif
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgfavg(src,
     @ dims(src),
     @ rc,
     @ dims(c),
     @ rf,
     @ dims(f),
     @ dims(creg),
     @ ir, jr, idim, idir, hx, irz, imax)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 rc(dimdec(c))
      real*8 rf(dimdec(f))
      real*8 hx
      integer ir, jr, idim, idir, irz, imax
      real*8 fac0, fac, r, rfac, rfac0, rfac1
      real*8 rfac0m, rfac1m, rfac0p, rfac1p
      integer i, j, irc, jrc, irf, jrf, m, n
      if (idim .eq. 0) then
         i = cregl0
         if (idir .eq. 1) then
            irc = i - 1
            irf = i * ir
         else
            irc = i
            irf = i * ir - 1
         endif
         fac0 = 0.5d0 * ir / (ir+1)
         do 10 j = cregl1, cregh1
 10         src(i*ir,j*jr) = src(i*ir,j*jr) + fac0 *
     @         (rc(irc,j) + rc(irc,j-1))
         if (irz .eq. 1) then
            r = (irc + 0.5d0) * (hx * ir)
            rfac = (hx * ir) / (6.d0 * r)
            do 20 j = cregl1, cregh1
 20            src(i*ir,j*jr) = src(i*ir,j*jr) + fac0 *
     @           (rfac * idir * (rc(irc,j) + rc(irc,j-1)))
            r = (irf + 0.5d0) * hx
            rfac = hx / (6.d0 * r)
         endif
         fac0 = fac0 / (ir * jr * jr)
         i = i * ir
         do 50 n = 0, jr-1
            fac = (jr-n) * fac0
            if (n .eq. 0) fac = 0.5d0 * fac
            do 30 j = jr*cregl1, jr*cregh1, jr
 30            src(i,j) = src(i,j) + fac *
     @           (rf(irf,j-n) + rf(irf,j-n-1) +
     @            rf(irf,j+n) + rf(irf,j+n-1))
            if (irz .eq. 1) then
               do 40 j = jr*cregl1, jr*cregh1, jr
 40               src(i,j) = src(i,j) - fac *
     @              (rfac * idir * (rf(irf,j-n) + rf(irf,j-n-1) +
     @                              rf(irf,j+n) + rf(irf,j+n-1)))
            endif
 50         continue
      else if (idim .eq. 1) then
         j = cregl1
         if (idir .eq. 1) then
            jrc = j - 1
            jrf = j * jr
         else
            jrc = j
            jrf = j * jr - 1
         endif
         fac0 = 0.5d0 * jr / (jr+1)
         do 60 i = cregl0, cregh0
 60         src(i*ir,j*jr) = src(i*ir,j*jr) + fac0 *
     @        (rc(i,jrc) + rc(i-1,jrc))
         if (irz .eq. 1 .and. cregh0 .lt. imax) then
            do 70 i = cregl0, cregh0
               r = (i + 0.5d0) * (hx * ir)
               rfac0 = (hx * ir) / (6.d0 * (r - hx * ir))
               rfac1 = (hx * ir) / (6.d0 * r)
 70            src(i*ir,j*jr) = src(i*ir,j*jr) + fac0 *
     @           (rfac0 * rc(i-1,jrc) - rfac1 * rc(i,jrc))
         else if (irz .eq. 1) then
c This should only occur with a corner at the outer boundary, which
c should be handled by the cavg routine instead:
            i = cregh0
            r = (i - 0.5d0) * (hx * ir)
            rfac0 = (hx * ir) / (6.d0 * r)
            rfac1 = -rfac0
            src(i*ir,j*jr) = src(i*ir,j*jr) + fac0 *
     @           (rfac0 * rc(i-1,jrc) - rfac1 * rc(i,jrc))
         endif
         fac0 = fac0 / (ir * ir * jr)
         j = j * jr
         do 100 m = 0, ir-1
            fac = (ir-m) * fac0
            if (m .eq. 0) fac = 0.5d0 * fac
            do 80 i = ir*cregl0, ir*cregh0, ir
 80            src(i,j) = src(i,j) + fac *
     @           (rf(i-m,jrf) + rf(i-m-1,jrf) +
     @            rf(i+m,jrf) + rf(i+m-1,jrf))
            if (irz .eq. 1 .and. cregh0 .lt. imax) then
               do 90 i = ir*cregl0, ir*cregh0, ir
                  r = (i + 0.5d0) * hx
                  rfac0m = hx / (6.d0 * (r - (m + 1) * hx))
                  rfac1m = hx / (6.d0 * (r - m * hx))
                  rfac0p = hx / (6.d0 * (r + (m - 1) * hx))
                  rfac1p = hx / (6.d0 * (r + m * hx))
 90               src(i,j) = src(i,j) + fac *
     @              (rfac0m * rf(i-m-1,jrf) - rfac1m * rf(i-m,jrf) +
     @               rfac0p * rf(i+m-1,jrf) - rfac1p * rf(i+m,jrf))
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
               endif
               rfac0p = -rfac1m
               rfac1p = -rfac0m
               src(i,j) = src(i,j) + fac *
     @              (rfac0m * rf(i-m-1,jrf) - rfac1m * rf(i-m,jrf) +
     @               rfac0p * rf(i+m-1,jrf) - rfac1p * rf(i+m,jrf))
            endif
 100        continue
      endif
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgcavg(src,
     @ dims(src),
     @ rc,
     @ dims(c),
     @ rf,
     @ dims(f),
     @ dims(creg),
     @ ir, jr, ga, hx, irz, imax)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 rc(dimdec(c))
      real*8 rf(dimdec(f))
      real*8 hx
      integer ir, jr, ga(0:1,0:1), irz, imax
      real*8 rm2, sum, center, fac, ffac, cfac, fac1
      real*8 r, rfac, rfac0, rfac1
      integer ic, jc, if, jf, ii, ji, irc, jrc, irf, jrf
      integer m, n, idir, jdir
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
      do 10 ji = 0, 1
         jrf = jf + ji - 1
         jrc = jc + ji - 1
         do 10 ii = 0, 1
            if (ga(ii,ji) .eq. 1) then
               irf = if + ii - 1
               center = center + ffac
               sum = sum + fac * rf(irf,jrf)
               if (irz .eq. 1) then
                  idir = 2 * ii - 1
                  r = (irf + 0.5d0) * hx
                  if (irf .lt. (ir * imax)) then
                     rfac =  hx / (6.d0 * r)
                  else
                     rfac = -hx / (6.d0 * (r - hx))
                  endif
                  sum = sum - fac * rfac * idir * rf(irf,jrf)
               endif
            else
               irc = ic + ii - 1
               center = center + cfac
               sum = sum + rc(irc,jrc)
               if (irz .eq. 1) then
                  idir = 2 * ii - 1
                  r = (irc + 0.5d0) * (hx * ir)
                  if (irc .lt. imax) then
                     rfac =  (hx * ir) / (6.d0 * r)
                  else
                     rfac = -(hx * ir) / (6.d0 * (r - (hx * ir)))
                  endif
                  sum = sum - rfac * idir * rc(irc,jrc)
               endif
            endif
 10         continue
c edges
      do 40 ji = 0, 1
         jdir = 2 * ji - 1
         jrf = jf + ji - 1
         do 40 ii = 0, 1
            idir = 2 * ii - 1
            irf = if + ii - 1
            if (ga(ii,ji) - ga(ii,1-ji) .eq. 1) then
               fac1 = rm2 / ir
               ffac = (ir-1) * rm2
               center = center + ffac
               do 20 m = idir, idir*(ir-1), idir
                  fac = (ir-abs(m)) * fac1
                  sum = sum + fac * (rf(if+m,jrf) + rf(if+m-1,jrf))
                  if (irz .eq. 1) then
                     r = (if + 0.5d0) * hx
                     if (irf .lt. (ir * imax)) then
                        rfac0 =  hx / (6.d0 * (r + (m - 1) * hx))
                        rfac1 =  hx / (6.d0 * (r + m * hx))
                     else
                        rfac0 = -hx / (6.d0 * (r - m * hx))
                        rfac1 = -hx / (6.d0 * (r - (m + 1) * hx))
                     endif
                     sum = sum + fac *
     @                 (rfac0 * rf(if+m-1,jrf) - rfac1 * rf(if+m,jrf))
                  endif
 20               continue
            endif
            if (ga(ii,ji) - ga(1-ii,ji) .eq. 1) then
               fac1 = rm2 / jr
               ffac = (jr-1) * rm2
               center = center + ffac
               do 30 n = jdir, jdir*(jr-1), jdir
                  fac = (jr-abs(n)) * fac1
                  sum = sum + fac * (rf(irf,jf+n) + rf(irf,jf+n-1))
                  if (irz .eq. 1) then
                     r = (irf + 0.5d0) * hx
                     rfac = hx / (6.d0 * r)
                     sum = sum - fac *
     @                 (rfac * idir * (rf(irf,jf+n) + rf(irf,jf+n-1)))
                  endif
 30               continue
            endif
 40         continue
c weighting
      src(if,jf) = src(if,jf) + sum / center
      end
