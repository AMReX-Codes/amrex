
#define dims(a) a l0,a h0,a l1,a h1,a l2,a h2
#define dimdec(a) a l0:a h0,a l1:a h1,a l2:a h2
c#define dimdecu(a) (a h0-a l0+1)*(a h1-a l1+1)*(a h2-a l2+1)
#define dimdecu(a) *

      subroutine hgavg(src,
     @ dims(src),
     @ rf,
     @ dims(f),
     @ dims(freg))
      integer dims(src)
      integer dims(f)
      integer dims(freg)
      real*8 src(dimdec(src))
      real*8 rf(dimdec(f))
      real*8 fac
      integer i, j, k
      fac = 0.125d0
      do 10 k = fregl2, fregh2
         do 10 j = fregl1, fregh1
            do 10 i = fregl0, fregh0
 10            src(i,j,k) = src(i,j,k) + fac *
     @              (rf(i-1,j-1,k-1) + rf(i-1,j-1,k) +
     @               rf(i-1,j,k-1)   + rf(i-1,j,k) +
     @               rf(i,j-1,k-1)   + rf(i,j-1,k) +
     @               rf(i,j,k-1)     + rf(i,j,k))
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgfavg(src,
     @ dims(src),
     @ rc,
     @ dims(c),
     @ rf,
     @ dims(f),
     @ dims(creg),
     @ ir, jr, kr, idim, idir)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 rc(dimdec(c))
      real*8 rf(dimdec(f))
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
         endif
         fac0 = 0.25d0 * ir / (ir+1)
         do 10 k = cregl2, cregh2
            do 10 j = cregl1, cregh1
 10            src(i*ir,j*jr,k*kr) =
     @            src(i*ir,j*jr,k*kr) + fac0 *
     @              (rc(irc,j,k)   + rc(irc,j,k-1) +
     @               rc(irc,j-1,k) + rc(irc,j-1,k-1))
         fac0 = fac0 / (ir * jr * kr * jr * kr)
         i = i * ir
         do 20 l = 0, kr-1
            fac1 = (kr-l) * fac0
            if (l .eq. 0) fac1 = 0.5d0 * fac1
            do 20 n = 0, jr-1
               fac = (jr-n) * fac1
               if (n .eq. 0) fac = 0.5d0 * fac
               do 20 k = kr*cregl2, kr*cregh2, kr
                  do 20 j = jr*cregl1, jr*cregh1, jr
 20                  src(i,j,k) = src(i,j,k) + fac *
     @                 (rf(irf,j-n,k-l) + rf(irf,j-n,k-l-1) +
     @                  rf(irf,j-n-1,k-l) + rf(irf,j-n-1,k-l-1) +
     @                  rf(irf,j-n,k+l) + rf(irf,j-n,k+l-1) +
     @                  rf(irf,j-n-1,k+l) + rf(irf,j-n-1,k+l-1) +
     @                  rf(irf,j+n,k-l) + rf(irf,j+n,k-l-1) +
     @                  rf(irf,j+n-1,k-l) + rf(irf,j+n-1,k-l-1) +
     @                  rf(irf,j+n,k+l) + rf(irf,j+n,k+l-1) +
     @                  rf(irf,j+n-1,k+l) + rf(irf,j+n-1,k+l-1))
      else if (idim .eq. 1) then
         j = cregl1
         if (idir .eq. 1) then
            jrc = j - 1
            jrf = j * jr
         else
            jrc = j
            jrf = j * jr - 1
         endif
         fac0 = 0.25d0 * jr / (jr+1)
         do 30 k = cregl2, cregh2
            do 30 i = cregl0, cregh0
 30            src(i*ir,j*jr,k*kr) =
     @            src(i*ir,j*jr,k*kr) + fac0 *
     @           (rc(i,jrc,k)   + rc(i-1,jrc,k) +
     @            rc(i,jrc,k-1) + rc(i-1,jrc,k-1))
         fac0 = fac0 / (ir * jr * kr * ir * kr)
         j = j * jr
         do 40 l = 0, kr-1
            fac1 = (kr-l) * fac0
            if (l .eq. 0) fac1 = 0.5d0 * fac1
            do 40 m = 0, ir-1
               fac = (ir-m) * fac1
               if (m .eq. 0) fac = 0.5d0 * fac
               do 40 k = kr*cregl2, kr*cregh2, kr
                  do 40 i = ir*cregl0, ir*cregh0, ir
 40                  src(i,j,k) = src(i,j,k) + fac *
     @                 (rf(i-m,jrf,k-l) + rf(i-m-1,jrf,k-l) +
     @                  rf(i-m,jrf,k-l-1) + rf(i-m-1,jrf,k-l-1) +
     @                  rf(i-m,jrf,k+l) + rf(i-m-1,jrf,k+l) +
     @                  rf(i-m,jrf,k+l-1) + rf(i-m-1,jrf,k+l-1) +
     @                  rf(i+m,jrf,k-l) + rf(i+m-1,jrf,k-l) +
     @                  rf(i+m,jrf,k-l-1) + rf(i+m-1,jrf,k-l-1) +
     @                  rf(i+m,jrf,k+l) + rf(i+m-1,jrf,k+l) +
     @                  rf(i+m,jrf,k+l-1) + rf(i+m-1,jrf,k+l-1))
      else
         k = cregl2
         if (idir .eq. 1) then
            krc = k - 1
            krf = k * kr
         else
            krc = k
            krf = k * kr - 1
         endif
         fac0 = 0.25d0 * kr / (kr+1)
         do 50 j = cregl1, cregh1
            do 50 i = cregl0, cregh0
 50            src(i*ir,j*jr,k*kr) =
     @            src(i*ir,j*jr,k*kr) + fac0 *
     @           (rc(i,j,krc)   + rc(i-1,j,krc) +
     @            rc(i,j-1,krc) + rc(i-1,j-1,krc))
         fac0 = fac0 / (ir * jr * kr * ir * jr)
         k = k * kr
         do 60 n = 0, jr-1
            fac1 = (jr-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do 60 m = 0, ir-1
               fac = (ir-m) * fac1
               if (m .eq. 0) fac = 0.5d0 * fac
               do 60 j = jr*cregl1, jr*cregh1, jr
                  do 60 i = ir*cregl0, ir*cregh0, ir
 60                  src(i,j,k) = src(i,j,k) + fac *
     @                 (rf(i-m,j-n,krf) + rf(i-m-1,j-n,krf) +
     @                  rf(i-m,j-n-1,krf) + rf(i-m-1,j-n-1,krf) +
     @                  rf(i-m,j+n,krf) + rf(i-m-1,j+n,krf) +
     @                  rf(i-m,j+n-1,krf) + rf(i-m-1,j+n-1,krf) +
     @                  rf(i+m,j-n,krf) + rf(i+m-1,j-n,krf) +
     @                  rf(i+m,j-n-1,krf) + rf(i+m-1,j-n-1,krf) +
     @                  rf(i+m,j+n,krf) + rf(i+m-1,j+n,krf) +
     @                  rf(i+m,j+n-1,krf) + rf(i+m-1,j+n-1,krf))
      endif
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgeavg(src,
     @ dims(src),
     @ rc,
     @ dims(c),
     @ rf,
     @ dims(f),
     @ dims(creg),
     @ ir, jr, kr, ivect, ga)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 rc(dimdec(c))
      real*8 rf(dimdec(f))
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
         do 10 ki = 0, 1
            do 10 ji = 0, 1
               if (ga(0,ji,ki) .eq. 1) then
                  center = center + ffac
               else
                  center = center + cfac
               endif
 10            continue
         do 20 ki = 0, 1
            do 20 ji = 0, 1
               if (ga(0,ji,ki) - ga(0,ji,1-ki) .eq. 1) then
                  ffac = 2.d0 * ir * (jr - 1) * rm3
                  center = center + ffac
               endif
               if (ga(0,ji,ki) - ga(0,1-ji,ki) .eq. 1) then
                  ffac = 2.d0 * ir * (kr - 1) * rm3
                  center = center + ffac
               endif
 20            continue
         center = 1.d0 / center
         fac1 = center * rm3 / ir
         do 50 ki = 0, 1
            do 50 ji = 0, 1
               if (ga(0,ji,ki) .eq. 1) then
                  krf = kf + ki - 1
                  jrf = jf + ji - 1
                  do 30 m = 0, ir-1
                     fac = (ir-m) * fac1
                     if (m .eq. 0) fac = 0.5d0 * fac
                     do 30 if = ir*cregl0, ir*cregh0, ir
 30                     src(if,jf,kf) = src(if,jf,kf) + fac *
     @                    (rf(if-m,jrf,krf) + rf(if-m-1,jrf,krf) +
     @                     rf(if+m,jrf,krf) + rf(if+m-1,jrf,krf))
               else
                  krc = kc + ki - 1
                  jrc = jc + ji - 1
                  do 40 ic = cregl0, cregh0
                     if = ic * ir
 40                  src(if,jf,kf) = src(if,jf,kf) + center *
     @                 (rc(ic,jrc,krc) + rc(ic-1,jrc,krc))
               endif
 50            continue
c faces
c each face is two faces and two sides of an edge
         do 80 ki = 0, 1
            kdir = 2 * ki - 1
            do 80 ji = 0, 1
               jdir = 2 * ji - 1
               if (ga(0,ji,ki) - ga(0,ji,1-ki) .eq. 1) then
                  krf = kf + ki - 1
                  fac0 = center * rm3 / (ir * jr)
                  do 60 n = jdir, jdir*(jr-1), jdir
                     fac1 = (jr-abs(n)) * fac0
                     do 60 m = 0, ir-1
                        fac = (ir-m) * fac1
                        if (m .eq. 0) fac = 0.5d0 * fac
                        do 60 if = ir*cregl0, ir*cregh0, ir
 60                        src(if,jf,kf) = src(if,jf,kf) + fac *
     @                       (rf(if-m,jf+n,krf) + rf(if-m-1,jf+n,krf) +
     @                        rf(if-m,jf+n-1,krf) + rf(if-m-1,jf+n-1,krf) +
     @                        rf(if+m,jf+n,krf) + rf(if+m-1,jf+n,krf) +
     @                        rf(if+m,jf+n-1,krf) + rf(if+m-1,jf+n-1,krf))
               endif
               if (ga(0,ji,ki) - ga(0,1-ji,ki) .eq. 1) then
                  jrf = jf + ji - 1
                  fac0 = center * rm3 / (ir * kr)
                  do 70 l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do 70 m = 0, ir-1
                        fac = (ir-m) * fac1
                        if (m .eq. 0) fac = 0.5d0 * fac
                        do 70 if = ir*cregl0, ir*cregh0, ir
 70                        src(if,jf,kf) = src(if,jf,kf) + fac *
     @                       (rf(if-m,jrf,kf+l) + rf(if-m-1,jrf,kf+l) +
     @                        rf(if-m,jrf,kf+l-1) + rf(if-m-1,jrf,kf+l-1) +
     @                        rf(if+m,jrf,kf+l) + rf(if+m-1,jrf,kf+l) +
     @                        rf(if+m,jrf,kf+l-1) + rf(if+m-1,jrf,kf+l-1))
               endif
 80            continue
      else if (ivect(1) .eq. 0) then
c quadrants
c each quadrant is two octants and their share of the two central edges
         ffac = 2.d0 * jr * rm3
         cfac = 2.d0
         do 110 ki = 0, 1
            do 110 ii = 0, 1
               if (ga(ii,0,ki) .eq. 1) then
                  center = center + ffac
               else
                  center = center + cfac
               endif
 110           continue
         do 120 ki = 0, 1
            do 120 ii = 0, 1
               if (ga(ii,0,ki) - ga(ii,0,1-ki) .eq. 1) then
                  ffac = 2.d0 * jr * (ir - 1) * rm3
                  center = center + ffac
               endif
               if (ga(ii,0,ki) - ga(1-ii,0,ki) .eq. 1) then
                  ffac = 2.d0 * jr * (kr - 1) * rm3
                  center = center + ffac
               endif
 120           continue
         center = 1.d0 / center
         fac1 = center * rm3 / jr
         do 150 ki = 0, 1
            do 150 ii = 0, 1
               if (ga(ii,0,ki) .eq. 1) then
                  krf = kf + ki - 1
                  irf = if + ii - 1
                  do 130 n = 0, jr-1
                     fac = (jr-n) * fac1
                     if (n .eq. 0) fac = 0.5d0 * fac
                     do 130 jf = jr*cregl1, jr*cregh1, jr
 130                    src(if,jf,kf) = src(if,jf,kf) + fac *
     @                    (rf(irf,jf-n,krf) + rf(irf,jf-n-1,krf) +
     @                     rf(irf,jf+n,krf) + rf(irf,jf+n-1,krf))
               else
                  krc = kc + ki - 1
                  irc = ic + ii - 1
                  do 140 jc = cregl1, cregh1
                     jf = jc * jr
 140                 src(if,jf,kf) = src(if,jf,kf) + center *
     @                 (rc(irc,jc,krc) + rc(irc,jc-1,krc))
               endif
 150           continue
c faces
c each face is two faces and two sides of an edge
         do 180 ki = 0, 1
            kdir = 2 * ki - 1
            do 180 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,0,ki) - ga(ii,0,1-ki) .eq. 1) then
                  krf = kf + ki - 1
                  fac0 = center * rm3 / (ir * jr)
                  do 160 n = 0, jr-1
                     fac1 = (jr-n) * fac0
                     if (n .eq. 0) fac1 = 0.5d0 * fac1
                     do 160 m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        do 160 jf = jr*cregl1, jr*cregh1, jr
 160                       src(if,jf,kf) = src(if,jf,kf) + fac *
     @                       (rf(if+m,jf-n,krf) + rf(if+m-1,jf-n,krf) +
     @                        rf(if+m,jf-n-1,krf) + rf(if+m-1,jf-n-1,krf) +
     @                        rf(if+m,jf+n,krf) + rf(if+m-1,jf+n,krf) +
     @                        rf(if+m,jf+n-1,krf) + rf(if+m-1,jf+n-1,krf))
               endif
               if (ga(ii,0,ki) - ga(1-ii,0,ki) .eq. 1) then
                  irf = if + ii - 1
                  fac0 = center * rm3 / (jr * kr)
                  do 170 l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do 170 n = 0, jr-1
                        fac = (jr-n) * fac1
                        if (n .eq. 0) fac = 0.5d0 * fac
                        do 170 jf = jr*cregl1, jr*cregh1, jr
 170                       src(if,jf,kf) = src(if,jf,kf) + fac *
     @                       (rf(irf,jf-n,kf+l) + rf(irf,jf-n-1,kf+l) +
     @                        rf(irf,jf-n,kf+l-1) + rf(irf,jf-n-1,kf+l-1) +
     @                        rf(irf,jf+n,kf+l) + rf(irf,jf+n-1,kf+l) +
     @                        rf(irf,jf+n,kf+l-1) + rf(irf,jf+n-1,kf+l-1))
               endif
 180           continue
      else
c quadrants
c each quadrant is two octants and their share of the two central edges
         ffac = 2.d0 * kr * rm3
         cfac = 2.d0
         do 210 ji = 0, 1
            do 210 ii = 0, 1
               if (ga(ii,ji,0) .eq. 1) then
                  center = center + ffac
               else
                  center = center + cfac
               endif
 210           continue
         do 220 ji = 0, 1
            do 220 ii = 0, 1
               if (ga(ii,ji,0) - ga(ii,1-ji,0) .eq. 1) then
                  ffac = 2.d0 * kr * (ir - 1) * rm3
                  center = center + ffac
               endif
               if (ga(ii,ji,0) - ga(1-ii,ji,0) .eq. 1) then
                  ffac = 2.d0 * kr * (jr - 1) * rm3
                  center = center + ffac
               endif
 220           continue
         center = 1.d0 / center
         fac1 = center * rm3 / kr
         do 250 ji = 0, 1
            do 250 ii = 0, 1
               if (ga(ii,ji,0) .eq. 1) then
                  jrf = jf + ji - 1
                  irf = if + ii - 1
                  do 230 l = 0, kr-1
                     fac = (kr-l) * fac1
                     if (l .eq. 0) fac = 0.5d0 * fac
                     do 230 kf = kr*cregl2, kr*cregh2, kr
 230                    src(if,jf,kf) = src(if,jf,kf) + fac *
     @                    (rf(irf,jrf,kf-l) + rf(irf,jrf,kf-l-1) +
     @                     rf(irf,jrf,kf+l) + rf(irf,jrf,kf+l-1))
               else
                  jrc = jc + ji - 1
                  irc = ic + ii - 1
                  do 240 kc = cregl2, cregh2
                     kf = kc * kr
 240                 src(if,jf,kf) = src(if,jf,kf) + center *
     @                 (rc(irc,jrc,kc) + rc(irc,jrc,kc-1))
               endif
 250           continue
c faces
c each face is two faces and two sides of an edge
         do 280 ji = 0, 1
            jdir = 2 * ji - 1
            do 280 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,0) - ga(ii,1-ji,0) .eq. 1) then
                  jrf = jf + ji - 1
                  fac0 = center * rm3 / (ir * kr)
                  do 260 l = 0, kr-1
                     fac1 = (kr-l) * fac0
                     if (l .eq. 0) fac1 = 0.5d0 * fac1
                     do 260 m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
                        do 260 kf = kr*cregl2, kr*cregh2, kr
 260                       src(if,jf,kf) = src(if,jf,kf) + fac *
     @                       (rf(if+m,jrf,kf-l) + rf(if+m-1,jrf,kf-l) +
     @                        rf(if+m,jrf,kf-l-1) + rf(if+m-1,jrf,kf-l-1) +
     @                        rf(if+m,jrf,kf+l) + rf(if+m-1,jrf,kf+l) +
     @                        rf(if+m,jrf,kf+l-1) + rf(if+m-1,jrf,kf+l-1))
               endif
               if (ga(ii,ji,0) - ga(1-ii,ji,0) .eq. 1) then
                  irf = if + ii - 1
                  fac0 = center * rm3 / (jr * kr)
                  do 270 l = 0, kr-1
                     fac1 = (kr-l) * fac0
                     if (l .eq. 0) fac1 = 0.5d0 * fac1
                     do 270 n = jdir, jdir*(jr-1), jdir
                        fac = (jr-abs(n)) * fac1
                        do 270 kf = kr*cregl2, kr*cregh2, kr
 270                       src(if,jf,kf) = src(if,jf,kf) + fac *
     @                       (rf(irf,jf+n,kf-l) + rf(irf,jf+n-1,kf-l) +
     @                        rf(irf,jf+n,kf-l-1) + rf(irf,jf+n-1,kf-l-1) +
     @                        rf(irf,jf+n,kf+l) + rf(irf,jf+n-1,kf+l) +
     @                        rf(irf,jf+n,kf+l-1) + rf(irf,jf+n-1,kf+l-1))
               endif
 280           continue
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
     @ ir, jr, kr, ga)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 rc(dimdec(c))
      real*8 rf(dimdec(f))
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
      do 10 ki = 0, 1
         krf = kf + ki - 1
         krc = kc + ki - 1
         do 10 ji = 0, 1
            jrf = jf + ji - 1
            jrc = jc + ji - 1
            do 10 ii = 0, 1
               if (ga(ii,ji,ki) .eq. 1) then
                  irf = if + ii - 1
                  center = center + ffac
                  sum = sum + fac * rf(irf,jrf,krf)
               else
                  irc = ic + ii - 1
                  center = center + cfac
                  sum = sum + rc(irc,jrc,krc)
               endif
 10            continue
c faces
      do 50 ki = 0, 1
         kdir = 2 * ki - 1
         krf = kf + ki - 1
         do 50 ji = 0, 1
            jdir = 2 * ji - 1
            jrf = jf + ji - 1
            do 50 ii = 0, 1
               idir = 2 * ii - 1
               irf = if + ii - 1
               if (ga(ii,ji,ki) - ga(ii,ji,1-ki) .eq. 1) then
                  fac0 = rm3 / (ir * jr)
                  ffac = (ir-1) * (jr-1) * rm3
                  center = center + ffac
                  do 20 n = jdir, jdir*(jr-1), jdir
                     fac1 = (jr-abs(n)) * fac0
                     do 20 m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
 20                     sum = sum + fac *
     @                    (rf(if+m,jf+n,krf) + rf(if+m-1,jf+n,krf) +
     @                     rf(if+m,jf+n-1,krf) + rf(if+m-1,jf+n-1,krf))
               endif
               if (ga(ii,ji,ki) - ga(ii,1-ji,ki) .eq. 1) then
                  fac0 = rm3 / (ir * kr)
                  ffac = (ir-1) * (kr-1) * rm3
                  center = center + ffac
                  do 30 l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do 30 m = idir, idir*(ir-1), idir
                        fac = (ir-abs(m)) * fac1
 30                     sum = sum + fac *
     @                    (rf(if+m,jrf,kf+l) + rf(if+m-1,jrf,kf+l) +
     @                     rf(if+m,jrf,kf+l-1) + rf(if+m-1,jrf,kf+l-1))
               endif
               if (ga(ii,ji,ki) - ga(1-ii,ji,ki) .eq. 1) then
                  fac0 = rm3 / (jr * kr)
                  ffac = (jr-1) * (kr-1) * rm3
                  center = center + ffac
                  do 40 l = kdir, kdir*(kr-1), kdir
                     fac1 = (kr-abs(l)) * fac0
                     do 40 n = jdir, jdir*(jr-1), jdir
                        fac = (jr-abs(n)) * fac1
 40                     sum = sum + fac *
     @                    (rf(irf,jf+n,kf+l) + rf(irf,jf+n,kf+l-1) +
     @                     rf(irf,jf+n-1,kf+l) + rf(irf,jf+n-1,kf+l-1))
               endif
 50            continue
c edges
      do 90 ki = 0, 1
         kdir = 2 * ki - 1
         krf = kf + ki - 1
         do 90 ji = 0, 1
            jdir = 2 * ji - 1
            jrf = jf + ji - 1
            do 90 ii = 0, 1
               idir = 2 * ii - 1
               irf = if + ii - 1
               if (ga(ii,ji,ki) -
     @             min(ga(ii,ji,1-ki), ga(ii,1-ji,ki), ga(ii,1-ji,1-ki))
     @             .eq. 1) then
                  fac1 = rm3 / ir
                  ffac = (ir-1) * rm3
                  center = center + ffac
                  do 60 m = idir, idir*(ir-1), idir
                     fac = (ir-abs(m)) * fac1
 60                  sum = sum + fac *
     @                 (rf(if+m,jrf,krf) + rf(if+m-1,jrf,krf))
               endif
               if (ga(ii,ji,ki) -
     @             min(ga(ii,ji,1-ki), ga(1-ii,ji,ki), ga(1-ii,ji,1-ki))
     @             .eq. 1) then
                  fac1 = rm3 / jr
                  ffac = (jr-1) * rm3
                  center = center + ffac
                  do 70 n = jdir, jdir*(jr-1), jdir
                     fac = (jr-abs(n)) * fac1
 70                  sum = sum + fac *
     @                 (rf(irf,jf+n,krf) + rf(irf,jf+n-1,krf))
               endif
               if (ga(ii,ji,ki) -
     @             min(ga(ii,1-ji,ki), ga(1-ii,ji,ki), ga(1-ii,1-ji,ki))
     @             .eq. 1) then
                  fac1 = rm3 / kr
                  ffac = (kr-1) * rm3
                  center = center + ffac
                  do 80 l = kdir, kdir*(kr-1), kdir
                     fac = (kr-abs(l)) * fac1
 80                  sum = sum + fac *
     @                 (rf(irf,jrf,kf+l) + rf(irf,jrf,kf+l-1))
               endif
 90            continue
c weighting
      src(if,jf,kf) = src(if,jf,kf) + sum / center
      end
