
#include "hg_version.H"

#define dims(a) a l0,a h0,a l1,a h1,a l2,a h2
#define dimdec(a) a l0:a h0,a l1:a h1,a l2:a h2
c#define dimdecu(a) (a h0-a l0+1)*(a h1-a l1+1)*(a h2-a l2+1)
#define dimdecu(a) *

      subroutine gvgrad(gp,
     @ dims(gp),
     @ pot,
     @ dims(pot),
     @ dims(reg),
     @ h)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(gp)
      integer dims(pot)
      integer dims(reg)
      real*8 gp(dimdec(gp), 3)
      real*8 pot(dimdec(pot))
      real*8 h(3)
      hxm1h = 0.25d0 / h(1)
      hym1h = 0.25d0 / h(2)
      hzm1h = 0.25d0 / h(3)
      do 10 k = regl2, regh2
         do 10 j = regl1, regh1
            do 10 i = regl0, regh0
               gp(i,j,k,1) = hxm1h * (pot(i+1,j,k) + pot(i+1,j+1,k) +
     @                                pot(i+1,j,k+1) + pot(i+1,j+1,k+1) -
     @                                pot(i,j,k) - pot(i,j+1,k) -
     @                                pot(i,j,k+1) - pot(i,j+1,k+1))
               gp(i,j,k,2) = hym1h * (pot(i,j+1,k) + pot(i+1,j+1,k) +
     @                                pot(i,j+1,k+1) + pot(i+1,j+1,k+1) -
     @                                pot(i,j,k) - pot(i+1,j,k) -
     @                                pot(i,j,k+1) - pot(i+1,j,k+1))
 10            gp(i,j,k,3) = hzm1h * (pot(i,j,k+1) + pot(i+1,j,k+1) +
     @                                pot(i,j+1,k+1) + pot(i+1,j+1,k+1) -
     @                                pot(i,j,k) - pot(i+1,j,k) -
     @                                pot(i,j+1,k) - pot(i+1,j+1,k))
      end

      subroutine gvavg(src,
     @ dims(src),
     @ rf,
     @ dims(f),
     @ dims(freg))
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(src)
      integer dims(f)
      integer dims(freg)
      real*8 src(dimdec(src))
      real*8 rf(dimdec(f))
      fac = 0.125d0
      do 10 k = fregl2, fregh2
         do 10 j = fregl1, fregh1
            do 10 i = fregl0, fregh0
 10            src(i,j,k) = fac *
     @              (rf(i-1,j-1,k-1) + rf(i-1,j-1,k) +
     @               rf(i-1,j,k-1)   + rf(i-1,j,k) +
     @               rf(i,j-1,k-1)   + rf(i,j-1,k) +
     @               rf(i,j,k-1)     + rf(i,j,k))
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine gvfavg(src,
     @ dims(src),
     @ rc,
     @ dims(c),
     @ rf,
     @ dims(f),
     @ dims(creg),
     @ irat, idim, idir)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 rc(dimdec(c))
      real*8 rf(dimdec(f))
      integer irat, idim, idir
      fac0 = 0.25d0 * irat / (irat+1)
      if (idim .eq. 0) then
         i = cregl0
         if (idir .eq. 1) then
            irc = i - 1
            irf = i * irat
         else
            irc = i
            irf = i * irat - 1
         endif
         do 10 k = cregl2, cregh2
            do 10 j = cregl1, cregh1
 10            src(i*irat,j*irat,k*irat) = fac0 *
     @            (rc(irc,j,k)   + rc(irc,j,k-1) +
     @             rc(irc,j-1,k) + rc(irc,j-1,k-1))
         fac0 = fac0 / (irat**5)
         i = i * irat
         do 20 l = 0, irat-1
            fac1 = (irat-l) * fac0
            if (l .eq. 0) fac1 = 0.5d0 * fac1
            do 20 n = 0, irat-1
               fac = (irat-n) * fac1
               if (n .eq. 0) fac = 0.5d0 * fac
               do 20 k = irat*cregl2, irat*cregh2, irat
                  do 20 j = irat*cregl1, irat*cregh1, irat
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
            jrf = j * irat
         else
            jrc = j
            jrf = j * irat - 1
         endif
         do 30 k = cregl2, cregh2
            do 30 i = cregl0, cregh0
 30            src(i*irat,j*irat,k*irat) = fac0 *
     @           (rc(i,jrc,k)   + rc(i-1,jrc,k) +
     @            rc(i,jrc,k-1) + rc(i-1,jrc,k-1))
         fac0 = fac0 / (irat**5)
         j = j * irat
         do 40 l = 0, irat-1
            fac1 = (irat-l) * fac0
            if (l .eq. 0) fac1 = 0.5d0 * fac1
            do 40 m = 0, irat-1
               fac = (irat-m) * fac1
               if (m .eq. 0) fac = 0.5d0 * fac
               do 40 k = irat*cregl2, irat*cregh2, irat
                  do 40 i = irat*cregl0, irat*cregh0, irat
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
            krf = k * irat
         else
            krc = k
            krf = k * irat - 1
         endif
         do 50 j = cregl1, cregh1
            do 50 i = cregl0, cregh0
 50            src(i*irat,j*irat,k*irat) = fac0 *
     @           (rc(i,j,krc)   + rc(i-1,j,krc) +
     @            rc(i,j-1,krc) + rc(i-1,j-1,krc))
         fac0 = fac0 / (irat**5)
         k = k * irat
         do 60 n = 0, irat-1
            fac1 = (irat-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do 60 m = 0, irat-1
               fac = (irat-m) * fac1
               if (m .eq. 0) fac = 0.5d0 * fac
               do 60 j = irat*cregl1, irat*cregh1, irat
                  do 60 i = irat*cregl0, irat*cregh0, irat
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
      subroutine gveavg(src,
     @ dims(src),
     @ rc,
     @ dims(c),
     @ rf,
     @ dims(f),
     @ dims(creg),
     @ irat, ivect, ga)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 rc(dimdec(c))
      real*8 rf(dimdec(f))
      integer irat, ivect(0:2), ga(0:1,0:1,0:1)
      rm1 = 1.d0 / irat
      rm2 = rm1 * rm1
      rm3 = rm1 * rm2
      ic = cregl0
      jc = cregl1
      kc = cregl2
      if = ic * irat
      jf = jc * irat
      kf = kc * irat
      center = 0.d0
      if (ivect(0) .eq. 0) then
         do 10 if = irat*cregl0, irat*cregh0, irat
 10         src(if,jf,kf) = 0.d0
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac1 = rm2 * rm2
         ffac = 2.d0 * rm2
         cfac = 2.d0
         do 40 ki = 0, 1
            kdir = 2 * ki - 1
            do 40 ji = 0, 1
               jdir = 2 * ji - 1
               if (ga(0,ji,ki) .eq. 1) then
                  krf = kf + ki - 1
                  jrf = jf + ji - 1
                  center = center + ffac
                  do 20 m = 0, irat-1
                     fac = (irat-m) * fac1
                     if (m .eq. 0) fac = 0.5d0 * fac
                     do 20 if = irat*cregl0, irat*cregh0, irat
 20                     src(if,jf,kf) = src(if,jf,kf) + fac *
     @                    (rf(if-m,jrf,krf) + rf(if-m-1,jrf,krf) +
     @                     rf(if+m,jrf,krf) + rf(if+m-1,jrf,krf))
               else
                  krc = kc + ki - 1
                  jrc = jc + ji - 1
                  center = center + cfac
                  do 30 ic = cregl0, cregh0
                     if = ic * irat
 30                  src(if,jf,kf) = src(if,jf,kf) +
     @                 (rc(ic,jrc,krc) + rc(ic-1,jrc,krc))
               endif
 40            continue
c faces
c each face is two faces and two sides of an edge
         fac0 = rm2 * rm3
         ffac = 2.d0 * (irat - 1) * rm2
         do 70 ki = 0, 1
            kdir = 2 * ki - 1
            do 70 ji = 0, 1
               jdir = 2 * ji - 1
               if (ga(0,ji,ki) - ga(0,ji,1-ki) .eq. 1) then
                  krf = kf + ki - 1
                  center = center + ffac
                  do 50 n = jdir, jdir*(irat-1), jdir
                     fac1 = (irat-abs(n)) * fac0
                     do 50 m = 0, irat-1
                        fac = (irat-m) * fac1
                        if (m .eq. 0) fac = 0.5d0 * fac
                        do 50 if = irat*cregl0, irat*cregh0, irat
 50                        src(if,jf,kf) = src(if,jf,kf) + fac *
     @                       (rf(if-m,jf+n,krf) + rf(if-m-1,jf+n,krf) +
     @                        rf(if-m,jf+n-1,krf) + rf(if-m-1,jf+n-1,krf) +
     @                        rf(if+m,jf+n,krf) + rf(if+m-1,jf+n,krf) +
     @                        rf(if+m,jf+n-1,krf) + rf(if+m-1,jf+n-1,krf))
               endif
               if (ga(0,ji,ki) - ga(0,1-ji,ki) .eq. 1) then
                  jrf = jf + ji - 1
                  center = center + ffac
                  do 60 l = kdir, kdir*(irat-1), kdir
                     fac1 = (irat-abs(l)) * fac0
                     do 60 m = 0, irat-1
                        fac = (irat-m) * fac1
                        if (m .eq. 0) fac = 0.5d0 * fac
                        do 60 if = irat*cregl0, irat*cregh0, irat
 60                        src(if,jf,kf) = src(if,jf,kf) + fac *
     @                       (rf(if-m,jrf,kf+l) + rf(if-m-1,jrf,kf+l) +
     @                        rf(if-m,jrf,kf+l-1) + rf(if-m-1,jrf,kf+l-1) +
     @                        rf(if+m,jrf,kf+l) + rf(if+m-1,jrf,kf+l) +
     @                        rf(if+m,jrf,kf+l-1) + rf(if+m-1,jrf,kf+l-1))
               endif
 70            continue
c weighting
         do 80 if = irat*cregl0, irat*cregh0, irat
 80         src(if,jf,kf) = src(if,jf,kf) / center
      else if (ivect(1) .eq. 0) then
         do 110 jf = irat*cregl1, irat*cregh1, irat
 110        src(if,jf,kf) = 0.d0
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac1 = rm2 * rm2
         ffac = 2.d0 * rm2
         cfac = 2.d0
         do 140 ki = 0, 1
            kdir = 2 * ki - 1
            do 140 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,0,ki) .eq. 1) then
                  krf = kf + ki - 1
                  irf = if + ii - 1
                  center = center + ffac
                  do 120 n = 0, irat-1
                     fac = (irat-n) * fac1
                     if (n .eq. 0) fac = 0.5d0 * fac
                     do 120 jf = irat*cregl1, irat*cregh1, irat
 120                    src(if,jf,kf) = src(if,jf,kf) + fac *
     @                    (rf(irf,jf-n,krf) + rf(irf,jf-n-1,krf) +
     @                     rf(irf,jf+n,krf) + rf(irf,jf+n-1,krf))
               else
                  krc = kc + ki - 1
                  irc = ic + ii - 1
                  center = center + cfac
                  do 130 jc = cregl1, cregh1
                     jf = jc * irat
 130                 src(if,jf,kf) = src(if,jf,kf) +
     @                 (rc(irc,jc,krc) + rc(irc,jc-1,krc))
               endif
 140           continue
c faces
c each face is two faces and two sides of an edge
         fac0 = rm2 * rm3
         ffac = 2.d0 * (irat - 1) * rm2
         do 170 ki = 0, 1
            kdir = 2 * ki - 1
            do 170 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,0,ki) - ga(ii,0,1-ki) .eq. 1) then
                  krf = kf + ki - 1
                  center = center + ffac
                  do 150 n = 0, irat-1
                     fac1 = (irat-n) * fac0
                     if (n .eq. 0) fac1 = 0.5d0 * fac1
                     do 150 m = idir, idir*(irat-1), idir
                        fac = (irat-abs(m)) * fac1
                        do 150 jf = irat*cregl1, irat*cregh1, irat
 150                       src(if,jf,kf) = src(if,jf,kf) + fac *
     @                       (rf(if+m,jf-n,krf) + rf(if+m-1,jf-n,krf) +
     @                        rf(if+m,jf-n-1,krf) + rf(if+m-1,jf-n-1,krf) +
     @                        rf(if+m,jf+n,krf) + rf(if+m-1,jf+n,krf) +
     @                        rf(if+m,jf+n-1,krf) + rf(if+m-1,jf+n-1,krf))
               endif
               if (ga(ii,0,ki) - ga(1-ii,0,ki) .eq. 1) then
                  irf = if + ii - 1
                  center = center + ffac
                  do 160 l = kdir, kdir*(irat-1), kdir
                     fac1 = (irat-abs(l)) * fac0
                     do 160 n = 0, irat-1
                        fac = (irat-n) * fac1
                        if (n .eq. 0) fac = 0.5d0 * fac
                        do 160 jf = irat*cregl1, irat*cregh1, irat
 160                       src(if,jf,kf) = src(if,jf,kf) + fac *
     @                       (rf(irf,jf-n,kf+l) + rf(irf,jf-n-1,kf+l) +
     @                        rf(irf,jf-n,kf+l-1) + rf(irf,jf-n-1,kf+l-1) +
     @                        rf(irf,jf+n,kf+l) + rf(irf,jf+n-1,kf+l) +
     @                        rf(irf,jf+n,kf+l-1) + rf(irf,jf+n-1,kf+l-1))
               endif
 170           continue
c weighting
         do 180 jf = irat*cregl1, irat*cregh1, irat
 180        src(if,jf,kf) = src(if,jf,kf) / center
      else
         do 210 kf = irat*cregl2, irat*cregh2, irat
 210        src(if,jf,kf) = 0.d0
c quadrants
c each quadrant is two octants and their share of the two central edges
         fac1 = rm2 * rm2
         ffac = 2.d0 * rm2
         cfac = 2.d0
         do 240 ji = 0, 1
            jdir = 2 * ji - 1
            do 240 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,0) .eq. 1) then
                  jrf = jf + ji - 1
                  irf = if + ii - 1
                  center = center + ffac
                  do 220 l = 0, irat-1
                     fac = (irat-l) * fac1
                     if (l .eq. 0) fac = 0.5d0 * fac
                     do 220 kf = irat*cregl2, irat*cregh2, irat
 220                    src(if,jf,kf) = src(if,jf,kf) + fac *
     @                    (rf(irf,jrf,kf-l) + rf(irf,jrf,kf-l-1) +
     @                     rf(irf,jrf,kf+l) + rf(irf,jrf,kf+l-1))
               else
                  jrc = jc + ji - 1
                  irc = ic + ii - 1
                  center = center + cfac
                  do 230 kc = cregl2, cregh2
                     kf = kc * irat
 230                 src(if,jf,kf) = src(if,jf,kf) +
     @                 (rc(irc,jrc,kc) + rc(irc,jrc,kc-1))
               endif
 240           continue
c faces
c each face is two faces and two sides of an edge
         fac0 = rm2 * rm3
         ffac = 2.d0 * (irat - 1) * rm2
         do 270 ji = 0, 1
            jdir = 2 * ji - 1
            do 270 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,0) - ga(ii,1-ji,0) .eq. 1) then
                  jrf = jf + ji - 1
                  center = center + ffac
                  do 250 l = 0, irat-1
                     fac1 = (irat-l) * fac0
                     if (l .eq. 0) fac1 = 0.5d0 * fac1
                     do 250 m = idir, idir*(irat-1), idir
                        fac = (irat-abs(m)) * fac1
                        do 250 kf = irat*cregl2, irat*cregh2, irat
 250                       src(if,jf,kf) = src(if,jf,kf) + fac *
     @                       (rf(if+m,jrf,kf-l) + rf(if+m-1,jrf,kf-l) +
     @                        rf(if+m,jrf,kf-l-1) + rf(if+m-1,jrf,kf-l-1) +
     @                        rf(if+m,jrf,kf+l) + rf(if+m-1,jrf,kf+l) +
     @                        rf(if+m,jrf,kf+l-1) + rf(if+m-1,jrf,kf+l-1))
               endif
               if (ga(ii,ji,0) - ga(1-ii,ji,0) .eq. 1) then
                  irf = if + ii - 1
                  center = center + ffac
                  do 260 l = 0, irat-1
                     fac1 = (irat-l) * fac0
                     if (l .eq. 0) fac1 = 0.5d0 * fac1
                     do 260 n = jdir, jdir*(irat-1), jdir
                        fac = (irat-abs(n)) * fac1
                        do 260 kf = irat*cregl2, irat*cregh2, irat
 260                       src(if,jf,kf) = src(if,jf,kf) + fac *
     @                       (rf(irf,jf+n,kf-l) + rf(irf,jf+n-1,kf-l) +
     @                        rf(irf,jf+n,kf-l-1) + rf(irf,jf+n-1,kf-l-1) +
     @                        rf(irf,jf+n,kf+l) + rf(irf,jf+n-1,kf+l) +
     @                        rf(irf,jf+n,kf+l-1) + rf(irf,jf+n-1,kf+l-1))
               endif
 270           continue
c weighting
         do 280 kf = irat*cregl2, irat*cregh2, irat
 280        src(if,jf,kf) = src(if,jf,kf) / center
      endif
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine gvcavg(src,
     @ dims(src),
     @ rc,
     @ dims(c),
     @ rf,
     @ dims(f),
     @ dims(creg),
     @ irat, ga)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 rc(dimdec(c))
      real*8 rf(dimdec(f))
      integer irat, ga(0:1,0:1,0:1)
      rm1 = 1.d0 / irat
      rm2 = rm1 * rm1
      rm3 = rm1 * rm2
      ic = cregl0
      jc = cregl1
      kc = cregl2
      if = ic * irat
      jf = jc * irat
      kf = kc * irat
      sum = 0.d0
      center = 0.d0
c octants
      fac = rm3
      ffac = rm3
      cfac = 1.d0
      do 10 ki = 0, 1
         kdir = 2 * ki - 1
         krf = kf + ki - 1
         krc = kc + ki - 1
         do 10 ji = 0, 1
            jdir = 2 * ji - 1
            jrf = jf + ji - 1
            jrc = jc + ji - 1
            do 10 ii = 0, 1
               idir = 2 * ii - 1
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
      fac0 = rm2 * rm3
      ffac = (irat-1) * (irat-1) * rm3
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
                  center = center + ffac
                  do 20 n = jdir, jdir*(irat-1), jdir
                     fac1 = (irat-abs(n)) * fac0
                     do 20 m = idir, idir*(irat-1), idir
                        fac = (irat-abs(m)) * fac1
 20                     sum = sum + fac *
     @                    (rf(if+m,jf+n,krf) + rf(if+m-1,jf+n,krf) +
     @                     rf(if+m,jf+n-1,krf) + rf(if+m-1,jf+n-1,krf))
               endif
               if (ga(ii,ji,ki) - ga(ii,1-ji,ki) .eq. 1) then
                  center = center + ffac
                  do 30 l = kdir, kdir*(irat-1), kdir
                     fac1 = (irat-abs(l)) * fac0
                     do 30 m = idir, idir*(irat-1), idir
                        fac = (irat-abs(m)) * fac1
 30                     sum = sum + fac *
     @                    (rf(if+m,jrf,kf+l) + rf(if+m-1,jrf,kf+l) +
     @                     rf(if+m,jrf,kf+l-1) + rf(if+m-1,jrf,kf+l-1))
               endif
               if (ga(ii,ji,ki) - ga(1-ii,ji,ki) .eq. 1) then
                  center = center + ffac
                  do 40 l = kdir, kdir*(irat-1), kdir
                     fac1 = (irat-abs(l)) * fac0
                     do 40 n = jdir, jdir*(irat-1), jdir
                        fac = (irat-abs(n)) * fac1
 40                     sum = sum + fac *
     @                    (rf(irf,jf+n,kf+l) + rf(irf,jf+n,kf+l-1) +
     @                     rf(irf,jf+n-1,kf+l) + rf(irf,jf+n-1,kf+l-1))
               endif
 50            continue
c edges
      fac1 = rm2 * rm2
      ffac = (irat-1) * rm3
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
                  center = center + ffac
                  do 60 m = idir, idir*(irat-1), idir
                     fac = (irat-abs(m)) * fac1
 60                  sum = sum + fac *
     @                 (rf(if+m,jrf,krf) + rf(if+m-1,jrf,krf))
               endif
               if (ga(ii,ji,ki) -
     @             min(ga(ii,ji,1-ki), ga(1-ii,ji,ki), ga(1-ii,ji,1-ki))
     @             .eq. 1) then
                  center = center + ffac
                  do 70 n = jdir, jdir*(irat-1), jdir
                     fac = (irat-abs(n)) * fac1
 70                  sum = sum + fac *
     @                 (rf(irf,jf+n,krf) + rf(irf,jf+n-1,krf))
               endif
               if (ga(ii,ji,ki) -
     @             min(ga(ii,1-ji,ki), ga(1-ii,ji,ki), ga(1-ii,1-ji,ki))
     @             .eq. 1) then
                  center = center + ffac
                  do 80 l = kdir, kdir*(irat-1), kdir
                     fac = (irat-abs(l)) * fac1
 80                  sum = sum + fac *
     @                 (rf(irf,jrf,kf+l) + rf(irf,jrf,kf+l-1))
               endif
 90            continue
c weighting
      src(if,jf,kf) = sum / center
      end
