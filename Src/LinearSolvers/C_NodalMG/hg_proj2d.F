
#define dims(a) a##l0, a##h0, a##l1, a##h1
#define dimdec(a) a##l0:a##h0,a##l1:a##h1

#ifdef HG_TERRAIN

      subroutine hggrad(gpx, gpy, dims(gp),
     @ dest, dims(dest),
     @ dims(freg))
      integer dims(gp)
      integer dims(dest)
      integer dims(freg)
      real*8 gpx(dimdec(gp))
      real*8 gpy(dimdec(gp))
      real*8 dest(dimdec(dest))
      integer i, j
      do 10 j = fregl1, fregh1
         do 10 i = fregl0, fregh0
            gpx(i,j) = 0.5d0 * (dest(i+1,j) + dest(i+1,j+1) -
     @                          dest(i  ,j) - dest(i  ,j+1))
 10         gpy(i,j) = 0.5d0 * (dest(i,j+1) + dest(i+1,j+1) -
     @                          dest(i,j  ) - dest(i+1,j  ))
      end

      subroutine hgdiv(src, dims(src),
     @ uf, vf, dims(f),
     @ dims(freg))
      integer dims(src)
      integer dims(f)
      integer dims(freg)
      real*8 src(dimdec(src))
      real*8 uf(dimdec(f))
      real*8 vf(dimdec(f))
      real*8 fac
      integer i, j
      fac = 0.5d0
      do 10 j = fregl1, fregh1
         do 10 i = fregl0, fregh0
 10         src(i,j) = fac *
     @        (uf(i,j-1) - uf(i-1,j-1) +
     @         uf(i,j  ) - uf(i-1,j) +
     @         vf(i-1,j) - vf(i-1,j-1) +
     @         vf(i  ,j) - vf(i,j-1))
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgfdiv(src, dims(src),
     @ uc, vc, dims(c),
     @ uf, vf, dims(f),
     @ dims(creg),
     @ ir, jr, idim, idir)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 uc(dimdec(c))
      real*8 vc(dimdec(c))
      real*8 uf(dimdec(f))
      real*8 vf(dimdec(f))
      integer ir, jr, idim, idir
      real*8 fac0, fac1
      integer i, j, iuf, juf, iuc, juc, m, n
      if (idim .eq. 0) then
         i = cregl0
         if (idir .eq. 1) then
            iuc = i - 1
            iuf = i * ir
         else
            iuc = i
            iuf = i * ir - 1
         end if
         fac0 = 0.5d0
         do 10 j = cregl1, cregh1
 10         src(i*ir,j*jr) = fac0 *
     @        ((vc(iuc,j) - vc(iuc,j-1)) -
     @         idir * (uc(iuc,j) + uc(iuc,j-1)))
         fac0 = fac0 / jr
         i = i * ir
         do 30 n = 0, jr-1
            fac1 = (jr-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do 30 j = jr*cregl1, jr*cregh1, jr
 30            src(i,j) = src(i,j) + fac1 *
     @           (idir * (uf(iuf,j-n) + uf(iuf,j-n-1) +
     @                    uf(iuf,j+n) + uf(iuf,j+n-1)) +
     @                   (vf(iuf,j-n) - vf(iuf,j-n-1) +
     @                    vf(iuf,j+n) - vf(iuf,j+n-1)))
      else
         j = cregl1
         if (idir .eq. 1) then
            juc = j - 1
            juf = j * jr
         else
            juc = j
            juf = j * jr - 1
         end if
         fac0 = 0.5d0
         do 60 i = cregl0, cregh0
 60         src(i*ir,j*jr) = fac0 *
     @        ((uc(i,juc) - uc(i-1,juc)) -
     @         idir * (vc(i,juc) + vc(i-1,juc)))
         fac0 = fac0 / ir
         j = j * jr
         do 80 m = 0, ir-1
            fac1 = (ir-m) * fac0
            if (m .eq. 0) fac1 = 0.5d0 * fac1
            do 80 i = ir*cregl0, ir*cregh0, ir
 80            src(i,j) = src(i,j) + fac1 *
     @           ((uf(i-m,juf) - uf(i-m-1,juf) +
     @             uf(i+m,juf) - uf(i+m-1,juf)) +
     @            idir * (vf(i-m,juf) + vf(i-m-1,juf) +
     @                    vf(i+m,juf) + vf(i+m-1,juf)))
      end if
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgcdiv(src, dims(src),
     @ uc, vc, dims(c),
     @ uf, vf, dims(f),
     @ dims(creg),
     @ ir, jr, ga)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 uc(dimdec(c))
      real*8 vc(dimdec(c))
      real*8 uf(dimdec(f))
      real*8 vf(dimdec(f))
      integer ir, jr, ga(0:1,0:1)
      real*8 sum, fac, fac1
      integer ic, jc, if, jf, ii, ji, idir, jdir, m, n
      ic = cregl0
      jc = cregl1
      if = ic * ir
      jf = jc * jr
      sum = 0.d0
c quadrants
      do 10 ji = 0, 1
         jdir = 2 * ji - 1
         do 10 ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,ji) .eq. 1) then
               sum = sum + idir * uf(if+ii-1,jf+ji-1)
     @                   + jdir * vf(if+ii-1,jf+ji-1)
            else
               sum = sum + idir * uc(ic+ii-1,jc+ji-1)
     @                   + jdir * vc(ic+ii-1,jc+ji-1)
            end if
 10         continue
c edges
      do 40 ji = 0, 1
         jdir = 2 * ji - 1
         do 40 ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,ji) - ga(ii,1-ji) .eq. 1) then
               fac1 = 1.d0 / ir
               do 20 m = idir, idir*(ir-1), idir
                  fac = (ir-abs(m)) * fac1
 20               sum = sum + fac * (       (uf(if+m,  jf+ji-1) -
     @                                       uf(if+m-1,jf+ji-1)) +
     @                               jdir * (vf(if+m,  jf+ji-1) +
     @                                       vf(if+m-1,jf+ji-1)))
            end if
            if (ga(ii,ji) - ga(1-ii,ji) .eq. 1) then
               fac1 = 1.d0 / jr
               do 30 n = jdir, jdir*(jr-1), jdir
                  fac = (jr-abs(n)) * fac1
 30               sum = sum + fac * (idir * (uf(if+ii-1,jf+n) +
     @                                       uf(if+ii-1,jf+n-1)) +
     @                                      (vf(if+ii-1,jf+n) -
     @                                       vf(if+ii-1,jf+n-1)))
            end if
 40         continue
c weighting
      src(if,jf) = 0.5d0 * sum
      end

#else

      subroutine hggrad(gpx, gpy, dims(gp),
     @ dest, dims(dest),
     @ dims(freg),
     @ hx, hy, irz)
      integer dims(gp)
      integer dims(dest)
      integer dims(freg)
      real*8 gpx(dimdec(gp))
      real*8 gpy(dimdec(gp))
      real*8 dest(dimdec(dest))
      real*8 hx, hy
      integer irz
      real*8 hxm1h, hym1h, fac, r
      integer i, j
      hxm1h = 0.5d0 / hx
      hym1h = 0.5d0 / hy
      do 10 j = fregl1, fregh1
         do 10 i = fregl0, fregh0
            gpx(i,j) = hxm1h * (dest(i+1,j) + dest(i+1,j+1) -
     @                          dest(i,j) - dest(i,j+1))
 10         gpy(i,j) = hym1h * (dest(i,j+1) + dest(i+1,j+1) -
     @                          dest(i,j) - dest(i+1,j))
      if (irz .eq. 1) then
         fac = hx / (12.d0 * hy)
         do 20 j = fregl1, fregh1
            r = (fregl0 - 0.5d0) * hx
            do 20 i = fregl0, fregh0
               r = r + hx
 20            gpy(i,j) = gpy(i,j) + (fac / r) *
     @                               (dest(i+1,j+1) - dest(i+1,j) -
     @                                dest(i,j+1) + dest(i,j))
      end if
      end

      subroutine hgdiv(src, dims(src),
     @ uf, vf, dims(f),
     @ dims(freg),
     @ hx, hy, irz, imax)
      integer dims(src)
      integer dims(f)
      integer dims(freg)
      real*8 src(dimdec(src))
      real*8 uf(dimdec(f))
      real*8 vf(dimdec(f))
      real*8 hx, hy
      integer irz, imax
      real*8 hxm1, hym1, fac, r1, r0m, r1m
      integer i, j
      hxm1 = 1.d0 / hx
      hym1 = 1.d0 / hy
      fac = 0.5d0
      do 10 j = fregl1, fregh1
         do 10 i = fregl0, fregh0
 10         src(i,j) = fac *
     @        (hxm1 * (uf(i,j-1) - uf(i-1,j-1) +
     @                 uf(i,j)   - uf(i-1,j)) +
     @         hym1 * (vf(i-1,j) - vf(i-1,j-1) +
     @                 vf(i,j)   - vf(i,j-1)))
      if (irz .eq. 1) then
         fac = hx / (12.d0 * hy)
         r1  = (fregl0 - 0.5d0) * hx
         r1m = 1.d0 / r1
         do 20 i = fregl0, fregh0
            r0m = r1m
            r1  = r1 + hx
            if (i .lt. imax) then
               r1m = 1.d0 / r1
            else
               r1m = -r0m
            end if
            do 20 j = fregl1, fregh1
 20            src(i,j) = src(i,j) + fac *
     @            (r0m * (vf(i-1,j) - vf(i-1,j-1)) +
     @             r1m * (vf(i,j-1) - vf(i,j)))
      end if
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgfdiv(src, dims(src),
     @ uc, vc, dims(c),
     @ uf, vf, dims(f),
     @ dims(creg),
     @ hx, hy, ir, jr, idim, idir, irz, imax)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 uc(dimdec(c))
      real*8 vc(dimdec(c))
      real*8 uf(dimdec(f))
      real*8 vf(dimdec(f))
      real*8 hx, hy
      integer ir, jr, idim, idir, irz, imax
      real*8 hxm1, hym1, fac0, fac1, r, rfac, rfac0, rfac1
      real*8 rfac0m, rfac1m, rfac0p, rfac1p
      integer i, j, iuf, juf, iuc, juc, m, n
      if (idim .eq. 0) then
         i = cregl0
         if (idir .eq. 1) then
            iuc = i - 1
            iuf = i * ir
         else
            iuc = i
            iuf = i * ir - 1
         end if
         fac0 = ir / (ir + 1.d0)
         hxm1 = 1.d0 / (ir * hx)
         hym1 = 1.d0 / (jr * hy)
         do 10 j = cregl1, cregh1
 10         src(i*ir,j*jr) = fac0 *
     @        (hym1 * (vc(iuc,j) - vc(iuc,j-1)) -
     @         hxm1 * idir * (uc(iuc,j) + uc(iuc,j-1)))
         if (irz .eq. 1) then
            r = (iuc + 0.5d0) * (hx * ir)
            rfac = (ir * hx) / (6.d0 * jr * hy * r)
            do 20 j = cregl1, cregh1
 20            src(i*ir,j*jr) = src(i*ir,j*jr) + fac0 *
     @           (rfac * idir * (vc(iuc,j) - vc(iuc,j-1)))
            r = (iuf + 0.5d0) * hx
            rfac = hx / (6.d0 * hy * r)
         end if
         fac0 = fac0 / (ir * jr * jr)
         hxm1 = ir * hxm1
         hym1 = jr * hym1
         i = i * ir
         do 50 n = 0, jr-1
            fac1 = (jr-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do 30 j = jr*cregl1, jr*cregh1, jr
 30            src(i,j) = src(i,j) + fac1 *
     @           (hxm1 * idir * (uf(iuf,j-n) + uf(iuf,j-n-1) +
     @                           uf(iuf,j+n) + uf(iuf,j+n-1)) +
     @            hym1 * (vf(iuf,j-n) - vf(iuf,j-n-1) +
     @                    vf(iuf,j+n) - vf(iuf,j+n-1)))
            if (irz .eq. 1) then
               do 40 j = jr*cregl1, jr*cregh1, jr
 40               src(i,j) = src(i,j) - fac1 *
     @              (rfac * idir * (vf(iuf,j-n) - vf(iuf,j-n-1) +
     @                              vf(iuf,j+n) - vf(iuf,j+n-1)))
            end if
 50         continue
      else
         j = cregl1
         if (idir .eq. 1) then
            juc = j - 1
            juf = j * jr
         else
            juc = j
            juf = j * jr - 1
         end if
         fac0 = jr / (jr + 1.d0)
         hxm1 = 1.d0 / (ir * hx)
         hym1 = 1.d0 / (jr * hy)
         do 60 i = cregl0, cregh0
 60         src(i*ir,j*jr) = fac0 *
     @        (hxm1 * (uc(i,juc) - uc(i-1,juc)) -
     @         hym1 * idir * (vc(i,juc) + vc(i-1,juc)))
         if (irz .eq. 1 .and. cregh0 .lt. imax) then
            do 70 i = cregl0, cregh0
               r = (i + 0.5d0) * (hx * ir)
               rfac0 = (ir * hx) / (6.d0 * jr * hy * (r - hx * ir))
               rfac1 = (ir * hx) / (6.d0 * jr * hy * r)
 70            src(i*ir,j*jr) = src(i*ir,j*jr) + fac0 *
     @           (idir * (rfac1 * vc(i,juc) - rfac0 * vc(i-1,juc)))
         else if (irz .eq. 1) then
c This should only occur with a corner at the outer boundary:
            i = cregh0
            r = (i - 0.5d0) * (hx * ir)
            rfac0 = (ir * hx) / (6.d0 * jr * hy * r)
            rfac1 = -rfac0
            src(i*ir,j*jr) = src(i*ir,j*jr) + fac0 *
     @           (idir * (rfac1 * vc(i,juc) - rfac0 * vc(i-1,juc)))
         end if
         fac0 = fac0 / (ir * ir * jr)
         hxm1 = ir * hxm1
         hym1 = jr * hym1
         j = j * jr
         do 100 m = 0, ir-1
            fac1 = (ir-m) * fac0
            if (m .eq. 0) fac1 = 0.5d0 * fac1
            do 80 i = ir*cregl0, ir*cregh0, ir
 80            src(i,j) = src(i,j) + fac1 *
     @           (hxm1 * (uf(i-m,juf) - uf(i-m-1,juf) +
     @                    uf(i+m,juf) - uf(i+m-1,juf)) +
     @            hym1 * idir * (vf(i-m,juf) + vf(i-m-1,juf) +
     @                           vf(i+m,juf) + vf(i+m-1,juf)))
            if (irz .eq. 1 .and. cregh0 .lt. imax) then
               do 90 i = ir*cregl0, ir*cregh0, ir
                  r = (i + 0.5d0) * hx
                  rfac0m = hx / (6.d0 * hy * (r - (m + 1) * hx))
                  rfac1m = hx / (6.d0 * hy * (r - m * hx))
                  rfac0p = hx / (6.d0 * hy * (r + (m - 1) * hx))
                  rfac1p = hx / (6.d0 * hy * (r + m * hx))
 90               src(i,j) = src(i,j) - fac1 *
     @              (idir * (rfac1m * vf(i-m,juf) - rfac0m * vf(i-m-1,juf) +
     @                       rfac1p * vf(i+m,juf) - rfac0p * vf(i+m-1,juf)))
            else if (irz .eq. 1) then
c This should only occur with a corner at the outer boundary:
               i = ir*cregh0
               r = (i + 0.5d0) * hx
               rfac0m = hx / (6.d0 * hy * (r - (m + 1) * hx))
               if (m .eq. 0) then
                  rfac1m = -rfac0m
               else
                  rfac1m = hx / (6.d0 * hy * (r - m * hx))
               end if
               rfac0p = -rfac1m
               rfac1p = -rfac0m
               src(i,j) = src(i,j) - fac1 *
     @              (idir * (rfac1m * vf(i-m,juf) - rfac0m * vf(i-m-1,juf) +
     @                       rfac1p * vf(i+m,juf) - rfac0p * vf(i+m-1,juf)))
            end if
 100        continue
      end if
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgodiv(src, dims(src),
     @ uc, vc, dims(c),
     @ uf, vf, dims(f),
     @ dims(creg),
     @ hx, hy, ir, jr, idir, jdir, irz)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 uc(dimdec(c))
      real*8 vc(dimdec(c))
      real*8 uf(dimdec(f))
      real*8 vf(dimdec(f))
      real*8 hx, hy
      integer ir, jr, idir, jdir, irz
      real*8 hxm1, hym1, fac0, fac1, r, rfac, rfac0, rfac1
      real*8 rfac0p, rfac1p
      integer i, j, iuf, juf, iuc, juc, m, n
      fac0 = (ir * jr) / (1.5d0*ir*jr + 0.5d0*ir + 0.5d0*jr - 0.5d0)
      hxm1 = 1.d0 / (ir * hx)
      hym1 = 1.d0 / (jr * hy)
      i = cregl0
      j = cregl1
      if (idir .eq. 1) then
         iuc = i - 1
         iuf = i * ir
      else
         iuc = i
         iuf = i * ir - 1
      end if
      if (jdir .eq. 1) then
         juc = j - 1
         juf = j * jr
      else
         juc = j
         juf = j * jr - 1
      end if
      src(i*ir,j*jr) = fac0 *
     @  (hxm1 * idir * (uc(iuc+idir,juc) - uc(iuc,juc) -
     @                  uc(iuc,juc+jdir)) +
     @   hym1 * jdir * (vc(iuc,juc+jdir) - vc(iuc,juc) -
     @                  vc(iuc+idir,juc)))
      if (irz .eq. 1) then
         r = (iuc + 0.5d0) * (hx * ir)
         rfac0 = (ir * hx) / (6.d0 * jr * hy * (r + idir * hx * ir))
         rfac1 = (ir * hx) / (6.d0 * jr * hy * r)
         src(i*ir,j*jr) = src(i*ir,j*jr) + fac0 *
     @     (idir * jdir * (rfac1 * (vc(iuc,juc+jdir) - vc(iuc,juc)) +
     @                     rfac0 *  vc(iuc+idir,juc)))
      end if
      i = i * ir
      j = j * jr
      fac0 = fac0 / (ir * jr)
      hxm1 = ir * hxm1
      hym1 = jr * hym1
      src(i,j) = src(i,j) + fac0 *
     @  (hxm1 * idir * uf(iuf,juf) + hym1 * jdir * vf(iuf,juf))
      if (irz .eq. 1) then
         r = (iuf + 0.5d0) * hx
         rfac = hx / (6.d0 * hy * r)
         src(i,j) = src(i,j) - fac0 *
     @     (idir * jdir * rfac * vf(iuf,juf))
      end if
      fac0 = fac0 / ir
      do 10 m = idir, idir*(ir-1), idir
         fac1 = (ir-abs(m)) * fac0
         src(i,j) = src(i,j) + fac1 *
     @     (hxm1 * idir * (uf(iuf+m,juf) - uf(iuf+m-idir,juf)) +
     @      hym1 * jdir * (vf(iuf+m,juf) + vf(iuf+m-idir,juf)))
         if (irz .eq. 1) then
            r = (iuf + 0.5d0) * hx
            rfac0p = hx / (6.d0 * hy * (r + (m - idir) * hx))
            rfac1p = hx / (6.d0 * hy * (r + m * hx))
            src(i,j) = src(i,j) - fac1 *
     @        (idir * jdir * (rfac1p * vf(iuf+m,juf) -
     @                        rfac0p * vf(iuf+m-idir,juf)))
         end if
 10      continue
      fac0 = ir * fac0 / jr
      do 20 n = jdir, jdir*(jr-1), jdir
         fac1 = (jr-abs(n)) * fac0
         src(i,j) = src(i,j) + fac1 *
     @     (hxm1 * idir * (uf(iuf,juf+n) + uf(iuf,juf+n-jdir)) +
     @      hym1 * jdir * (vf(iuf,juf+n) - vf(iuf,juf+n-jdir)))
         if (irz .eq. 1) then
            r = (iuf + 0.5d0) * hx
            rfac1  = hx / (6.d0 * hy * r)
            src(i,j) = src(i,j) - fac1 *
     @        idir * jdir * rfac1 * (vf(iuf,juf+n) - vf(iuf,juf+n-jdir))
         end if
 20      continue
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgidiv(src, dims(src),
     @ uc, vc, dims(c),
     @ uf, vf, dims(f),
     @ dims(creg),
     @ hx, hy, ir, jr, idir, jdir, irz)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 uc(dimdec(c))
      real*8 vc(dimdec(c))
      real*8 uf(dimdec(f))
      real*8 vf(dimdec(f))
      real*8 hx, hy
      integer ir, jr, idir, jdir, irz
      real*8 hxm1, hym1, fac0, fac1, r, rfac, rfac0, rfac1
      real*8 rfac0m, rfac1m
      integer i, j, iuf, juf, iuc, juc, m, n
      fac0 = (ir * jr) / (0.5d0*ir*jr + 0.5d0*ir + 0.5d0*jr + 0.5d0)
      hxm1 = 1.d0 / (ir * hx)
      hym1 = 1.d0 / (jr * hy)
      i = cregl0
      j = cregl1
      if (idir .eq. 1) then
         iuc = i - 1
         iuf = i * ir
      else
         iuc = i
         iuf = i * ir - 1
      end if
      if (jdir .eq. 1) then
         juc = j - 1
         juf = j * jr
      else
         juc = j
         juf = j * jr - 1
      end if
      src(i*ir,j*jr) = -fac0 *
     @  (hxm1 * idir * uc(iuc,juc) + hym1 * jdir * vc(iuc,juc))
      if (irz .eq. 1) then
         r = (iuc + 0.5d0) * (hx * ir)
         rfac = (ir * hx) / (6.d0 * jr * hy * r)
         src(i*ir,j*jr) = src(i*ir,j*jr) - fac0 *
     @     (idir * jdir * rfac * vc(iuc,juc))
      end if
      i = i * ir
      j = j * jr
      fac0 = fac0 / (ir * jr)
      hxm1 = ir * hxm1
      hym1 = jr * hym1
      src(i,j) = src(i,j) + fac0 *
     @  (hxm1 * idir * (uf(iuf,juf) - uf(iuf-idir,juf) +
     @                  uf(iuf,juf-jdir)) +
     @   hym1 * jdir * (vf(iuf,juf) + vf(iuf-idir,juf) -
     @                  vf(iuf,juf-jdir)))
      if (irz .eq. 1) then
         r = (iuf + 0.5d0) * hx
         rfac0 = hx / (6.d0 * hy * (r - idir * hx))
         rfac1 = hx / (6.d0 * hy * r)
         src(i,j) = src(i,j) + fac0 *
     @     (idir * jdir * (rfac1 * (vf(iuf,juf-jdir) - vf(iuf,juf)) +
     @                     rfac0 *  vf(iuf-idir,juf)))
      end if
      fac0 = fac0 / ir
      do 10 m = idir, idir*(ir-1), idir
         fac1 = (ir-abs(m)) * fac0
         src(i,j) = src(i,j) + fac1 *
     @     (hxm1 * idir * (uf(iuf-m,juf) - uf(iuf-m-idir,juf)) +
     @      hym1 * jdir * (vf(iuf-m,juf) + vf(iuf-m-idir,juf)))
         if (irz .eq. 1) then
            r = (iuf + 0.5d0) * hx
            rfac0m = hx / (6.d0 * hy * (r - (m + idir) * hx))
            rfac1m = hx / (6.d0 * hy * (r - m * hx))
            src(i,j) = src(i,j) - fac1 *
     @        (idir * jdir * (rfac1m * vf(iuf-m,juf) -
     @                        rfac0m * vf(iuf-m-idir,juf)))
         end if
 10      continue
      fac0 = ir * fac0 / jr
      do 20 n = jdir, jdir*(jr-1), jdir
         fac1 = (jr-abs(n)) * fac0
         src(i,j) = src(i,j) + fac1 *
     @     (hxm1 * idir * (uf(iuf,juf-n) + uf(iuf,juf-n-jdir)) +
     @      hym1 * jdir * (vf(iuf,juf-n) - vf(iuf,juf-n-jdir)))
         if (irz .eq. 1) then
            r = (iuf + 0.5d0) * hx
            rfac1  = hx / (6.d0 * hy * r)
            src(i,j) = src(i,j) - fac1 *
     @        idir * jdir * rfac1 * (vf(iuf,juf-n) - vf(iuf,juf-n-jdir))
         end if
 20      continue
      end

c Note---only generates values at coarse points along edge of fine grid
      subroutine hgddiv(src, dims(src),
     @ uc, vc, dims(c),
     @ uf, vf, dims(f),
     @ dims(creg),
     @ hx, hy, ir, jr, jdir, irz)
      integer dims(src)
      integer dims(c)
      integer dims(f)
      integer dims(creg)
      real*8 src(dimdec(src))
      real*8 uc(dimdec(c))
      real*8 vc(dimdec(c))
      real*8 uf(dimdec(f))
      real*8 vf(dimdec(f))
      real*8 hx, hy
      integer ir, jr, jdir, irz
      real*8 hxm1, hym1, fac0, fac1, r, rfac, rfac0, rfac1
      real*8 rfac0m, rfac1m, rfac0p, rfac1p
      integer i, j, juf, juc, m, n
      fac0 = (ir * jr) / (ir*jr + ir + jr - 1.d0)
      hxm1 = 1.d0 / (ir * hx)
      hym1 = 1.d0 / (jr * hy)
      i = cregl0
      j = cregl1
      if (jdir .eq. 1) then
         juc = j
         juf = j * jr
      else
         juc = j - 1
         juf = j * jr - 1
      end if
      src(i*ir,j*jr) = fac0 *
     @  (hxm1 * (uc(i,juc-jdir) - uc(i-1,juc)) +
     @   hym1 * jdir * (vc(i-1,juc) - vc(i,juc-jdir)))
      if (irz .eq. 1) then
         r = (i + 0.5d0) * (hx * ir)
         rfac0 = (ir * hx) / (6.d0 * jr * hy * (r - hx * ir))
         rfac1 = (ir * hx) / (6.d0 * jr * hy * r)
         src(i*ir,j*jr) = src(i*ir,j*jr) + fac0 *
     @     (jdir * (rfac1 * vc(i,juc-jdir) +
     @              rfac0 * vc(i-1,juc)))
      end if
      i = i * ir
      j = j * jr
      fac0 = fac0 / (ir * jr)
      hxm1 = ir * hxm1
      hym1 = jr * hym1
      src(i,j) = src(i,j) + fac0 *
     @  (hxm1 * (uf(i,juf) - uf(i-1,juf-jdir)) +
     @   hym1 * jdir * (vf(i,juf) - vf(i-1,juf-jdir)))
      if (irz .eq. 1) then
         r = (i + 0.5d0) * hx
         rfac0 = hx / (6.d0 * hy * (r - hx))
         rfac1 = hx / (6.d0 * hy * r)
         src(i,j) = src(i,j) - fac0 *
     @     (jdir * (rfac1 * vf(i,juf) +
     @              rfac0 * vf(i-1,juf-jdir)))
      end if
      fac0 = fac0 / ir
      do 10 m = 1, ir-1
         fac1 = (ir-m) * fac0
         src(i,j) = src(i,j) + fac1 *
     @     (hxm1 * (uf(i+m,juf) - uf(i+m-1,juf) +
     @              uf(i-m,juf-jdir) - uf(i-m-1,juf-jdir)) +
     @      hym1 * jdir * (vf(i+m,juf) + vf(i+m-1,juf) -
     @                     vf(i-m,juf-jdir) - vf(i-m-1,juf-jdir)))
         if (irz .eq. 1) then
            r = (i + 0.5d0) * hx
            rfac0m = hx / (6.d0 * hy * (r - (m + 1) * hx))
            rfac1m = hx / (6.d0 * hy * (r - m * hx))
            rfac0p = hx / (6.d0 * hy * (r + (m - 1) * hx))
            rfac1p = hx / (6.d0 * hy * (r + m * hx))
            src(i,j) = src(i,j) - fac1 *
     @        (jdir * (rfac1p * vf(i+m,juf) -
     @                 rfac0p * vf(i+m-1,juf) -
     @                 rfac1m * vf(i-m,juf-jdir) +
     @                 rfac0m * vf(i-m-1,juf-jdir)))
         end if
 10      continue
      fac0 = ir * fac0 / jr
      do 20 n = jdir, jdir*(jr-1), jdir
         fac1 = (jr-abs(n)) * fac0
         src(i,j) = src(i,j) + fac1 *
     @     (hxm1 * (uf(i,juf+n) + uf(i,juf+n-jdir) -
     @              uf(i-1,juf-n) - uf(i-1,juf-n-jdir)) +
     @      hym1 * jdir * (vf(i,juf+n) - vf(i,juf+n-jdir) +
     @                     vf(i-1,juf-n) - vf(i-1,juf-n-jdir)))
         if (irz .eq. 1) then
            r = (i + 0.5d0) * hx
            rfac0  = hx / (6.d0 * hy * (r - hx))
            rfac1  = hx / (6.d0 * hy * r)
            src(i,j) = src(i,j) - fac1 *
     @        (jdir * (rfac1 * (vf(i,juf+n) - vf(i,juf+n-jdir)) -
     @                 rfac0 * (vf(i-1,juf-n) - vf(i-1,juf-n-jdir))))
         end if
 20      continue
      end

#endif
