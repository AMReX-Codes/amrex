
#define dims(a) a l0, a h0, a l1, a h1
#define dimdec(a) a l0:a h0,a l1:a h1
c#define dimdecu(a) (a h0-a l0+1)*(a h1-a l1+1)
#define dimdecu(a) *

#include <REAL.H>

c Works for CELL- or NODE-based data.
      subroutine fcpy(dest, dims(dest),
     @ dims(reg),
     @ src,  dims(src), ncomp)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer ncomp
      real*8 dest(dimdec(dest), ncomp)
      real*8 src(dimdec(src), ncomp)
      integer i, j, k
      if (regh1 - regl1 .gt. regh0 - regl0) then
         do 10 k = 1, ncomp
            do 10 i = regl0, regh0
               do 10 j = regl1, regh1
 10               dest(i,j,k) = src(i,j,k)
      else
         do 20 k = 1, ncomp
            do 20 j = regl1, regh1
               do 20 i = regl0, regh0
 20               dest(i,j,k) = src(i,j,k)
      endif
      end

c Unrolled copy---works for CELL- or NODE-based data.
      subroutine fcpyu(dest, src,
     @ dims(reg),
     @ ncomp)
      integer dims(reg)
      real*8 dest(*)
      real*8 src(*)
      integer ncomp
      integer nvals, i
      nvals = (regh0-regl0+1) * (regh1-regl1+1) * ncomp
      do 10 i = 1, nvals
 10      dest(i) = src(i)
      end

c Works for NODE-based data.
      subroutine fcpy2(da, dims(da),
     @ db,  dims(db),
     @ dims(reg), w, ncomp)
      integer dims(da)
      integer dims(db)
      integer dims(reg)
      integer w, ncomp
      real*8 da(dimdec(da), ncomp)
      real*8 db(dimdec(db), ncomp)
      integer i, j, k
      if (regh1 - regl1 .gt. regh0 - regl0) then
         do 10 k = 1, ncomp
            do 10 i = 1, w
               do 10 j = regl1, regh1
                  da(regl0+i,j,k) = db(regl0+i,j,k)
 10               db(regl0-i,j,k) = da(regl0-i,j,k)
      else
         do 20 k = 1, ncomp
            do 20 j = 1, w
               do 20 i = regl0 - w, regh0 + w
                  da(i,regl1+j,k) = db(i,regl1+j,k)
 20               db(i,regl1-j,k) = da(i,regl1-j,k)
      endif
      end

      subroutine fccpy1(dest, src, nsets, dstart, sstart,
     @                  dstrid, sstrid, nvals)
      real*8 dest(0:1), src(0:1)
      integer nsets
      integer dstart(nsets), sstart(nsets)
      integer dstrid(nsets), sstrid(nsets), nvals(nsets)
      integer i, k
      do 10 k = 1, nsets
         do 10 i = 0, nvals(k) - 1
 10         dest(dstart(k) + i * dstrid(k)) =
     @       src(sstart(k) + i * sstrid(k))
      end

c Works for CELL- or NODE-based data.
      subroutine mxnorm(src, dims(src),
     @ dims(reg),
     @ retval)
      integer dims(src)
      integer dims(reg)
      real*8 src(dimdec(src))
      real*8 retval
      integer bufsiz
      parameter(bufsiz = 128)
      real*8 buf(0:bufsiz-1)
      integer i, j, k, ii, jj, idif, jdif
      idif = regh0 - regl0
      jdif = regh1 - regl1
      if (jdif .gt. idif) then
         k = jdif - mod(jdif, bufsiz)
         do 10 j = 0, min(bufsiz - 1, jdif)
 10         buf(j) = 0.d0
         do 20 i = regl0, regh0
            do 20 jj = 0, k, bufsiz
               do 20 j = 0, min(bufsiz - 1, jdif - jj)
 20               buf(j) = max(buf(j), abs(src(i, regl1 + jj + j)))
         do 30 j = 0, min(bufsiz - 1, jdif)
 30         retval = max(retval, buf(j))
      else
         k = idif - mod(idif, bufsiz)
         do 40 i = 0, min(bufsiz - 1, idif)
 40         buf(i) = 0.d0
         do 50 j = regl1, regh1
            do 50 ii = 0, k, bufsiz
               do 50 i = 0, min(bufsiz - 1, idif - ii)
 50               buf(i) = max(buf(i), abs(src(regl0 + ii + i, j)))
         do 60 i = 0, min(bufsiz - 1, idif)
 60         retval = max(retval, buf(i))
      endif
      end

c Works for CELL- or NODE-based data.
      subroutine fassnu(v, x,
     @ dims(reg),
     @ ncomp)
      integer dims(reg)
      real*8 v(dimdecu(reg))
      real*8 x
      integer ncomp
      integer nvals, i
      nvals = (regh0-regl0+1) * (regh1-regl1+1) * ncomp
      do 10 i = 1, nvals
 10      v(i) = x
      end

c Works for CELL- or NODE-based data.
      subroutine fplusu(v0, v1,
     @ dims(reg))
      integer dims(reg)
      real*8 v0(dimdecu(reg))
      real*8 v1(dimdecu(reg))
      integer i
      do 10 i = 1, (regh0 - regl0 + 1) * (regh1 - regl1 + 1)
 10      v0(i) = v0(i) + v1(i)
      end

c Works for CELL-based data.
      subroutine iprodc(v0, dims(v0),
     @ v1, dims(v1),
     @ dims(reg), sum)
      integer dims(v0)
      integer dims(v1)
      integer dims(reg)
      real*8 v0(dimdec(v0))
      real*8 v1(dimdec(v1))
      real*8 sum
      integer i, j
      if (regh1 - regl1 .gt. regh0 - regl0) then
         do 10 i = regl0, regh0
            do 10 j = regl1, regh1
 10            sum = sum + v0(i,j) * v1(i,j)
      else
         do 20 j = regl1, regh1
            do 20 i = regl0, regh0
 20            sum = sum + v0(i,j) * v1(i,j)
      endif
      end

c Works for NODE-based data.
      subroutine iprodn(v0, dims(v0),
     @ v1, dims(v1),
     @ dims(reg), sum)
      integer dims(v0)
      integer dims(v1)
      integer dims(reg)
      real*8 v0(dimdec(v0))
      real*8 v1(dimdec(v1))
      real*8 sum, sum0
      integer i, j
      sum0 = 0.5d0 * (v0(regl0,regl1) * v1(regl0,regl1) +
     @                v0(regl0,regh1) * v1(regl0,regh1) +
     @                v0(regh0,regl1) * v1(regh0,regl1) +
     @                v0(regh0,regh1) * v1(regh0,regh1))
      do 10 i = regl0 + 1, regh0 - 1
 10      sum0 = sum0 + (v0(i,regl1) * v1(i,regl1) +
     @                  v0(i,regh1) * v1(i,regh1))
      do 20 j = regl1 + 1, regh1 - 1
 20      sum0 = sum0 + (v0(regl0,j) * v1(regl0,j) +
     @                  v0(regh0,j) * v1(regh0,j))
      sum = sum + 0.5d0 * sum0
      if (regh1 - regl1 .gt. regh0 - regl0) then
         do 30 i = regl0 + 1, regh0 - 1
            do 30 j = regl1 + 1, regh1 - 1
 30            sum = sum + v0(i,j) * v1(i,j)
      else
         do 40 j = regl1 + 1, regh1 - 1
            do 40 i = regl0 + 1, regh0 - 1
 40            sum = sum + v0(i,j) * v1(i,j)
      endif
      end

c Works for CELL- or NODE-based data.
      subroutine bref(dest, dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ dims(bb), idir)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer dims(bb)
      integer idir
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      integer i, j
      if (idir .eq. 0) then
         do 10 i = regl0, regh0
            do 10 j = regl1, regh1
 10            dest(i,j) = src(bbh0-(i-regl0),j)
      else
         do 20 j = regl1, regh1
            do 20 i = regl0, regh0
 20            dest(i,j) = src(i,bbh1-(j-regl1))
      endif
      end

c Works for CELL- or NODE-based data.
      subroutine brefm(dest, dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ dims(bb), ra)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer dims(bb)
      integer ra(0:1)
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      integer i, j
      if (ra(0) .eq. 0 .and. ra(1) .eq. 0) then
         if (regh1 - regl1 .gt. regh0 - regl0) then
            do 10 i = regl0, regh0
               do 10 j = regl1, regh1
 10               dest(i,j) = src(bbl0+(i-regl0),bbl1+(j-regl1))
         else
            do 12 j = regl1, regh1
               do 12 i = regl0, regh0
 12               dest(i,j) = src(bbl0+(i-regl0),bbl1+(j-regl1))
         endif
      else if (ra(0) .eq. 0 .and. ra(1) .eq. 1) then
         do 20 j = regl1, regh1
            do 20 i = regl0, regh0
 20            dest(i,j) = src(bbl0+(i-regl0),bbh1-(j-regl1))
      else if (ra(0) .eq. 1 .and. ra(1) .eq. 0) then
         do 30 i = regl0, regh0
            do 30 j = regl1, regh1
 30            dest(i,j) = src(bbh0-(i-regl0),bbl1+(j-regl1))
      else
         do 40 j = regl1, regh1
            do 40 i = regl0, regh0
 40            dest(i,j) = src(bbh0-(i-regl0),bbh1-(j-regl1))
      endif
      end

c Works for CELL- or NODE-based data.
      subroutine bneg(dest, dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ dims(bb), idir)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer dims(bb)
      integer idir
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      integer i, j
      if (idir .eq. 0) then
         do 10 i = regl0, regh0
            do 10 j = regl1, regh1
 10            dest(i,j) = -src(bbh0-(i-regl0),j)
      else
         do 20 j = regl1, regh1
            do 20 i = regl0, regh0
 20            dest(i,j) = -src(i,bbh1-(j-regl1))
      endif
      end

c Works for CELL- or NODE-based data.
      subroutine bnegm(dest, dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ dims(bb), ra)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer dims(bb)
      integer ra(0:1)
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      integer i, j
      if (ra(0) .eq. 0 .and. ra(1) .eq. 0) then
         if (regh1 - regl1 .gt. regh0 - regl0) then
            do 10 i = regl0, regh0
               do 10 j = regl1, regh1
 10               dest(i,j) = -src(bbl0+(i-regl0),bbl1+(j-regl1))
         else
            do 12 j = regl1, regh1
               do 12 i = regl0, regh0
 12               dest(i,j) = -src(bbl0+(i-regl0),bbl1+(j-regl1))
         endif
      else if (ra(0) .eq. 0 .and. ra(1) .eq. 1) then
         do 20 j = regl1, regh1
            do 20 i = regl0, regh0
 20            dest(i,j) = -src(bbl0+(i-regl0),bbh1-(j-regl1))
      else if (ra(0) .eq. 1 .and. ra(1) .eq. 0) then
         do 30 i = regl0, regh0
            do 30 j = regl1, regh1
 30            dest(i,j) = -src(bbh0-(i-regl0),bbl1+(j-regl1))
      else
         do 40 i = regl0, regh0
            do 40 j = regl1, regh1
 40            dest(i,j) = -src(bbh0-(i-regl0),bbh1-(j-regl1))
      endif
      end

c Works for CELL-based velocity data.
c This routine assumes that the inflow face velocity data has not yet
c been altered.  Running fill_borders should call this routine on every
c inflow face, so that binfil can be run for subsequent fills
      subroutine binflo(dest, dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ dims(bb),
     @ idir)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer dims(bb)
      integer idir
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      integer i, j
      if (idir .eq. 0) then
         if (regl0 .lt. bbh0) then
            do 10 i = regl0, regh0
               do 10 j = regl1, regh1
 10               dest(i,j) = 2.d0 * dest(regh0,j) - src(bbh0-(i-regl0),j)
         else
            do 20 i = regh0, regl0, -1
               do 20 j = regl1, regh1
 20               dest(i,j) = 2.d0 * dest(regl0,j) - src(bbh0-(i-regl0),j)
         endif
      else
         if (regl1 .lt. bbh1) then
            do 30 j = regl1, regh1
               do 30 i = regl0, regh0
 30               dest(i,j) = 2.d0 * dest(i,regh1) - src(i,bbh1-(j-regl1))
         else
            do 40 j = regh1, regl1, -1
               do 40 i = regl0, regh0
 40               dest(i,j) = 2.d0 * dest(i,regl1) - src(i,bbh1-(j-regl1))
         endif
      endif
      end

c Works for CELL-based velocity data.
c This routine is called when the inflow face velocity data has already
c been altered by a call to fill_borders.  The box bb must have been
c extended to one cell past the boundary by the boundary::box routine.

      subroutine binfil(dest, dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ dims(bb),
     @ idir)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer dims(bb)
      integer idir
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      integer i, j
      if (idir .eq. 0) then
         if (regl0 .lt. bbh0) then
            do 10 i = regl0, regh0
               do 10 j = regl1, regh1
 10               dest(i,j) = src(bbl0,j) + src(bbl0+1,j) -
     @                        src(bbh0-(i-regl0),j)
         else
            do 20 i = regh0, regl0, -1
               do 20 j = regl1, regh1
 20               dest(i,j) = src(bbh0-1,j) + src(bbh0,j) -
     @                        src(bbl0+(regh0-i),j)
         endif
      else
         if (regl1 .lt. bbh1) then
            do 30 j = regl1, regh1
               do 30 i = regl0, regh0
 30               dest(i,j) = src(i,bbl1) + src(i,bbl1+1) -
     @                        src(i,bbh1-(j-regl1))
         else
            do 40 j = regh1, regl1, -1
               do 40 i = regl0, regh0
 40               dest(i,j) = src(i,bbh1-1) + src(i,bbh1) -
     @                        src(i,bbl1+(regh1-j))
         endif
      endif
      end

c CELL-based data only.
      subroutine acint2(dest, dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ dims(bb), ir, jr)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer dims(bb)
      integer ir, jr
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      real*8 xoff, yoff, sx, sy
      integer i, j, ic, jc
      do 10 j = regl1, regh1
         jc = j/jr
         yoff = (mod(j,jr) + 0.5) / jr - 0.5
         do 10 i = regl0, regh0
            ic = i/ir
            xoff = (mod(i,ir) + 0.5) / ir - 0.5
            sy = 0.5 * (src(ic,jc+1) - src(ic,jc-1))
            sx = 0.5 * (src(ic+1,jc) - src(ic-1,jc))
 10         dest(i,j) = src(ic,jc) + xoff * sx + yoff * sy
      end

c NODE-based data only.
      subroutine anint2(dest, dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ dims(bb), ir, jr)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer dims(bb)
      integer ir, jr
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      real*8 p, q
      integer ic, jc, j, m
      do 10 jc = bbl1, bbh1
         do 10 ic = bbl0, bbh0
 10         dest(ir*ic,jr*jc) = src(ic,jc)
      if (jr .eq. 2) then
         do 20 jc = bbl1, bbh1-1
            do 20 ic = bbl0, bbh0
 20            dest(ir*ic,2*jc+1) = 0.5d0 * (src(ic,jc) +
     @                                       src(ic,jc+1))
      else if (jr .eq. 4) then
         do 120 jc = bbl1, bbh1-1
cdir$ ivdep
            do 120 ic = bbl0, bbh0
               dest(ir*ic,4*jc+2) = 0.5d0 *
     @            (src(ic,jc) + src(ic,jc+1))
               dest(ir*ic,4*jc+1) = 0.5d0 *
     @            (dest(ir*ic,4*jc)   + dest(ir*ic,4*jc+2))
 120           dest(ir*ic,4*jc+3) = 0.5d0 *
     @            (dest(ir*ic,4*jc+2) + dest(ir*ic,4*jc+4))
      else if (jr .ne. 1) then
         do 220 m = 1, jr-1
            q = dble(m)/jr
            p = 1.0 - q
            do 220 jc = bbl1, bbh1-1
               do 220 ic = bbl0, bbh0
 220              dest(ir*ic,jr*jc+m) = p * src(ic,jc) +
     @                                  q * src(ic,jc+1)
      endif
      if (ir .eq. 2) then
         do 30 ic = bbl0, bbh0-1
cdir$ ivdep
            do 30 j = regl1, regh1
 30            dest(2*ic+1,j) = 0.5d0 * (dest(2*ic,j) +
     @                                   dest(2*ic+2,j))
      else if (ir .eq. 4) then
         do 130 ic = bbl0, bbh0-1
cdir$ ivdep
            do 130 j = regl1, regh1
               dest(4*ic+2,j) = 0.5d0 *
     @            (dest(4*ic,j)   + dest(4*ic+4,j))
               dest(4*ic+1,j) = 0.5d0 *
     @            (dest(4*ic,j)   + dest(4*ic+2,j))
 130           dest(4*ic+3,j) = 0.5d0 *
     @            (dest(4*ic+2,j) + dest(4*ic+4,j))
      else if (ir .ne. 1) then
         do 230 m = 1, ir-1
            q = dble(m)/ir
            p = 1.0 - q
            do 230 ic = bbl0, bbh0-1
cdir$ ivdep
               do 230 j = regl1, regh1
 230              dest(ir*ic+m,j) = p * dest(ir*ic,j) +
     @                              q * dest(ir*(ic+1),j)
      endif
      end

c CELL-based data only.
      subroutine acrst1(dest, dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ ir, jr, integ)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer ir, jr, integ
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      real*8 fac
      integer i, j, m, n
      do 10 j = regl1, regh1
         do 10 i = regl0, regh0
 10         dest(i,j) = 0.d0
      do 20 n = 0, jr-1
         do 20 m = 0, ir-1
            do 20 j = regl1, regh1
               do 20 i = regl0, regh0
 20               dest(i,j) = dest(i,j) + src(i*ir+m,j*jr+n)
      if (integ .eq. 0) then
         fac = 1.d0 / (ir*jr)
         do 30 j = regl1, regh1
            do 30 i = regl0, regh0
 30            dest(i,j) = dest(i,j) * fac
      endif
      end

c NODE-based data only.
      subroutine anrst1(dest, dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ ir, jr)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer ir, jr
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      integer i, j
      do 10 j = regl1, regh1
         do 10 i = regl0, regh0
 10         dest(i,j) = src(i*ir,j*jr)
      end

#define RESTR4 0.0625d0 * src(4*i,4*j) + \
      0.046875d0  * (src(4*i,4*j-1) + src(4*i,4*j+1) + \
                     src(4*i-1,4*j) + src(4*i+1,4*j)) + \
      0.03125d0   * (src(4*i,4*j-2) + src(4*i,4*j+2) + \
                     src(4*i-2,4*j) + src(4*i+2,4*j)) + \
      0.015625d0  * (src(4*i,4*j-3) + src(4*i,4*j+3) + src(4*i-3,4*j) + \
                     src(4*i+3,4*j) + src(4*i-2,4*j-2) + src(4*i-2,4*j+2) + \
                     src(4*i+2,4*j-2) + src(4*i+2,4*j+2)) + \
      0.03515625d0* (src(4*i-1,4*j-1) + src(4*i-1,4*j+1) + \
                     src(4*i+1,4*j-1) + src(4*i+1,4*j+1)) + \
      0.0234375d0 * (src(4*i-2,4*j-1) + src(4*i-2,4*j+1) + src(4*i-1,4*j-2) + \
                     src(4*i-1,4*j+2) + src(4*i+1,4*j-2) + src(4*i+1,4*j+2) + \
                     src(4*i+2,4*j-1) + src(4*i+2,4*j+1)) + \
      0.01171875d0* (src(4*i-3,4*j-1) + src(4*i-3,4*j+1) + src(4*i-1,4*j-3) + \
                     src(4*i-1,4*j+3) + src(4*i+1,4*j-3) + src(4*i+1,4*j+3) + \
                     src(4*i+3,4*j-1) + src(4*i+3,4*j+1)) + \
      0.0078125d0 * (src(4*i-3,4*j-2) + src(4*i-3,4*j+2) + src(4*i-2,4*j-3) + \
                     src(4*i-2,4*j+3) + src(4*i+2,4*j-3) + src(4*i+2,4*j+3) + \
                     src(4*i+3,4*j-2) + src(4*i+3,4*j+2)) + \
      0.00390625d0* (src(4*i-3,4*j-3) + src(4*i-3,4*j+3) + \
                     src(4*i+3,4*j-3) + src(4*i+3,4*j+3))

c NODE-based data only.
c Fills coarse region defined by reg, which must be smaller than
c fine region by at least one coarse cell on all sides.
      subroutine anrst2(dest, dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ ir, jr, integ)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer ir, jr, integ
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      real*8 fac0, fac1, fac
      integer i, j, m, n
      if (ir .eq. 2 .and. jr .eq. 2) then
         if (regh1 - regl1 .gt. regh0 - regl0) then
            do 10 i = regl0, regh0
               do 10 j = regl1, regh1
 10               dest(i,j) = 0.25d0   *  src(i*ir,j*jr) +
     @                        0.125d0  * (src(i*ir,j*jr-1) +
     @                                    src(i*ir,j*jr+1) +
     @                                    src(i*ir-1,j*jr) +
     @                                    src(i*ir+1,j*jr)) +
     @                        0.0625d0 * (src(i*ir-1,j*jr-1) +
     @                                    src(i*ir-1,j*jr+1) +
     @                                    src(i*ir+1,j*jr-1) +
     @                                    src(i*ir+1,j*jr+1))
         else
            do 20 j = regl1, regh1
               do 20 i = regl0, regh0
 20               dest(i,j) = 0.25d0   *  src(i*ir,j*jr) +
     @                        0.125d0  * (src(i*ir,j*jr-1) +
     @                                    src(i*ir,j*jr+1) +
     @                                    src(i*ir-1,j*jr) +
     @                                    src(i*ir+1,j*jr)) +
     @                        0.0625d0 * (src(i*ir-1,j*jr-1) +
     @                                    src(i*ir-1,j*jr+1) +
     @                                    src(i*ir+1,j*jr-1) +
     @                                    src(i*ir+1,j*jr+1))
         endif
      else if (ir .eq. 4 .and. jr .eq. 4) then
         if (regh1 - regl1 .gt. regh0 - regl0) then
            do 110 i = regl0, regh0
               do 110 j = regl1, regh1
 110              dest(i,j) = RESTR4
         else
            do 120 j = regl1, regh1
               do 120 i = regl0, regh0
 120              dest(i,j) = RESTR4
         endif
      else
         if (regh1 - regl1 .gt. regh0 - regl0) then
            do 210 i = regl0, regh0
               do 210 j = regl1, regh1
 210              dest(i,j) = 0.d0
            fac0 = 1.d0 / (ir*ir * jr*jr)
            do 220 n = 0, jr-1
               fac1 = (jr-n) * fac0
               if (n .eq. 0) fac1 = 0.5d0 * fac1
               do 220 m = 0, ir-1
                  fac = (ir-m) * fac1
                  if (m .eq. 0) fac = 0.5d0 * fac
                  do 220 i = regl0, regh0
                     do 220 j = regl1, regh1
 220                    dest(i,j) = dest(i,j) +
     @                              fac * (src(i*ir-m,j*jr-n)+
     @                                     src(i*ir-m,j*jr+n)+
     @                                     src(i*ir+m,j*jr-n)+
     @                                     src(i*ir+m,j*jr+n))
         else
            do 230 j = regl1, regh1
               do 230 i = regl0, regh0
 230              dest(i,j) = 0.d0
            fac0 = 1.d0 / (ir*ir * jr*jr)
            do 240 n = 0, jr-1
               fac1 = (jr-n) * fac0
               if (n .eq. 0) fac1 = 0.5d0 * fac1
               do 240 m = 0, ir-1
                  fac = (ir-m) * fac1
                  if (m .eq. 0) fac = 0.5d0 * fac
                  do 240 j = regl1, regh1
                     do 240 i = regl0, regh0
 240                    dest(i,j) = dest(i,j) +
     @                              fac * (src(i*ir-m,j*jr-n)+
     @                                     src(i*ir-m,j*jr+n)+
     @                                     src(i*ir+m,j*jr-n)+
     @                                     src(i*ir+m,j*jr+n))
         endif
      endif
      if (integ .eq. 1) then
         fac = ir * jr
         if (regh1 - regl1 .gt. regh0 - regl0) then
            do 310 i = regl0, regh0
               do 310 j = regl1, regh1
 310              dest(i,j) = fac * dest(i,j)
         else
            do 320 j = regl1, regh1
               do 320 i = regl0, regh0
 320              dest(i,j) = fac * dest(i,j)
         endif
      endif
      end

c NODE-based data only.
c Fills coarse region defined by reg.
      subroutine anfr2(dest, dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ ir, jr, idim, idir, integ)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer ir, jr, idim, idir, integ
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      real*8 fac0, fac1, fac
      integer i, j, m, n
      if (idim .eq. 0) then
         if (integ .eq. 0) then
            fac = (0.5d0 + 0.5d0 / ir)
            fac0 = 1.d0 / (ir*ir * jr*jr)
         else
            fac = 1.d0
            fac0 = 1.d0 / (ir * jr)
         endif
         i = regl0
         do 10 j = regl1, regh1
 10         dest(i,j) = fac * src(i*ir,j*jr)
         do 20 n = 0, jr-1
            fac1 = (jr-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do 20 m = idir, idir*(ir-1), idir
               fac = (ir-abs(m)) * fac1
               do 20 j = regl1, regh1
 20               dest(i,j) = dest(i,j) + fac * (src(i*ir+m,j*jr-n)+
     @                                           src(i*ir+m,j*jr+n))
      else
         if (integ .eq. 0) then
            fac = (0.5d0 + 0.5d0 / jr)
            fac0 = 1.d0 / (ir*ir * jr*jr)
         else
            fac = 1.d0
            fac0 = 1.d0 / (ir * jr)
         endif
         j = regl1
         do 30 i = regl0, regh0
 30         dest(i,j) = fac * src(i*ir,j*jr)
         do 40 n = idir, idir*(jr-1), idir
            fac1 = (jr-abs(n)) * fac0
            do 40 m = 0, ir-1
               fac = (ir-m) * fac1
               if (m .eq. 0) fac = 0.5d0 * fac
               do 40 i = regl0, regh0
 40               dest(i,j) = dest(i,j) + fac * (src(i*ir-m,j*jr+n)+
     @                                           src(i*ir+m,j*jr+n))
      endif
      end

c NODE-based data only.
c Fills coarse point defined by reg.
      subroutine anor2(dest, dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ ir, jr, idir0, idir1, integ)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer ir, jr, idir0, idir1, integ
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      real*8 fac0, fac1, fac
      integer i, j, m, n
      if (integ .eq. 0) then
         fac = (0.75d0 + 0.25d0 / ir + 0.25d0 / jr - 0.25d0 / (ir*jr))
         fac0 = 1.d0 / (ir*ir * jr*jr)
      else
         fac = 1.d0
         fac0 = 1.d0 / (ir * jr)
      endif
      i = regl0
      j = regl1
      dest(i,j) = fac * src(i*ir,j*jr)
      do 20 n = idir1, idir1*(jr-1), idir1
         fac1 = (jr-abs(n)) * fac0
         do 20 m = idir0, idir0*(ir-1), idir0
            fac = (ir-abs(m)) * fac1
 20         dest(i,j) = dest(i,j) + fac * src(i*ir+m,j*jr+n)
      end

c NODE-based data only.
c Fills coarse point defined by reg.
      subroutine anir2(dest, dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ ir, jr, idir0, idir1, integ)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer ir, jr, idir0, idir1, integ
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      real*8 fac0, fac1, fac
      integer i, j, m, n
      if (integ .eq. 0) then
         fac = (0.25d0 + 0.25d0 / ir + 0.25d0 / jr + 0.25d0 / (ir*jr))
         fac0 = 1.d0 / (ir*ir * jr*jr)
      else
         fac = 1.d0
         fac0 = 1.d0 / (ir * jr)
      endif
      i = regl0
      j = regl1
      dest(i,j) = fac * src(i*ir,j*jr)
      do 10 m = idir0, idir0*(ir-1), idir0
         fac1 = (ir-abs(m)) * fac0
         fac = jr * fac1
 10      dest(i,j) = dest(i,j) + fac * src(i*ir+m,j*jr)
      do 20 n = idir1, idir1*(jr-1), idir1
         fac1 = (jr-abs(n)) * fac0
         fac = ir * fac1
         dest(i,j) = dest(i,j) + fac * src(i*ir,j*jr+n)
         do 20 m = idir0, idir0*(ir-1), idir0
            fac = (ir-abs(m)) * fac1
 20         dest(i,j) = dest(i,j) + fac * (src(i*ir+m,j*jr+n) +
     @                                     src(i*ir-m,j*jr+n) +
     @                                     src(i*ir+m,j*jr-n))
      end

c NODE-based data only.
c Fills coarse point defined by reg.
      subroutine andr2(dest, dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ ir, jr, idir1, integ)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer ir, jr, idir1, integ
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      real*8 fac0, fac1, fac
      integer i, j, m, n
      if (integ .eq. 0) then
         fac = (0.5d0 + 0.5d0 / ir + 0.5d0 / jr - 0.5d0 / (ir*jr))
         fac0 = 1.d0 / (ir*ir * jr*jr)
      else
         fac = 1.d0
         fac0 = 1.d0 / (ir * jr)
      endif
      i = regl0
      j = regl1
      dest(i,j) = fac * src(i*ir,j*jr)
      do 20 n = idir1, idir1*(jr-1), idir1
         fac1 = (jr-abs(n)) * fac0
         do 20 m = 1, ir-1
            fac = (ir-m) * fac1
 20         dest(i,j) = dest(i,j) + fac * (src(i*ir+m,j*jr+n) +
     @                                     src(i*ir-m,j*jr-n))
      end
