
#define dims(a) a l0,a h0,a l1,a h1,a l2,a h2
#define dimdec(a) a l0:a h0,a l1:a h1,a l2:a h2
c#define dimdecu(a) (a h0-a l0+1)*(a h1-a l1+1)*(a h2-a l2+1)
#define dimdecu(a) *

c Works for CELL- or NODE-based data.
      subroutine fcpy(dest,
     @ dims(dest),
     @ dims(reg),
     @ src,  dims(src), ncomp)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
c      real*8 dest(dimdec(dest), ncomp)
c      real*8 src(dimdec(src), ncomp)
      real*8 dest(dimdec(dest), *)
      real*8 src(dimdec(src), *)
      integer ncomp
      integer i, j, k, ic
      if (regh2 - regl2 .gt. regh1 - regl1 .and.
     @    regh2 - regl2 .gt. regh0 - regl0) then
         do 10 ic = 1, ncomp
            do 10 j = regl1, regh1
               do 10 i = regl0, regh0
                  do 10 k = regl2, regh2
 10                  dest(i,j,k,ic) = src(i,j,k,ic)
      else if (regh1 - regl1 .gt. regh0 - regl0) then
         do 20 ic = 1, ncomp
            do 20 k = regl2, regh2
               do 20 i = regl0, regh0
                  do 20 j = regl1, regh1
 20                  dest(i,j,k,ic) = src(i,j,k,ic)
      else
         do 30 ic = 1, ncomp
            do 30 k = regl2, regh2
               do 30 j = regl1, regh1
                  do 30 i = regl0, regh0
 30                  dest(i,j,k,ic) = src(i,j,k,ic)
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
      nvals = (regh0-regl0+1) * (regh1-regl1+1) * (regh2-regl2+1) * ncomp
      do 10 i = 1, nvals
 10      dest(i) = src(i)
      end

#if 0
c Old version, failed at certain coarse-fine corners in 3D.
c Works for NODE-based data.
      subroutine fcpy2(da,
     @ dims(da),
     @ db,  dims(db),
     @ dims(reg), w, ncomp)
      implicit real*8(a-h,o-z), integer(i-n)
      integer dims(da)
      integer dims(db)
      integer dims(reg)
      real*8 da(dimdec(da), ncomp)
      real*8 db(dimdec(db), ncomp)
      integer w, ncomp
      if (regl0 .eq. regh0) then
         do 10 ic = 1, ncomp
            do 10 i = 1, w
               do 10 k = regl2, regh2
                  do 10 j = regl1, regh1
c               do 10 k = regl2 - w, regh2 + w
c                  do 10 j = regl1 - w, regh1 + w
                     da(regl0+i,j,k,ic) = db(regl0+i,j,k,ic)
 10                  db(regl0-i,j,k,ic) = da(regl0-i,j,k,ic)
      else if (regl1 .eq. regh1) then
         do 20 ic = 1, ncomp
            do 20 j = 1, w
               do 20 k = regl2, regh2
c               do 20 k = regl2 - w, regh2 + w
                  do 20 i = regl0 - w, regh0 + w
                     da(i,regl1+j,k,ic) = db(i,regl1+j,k,ic)
 20                  db(i,regl1-j,k,ic) = da(i,regl1-j,k,ic)
      else
         do 30 ic = 1, ncomp
            do 30 k = 1, w
               do 30 j = regl1 - w, regh1 + w
                  do 30 i = regl0 - w, regh0 + w
                     da(i,j,regl2+k,ic) = db(i,j,regl2+k,ic)
 30                  db(i,j,regl2-k,ic) = da(i,j,regl2-k,ic)
      endif
      end
#endif

c Works for NODE-based data.
      subroutine fcpy2(da,
     @ dims(da),
     @ db,  dims(db),
     @ dims(reg), w, b, ncomp)
      integer dims(da)
      integer dims(db)
      integer dims(reg)
      integer w, b, ncomp
      real*8 da(dimdec(da), ncomp)
      real*8 db(dimdec(db), ncomp)
      integer i, j, k, ic, ial, iah, ibl, ibh, jal, jah, jbl, jbh
      if (regl0 .eq. regh0) then
         do 10 ic = 1, ncomp
            do 10 i = 1, w
               do 10 k = regl2, regh2
                  do 10 j = regl1, regh1
                     da(regl0+i,j,k,ic) = db(regl0+i,j,k,ic)
 10                  db(regl0-i,j,k,ic) = da(regl0-i,j,k,ic)
      else if (regl1 .eq. regh1) then
         ial = 0
         iah = 0
         ibl = 0
         ibh = 0
         if (dal0 + b .eq. regl0) ial = w
         if (dah0 - b .eq. regh0) iah = w
         if (dbl0 + b .eq. regl0) ibl = w
         if (dbh0 - b .eq. regh0) ibh = w
         do 30 ic = 1, ncomp
            do 30 j = 1, w
               do 30 k = regl2, regh2
                  do 20 i = regl0 - ial, regh0 + iah
 20                  da(i,regl1+j,k,ic) = db(i,regl1+j,k,ic)
                  do 30 i = regl0 - ibl, regh0 + ibh
 30                  db(i,regl1-j,k,ic) = da(i,regl1-j,k,ic)
      else
         ial = 0
         iah = 0
         ibl = 0
         ibh = 0
         if (dal0 + b .eq. regl0) ial = w
         if (dah0 - b .eq. regh0) iah = w
         if (dbl0 + b .eq. regl0) ibl = w
         if (dbh0 - b .eq. regh0) ibh = w
         jal = 0
         jah = 0
         jbl = 0
         jbh = 0
         if (dal1 + b .eq. regl1) jal = w
         if (dah1 - b .eq. regh1) jah = w
         if (dbl1 + b .eq. regl1) jbl = w
         if (dbh1 - b .eq. regh1) jbh = w
         do 50 ic = 1, ncomp
            do 50 k = 1, w
               do 40 j = regl1 - jal, regh1 + jah
                  do 40 i = regl0 - ial, regh0 + iah
 40                  da(i,j,regl2+k,ic) = db(i,j,regl2+k,ic)
               do 50 j = regl1 - jbl, regh1 + jbh
                  do 50 i = regl0 - ibl, regh0 + ibh
 50                  db(i,j,regl2-k,ic) = da(i,j,regl2-k,ic)
      endif
      end

      subroutine fccpy2(dest, src, nsets, dstart, sstart,
     @                  dstrd1, dstrd2, sstrd1, sstrd2,
     @                  nvals1, nvals2)
      real*8 dest(0:1), src(0:1)
      integer nsets
      integer dstart(nsets), sstart(nsets)
      integer dstrd1(nsets), dstrd2(nsets)
      integer sstrd1(nsets), sstrd2(nsets)
      integer nvals1(nsets), nvals2(nsets)
      integer i, j, k
      do 10 k = 1, nsets
         do 10 j = 0, nvals2(k) - 1
            do 10 i = 0, nvals1(k) - 1
 10            dest(dstart(k) + j * dstrd2(k) + i * dstrd1(k)) =
     @          src(sstart(k) + j * sstrd2(k) + i * sstrd1(k))
      end

c Works for CELL- or NODE-based data.
      subroutine mxnorm(src,
     @ dims(src),
     @ dims(reg),
     @ retval)
      integer dims(src)
      integer dims(reg)
      real*8 src(dimdec(src))
      real*8 retval
      integer bufsiz
      parameter(bufsiz = 128)
      real*8 buf(0:bufsiz-1)
      integer i, j, k, ii, jj, kk, idif, jdif, kdif, nbuf
      idif = regh0 - regl0
      jdif = regh1 - regl1
      kdif = regh2 - regl2
      if (kdif .gt. jdif .and. kdif .gt. idif) then
         nbuf = kdif - mod(kdif, bufsiz)
         do 10 k = 0, min(bufsiz - 1, kdif)
 10         buf(k) = 0.d0
         do 20 j = regl1, regh1
            do 20 i = regl0, regh0
               do 20 kk = 0, nbuf, bufsiz
                  do 20 k = 0, min(bufsiz - 1, kdif - kk)
 20                  buf(k) = max(buf(k), abs(src(i, j, regl2+kk+k)))
         do 30 k = 0, min(bufsiz - 1, kdif)
 30         retval = max(retval, buf(k))
      else if (jdif .gt. idif) then
         nbuf = jdif - mod(jdif, bufsiz)
         do 40 j = 0, min(bufsiz - 1, jdif)
 40         buf(j) = 0.d0
         do 50 k = regl2, regh2
            do 50 i = regl0, regh0
               do 50 jj = 0, nbuf, bufsiz
                  do 50 j = 0, min(bufsiz - 1, jdif - jj)
 50                  buf(j) = max(buf(j), abs(src(i, regl1+jj+j, k)))
         do 60 j = 0, min(bufsiz - 1, jdif)
 60         retval = max(retval, buf(j))
      else
         nbuf = idif - mod(idif, bufsiz)
         do 70 i = 0, min(bufsiz - 1, idif)
 70         buf(i) = 0.d0
         do 80 k = regl2, regh2
            do 80 j = regl1, regh1
               do 80 ii = 0, nbuf, bufsiz
                  do 80 i = 0, min(bufsiz - 1, idif - ii)
 80                  buf(i) = max(buf(i), abs(src(regl0+ii+i, j, k)))
         do 90 i = 0, min(bufsiz - 1, idif)
 90         retval = max(retval, buf(i))
      endif
      end

c Works for CELL- or NODE-based data.
      subroutine fassnu(v, x,
     @ dims(reg),
     @ ncomp)
      integer dims(reg)
      real*8 v(*)
      real*8 x
      integer ncomp
      integer nvals, i
      nvals = (regh0-regl0+1) * (regh1-regl1+1) * (regh2-regl2+1) * ncomp
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
      do 10 i = 1, (regh0-regl0+1) * (regh1-regl1+1) * (regh2-regl2+1)
 10      v0(i) = v0(i) + v1(i)
      end

c Works for CELL-based data.
      subroutine iprodc(v0,
     @ dims(v0),
     @ v1, dims(v1),
     @ dims(reg), sum)
      integer dims(v0)
      integer dims(v1)
      integer dims(reg)
      real*8 v0(dimdec(v0))
      real*8 v1(dimdec(v1))
      real*8 sum
      integer i, j, k
      if (regh2 - regl2 .gt. regh1 - regl1 .and.
     @    regh2 - regl2 .gt. regh0 - regl0) then
         do 10 j = regl1, regh1
            do 10 i = regl0, regh0
               do 10 k = regl2, regh2
 10               sum = sum + v0(i,j,k) * v1(i,j,k)
      else if (regh1 - regl1 .gt. regh0 - regl0) then
         do 20 k = regl2, regh2
            do 20 i = regl0, regh0
               do 20 j = regl1, regh1
 20               sum = sum + v0(i,j,k) * v1(i,j,k)
      else
         do 30 k = regl2, regh2
            do 30 j = regl1, regh1
               do 30 i = regl0, regh0
 30               sum = sum + v0(i,j,k) * v1(i,j,k)
      endif
      end

c Works for CELL- or NODE-based data.
      subroutine iprodn(v0,
     @ dims(v0),
     @ v1, dims(v1),
     @ dims(reg), sum)
      integer dims(v0)
      integer dims(v1)
      integer dims(reg)
      real*8 v0(dimdec(v0))
      real*8 v1(dimdec(v1))
      real*8 sum
      integer i, j, k
      sum = sum + 0.125d0 *
     @            (v0(regl0,regl1,regl2) * v1(regl0,regl1,regl2) +
     @             v0(regl0,regl1,regh2) * v1(regl0,regl1,regh2) +
     @             v0(regl0,regh1,regl2) * v1(regl0,regh1,regl2) +
     @             v0(regl0,regh1,regh2) * v1(regl0,regh1,regh2) +
     @             v0(regh0,regl1,regl2) * v1(regh0,regl1,regl2) +
     @             v0(regh0,regl1,regh2) * v1(regh0,regl1,regh2) +
     @             v0(regh0,regh1,regl2) * v1(regh0,regh1,regl2) +
     @             v0(regh0,regh1,regh2) * v1(regh0,regh1,regh2))
      do 10 i = regl0 + 1, regh0 - 1
 10      sum = sum + 0.25d0 *
     @               (v0(i,regl1,regl2) * v1(i,regl1,regl2) +
     @                v0(i,regl1,regh2) * v1(i,regl1,regh2) +
     @                v0(i,regh1,regl2) * v1(i,regh1,regl2) +
     @                v0(i,regh1,regh2) * v1(i,regh1,regh2))
      do 20 j = regl1 + 1, regh1 - 1
 20      sum = sum + 0.25d0 *
     @               (v0(regl0,j,regl2) * v1(regl0,j,regl2) +
     @                v0(regl0,j,regh2) * v1(regl0,j,regh2) +
     @                v0(regh0,j,regl2) * v1(regh0,j,regl2) +
     @                v0(regh0,j,regh2) * v1(regh0,j,regh2))
      do 30 k = regl2 + 1, regh2 - 1
 30      sum = sum + 0.25d0 *
     @               (v0(regl0,regl1,k) * v1(regl0,regl1,k) +
     @                v0(regl0,regh1,k) * v1(regl0,regh1,k) +
     @                v0(regh0,regl1,k) * v1(regh0,regl1,k) +
     @                v0(regh0,regh1,k) * v1(regh0,regh1,k))
      do 40 j = regl1 + 1, regh1 - 1
         do 40 i = regl0 + 1, regh0 - 1
 40         sum = sum + 0.5d0 *
     @                  (v0(i,j,regl2) * v1(i,j,regl2) +
     @                   v0(i,j,regh2) * v1(i,j,regh2))
      do 50 k = regl2 + 1, regh2 - 1
         do 50 i = regl0 + 1, regh0 - 1
 50         sum = sum + 0.5d0 *
     @                  (v0(i,regl1,k) * v1(i,regl1,k) +
     @                   v0(i,regh1,k) * v1(i,regh1,k))
      do 60 k = regl2 + 1, regh2 - 1
         do 60 j = regl1 + 1, regh1 - 1
 60         sum = sum + 0.5d0 *
     @                  (v0(regl0,j,k) * v1(regl0,j,k) +
     @                   v0(regh0,j,k) * v1(regh0,j,k))
      if (regh2 - regl2 .gt. regh1 - regl1 .and.
     @    regh2 - regl2 .gt. regh0 - regl0) then
         do 70 j = regl1 + 1, regh1 - 1
            do 70 i = regl0 + 1, regh0 - 1
               do 70 k = regl2 + 1, regh2 - 1
 70               sum = sum + v0(i,j,k) * v1(i,j,k)
      else if (regh1 - regl1 .gt. regh0 - regl0) then
         do 80 k = regl2 + 1, regh2 - 1
            do 80 i = regl0 + 1, regh0 - 1
               do 80 j = regl1 + 1, regh1 - 1
 80               sum = sum + v0(i,j,k) * v1(i,j,k)
      else
         do 90 k = regl2 + 1, regh2 - 1
            do 90 j = regl1 + 1, regh1 - 1
               do 90 i = regl0 + 1, regh0 - 1
 90               sum = sum + v0(i,j,k) * v1(i,j,k)
      endif
      end

c Works for CELL- or NODE-based data.
      subroutine bref(dest,
     @ dims(dest),
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
      integer i, j, k
      if (idir .eq. 0) then
         do 10 i = regl0, regh0
            do 10 k = regl2, regh2
               do 10 j = regl1, regh1
 10               dest(i,j,k) = src(bbh0-(i-regl0),j,k)
      else if (idir .eq. 1) then
         do 20 j = regl1, regh1
            do 20 k = regl2, regh2
               do 20 i = regl0, regh0
 20               dest(i,j,k) = src(i,bbh1-(j-regl1),k)
      else
         do 30 k = regl2, regh2
            do 30 j = regl1, regh1
               do 30 i = regl0, regh0
 30               dest(i,j,k) = src(i,j,bbh2-(k-regl2))
      endif
      end

c Works for CELL- or NODE-based data.
      subroutine brefm(dest,
     @ dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ dims(bb), ra)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer dims(bb)
      integer ra(0:2)
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      integer i, j, k
      if (ra(0) .eq. 0 .and. ra(1) .eq. 0 .and. ra(2) .eq. 0) then
         if (regh2 - regl2 .gt. max(regh0 - regl0, regh1 - regl1)) then
            do 10 j = regl1, regh1
               do 10 i = regl0, regh0
                  do 10 k = regl2, regh2
 10                  dest(i,j,k) = src(bbl0+(i-regl0),
     @                                 bbl1+(j-regl1), bbl2+(k-regl2))
         else if (regh1 - regl1 .gt. regh0 - regl0) then
            do 12 k = regl2, regh2
               do 12 i = regl0, regh0
                  do 12 j = regl1, regh1
 12                  dest(i,j,k) = src(bbl0+(i-regl0),
     @                                 bbl1+(j-regl1), bbl2+(k-regl2))
         else
            do 14 k = regl2, regh2
               do 14 j = regl1, regh1
                  do 14 i = regl0, regh0
 14                  dest(i,j,k) = src(bbl0+(i-regl0),
     @                                 bbl1+(j-regl1), bbl2+(k-regl2))
         endif
      else if (ra(0) .eq. 0 .and. ra(1) .eq. 0 .and. ra(2) .eq. 1) then
         do 20 k = regl2, regh2
            do 20 j = regl1, regh1
               do 20 i = regl0, regh0
 20               dest(i,j,k) = src(bbl0+(i-regl0),
     @                              bbl1+(j-regl1), bbh2-(k-regl2))
      else if (ra(0) .eq. 0 .and. ra(1) .eq. 1 .and. ra(2) .eq. 0) then
         do 30 j = regl1, regh1
            do 30 k = regl2, regh2
               do 30 i = regl0, regh0
 30               dest(i,j,k) = src(bbl0+(i-regl0),
     @                              bbh1-(j-regl1), bbl2+(k-regl2))
      else if (ra(0) .eq. 0 .and. ra(1) .eq. 1 .and. ra(2) .eq. 1) then
         do 40 k = regl2, regh2
            do 40 j = regl1, regh1
               do 40 i = regl0, regh0
 40               dest(i,j,k) = src(bbl0+(i-regl0),
     @                              bbh1-(j-regl1), bbh2-(k-regl2))
      else if (ra(0) .eq. 1 .and. ra(1) .eq. 0 .and. ra(2) .eq. 0) then
         do 50 i = regl0, regh0
            do 50 k = regl2, regh2
               do 50 j = regl1, regh1
 50               dest(i,j,k) = src(bbh0-(i-regl0),
     @                              bbl1+(j-regl1), bbl2+(k-regl2))
      else if (ra(0) .eq. 1 .and. ra(1) .eq. 0 .and. ra(2) .eq. 1) then
         do 60 k = regl2, regh2
            do 60 i = regl0, regh0
               do 60 j = regl1, regh1
 60               dest(i,j,k) = src(bbh0-(i-regl0),
     @                              bbl1+(j-regl1), bbh2-(k-regl2))
      else if (ra(0) .eq. 1 .and. ra(1) .eq. 1 .and. ra(2) .eq. 0) then
         do 70 j = regl1, regh1
            do 70 i = regl0, regh0
               do 70 k = regl2, regh2
 70               dest(i,j,k) = src(bbh0-(i-regl0),
     @                              bbh1-(j-regl1), bbl2+(k-regl2))
      else
         do 80 k = regl2, regh2
            do 80 j = regl1, regh1
               do 80 i = regl0, regh0
 80               dest(i,j,k) = src(bbh0-(i-regl0),
     @                              bbh1-(j-regl1), bbh2-(k-regl2))
      endif
      end

c Works for CELL- or NODE-based data.
      subroutine bneg(dest,
     @ dims(dest),
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
      integer i, j, k
      if (idir .eq. 0) then
         do 10 i = regl0, regh0
            do 10 k = regl2, regh2
               do 10 j = regl1, regh1
 10               dest(i,j,k) = -src(bbh0-(i-regl0),j,k)
      else if (idir .eq. 1) then
         do 20 j = regl1, regh1
            do 20 k = regl2, regh2
               do 20 i = regl0, regh0
 20               dest(i,j,k) = -src(i,bbh1-(j-regl1),k)
      else
         do 30 k = regl2, regh2
            do 30 j = regl1, regh1
               do 30 i = regl0, regh0
 30               dest(i,j,k) = -src(i,j,bbh2-(k-regl2))
      endif
      end

c Works for CELL- or NODE-based data.
      subroutine bnegm(dest,
     @ dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ dims(bb), ra)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer dims(bb)
      integer ra(0:2)
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      integer i, j, k
      if (ra(0) .eq. 0 .and. ra(1) .eq. 0 .and. ra(2) .eq. 0) then
         if (regh2 - regl2 .gt. max(regh0 - regl0, regh1 - regl1)) then
            do 10 j = regl1, regh1
               do 10 i = regl0, regh0
                  do 10 k = regl2, regh2
 10                  dest(i,j,k) = -src(bbl0+(i-regl0),
     @                                  bbl1+(j-regl1), bbl2+(k-regl2))
         else if (regh1 - regl1 .gt. regh0 - regl0) then
            do 12 k = regl2, regh2
               do 12 i = regl0, regh0
                  do 12 j = regl1, regh1
 12                  dest(i,j,k) = -src(bbl0+(i-regl0),
     @                                  bbl1+(j-regl1), bbl2+(k-regl2))
         else
            do 14 k = regl2, regh2
               do 14 j = regl1, regh1
                  do 14 i = regl0, regh0
 14                  dest(i,j,k) = -src(bbl0+(i-regl0),
     @                                  bbl1+(j-regl1), bbl2+(k-regl2))
         endif
      else if (ra(0) .eq. 0 .and. ra(1) .eq. 0 .and. ra(2) .eq. 1) then
         do 20 k = regl2, regh2
            do 20 j = regl1, regh1
               do 20 i = regl0, regh0
 20               dest(i,j,k) = -src(bbl0+(i-regl0),
     @                               bbl1+(j-regl1), bbh2-(k-regl2))
      else if (ra(0) .eq. 0 .and. ra(1) .eq. 1 .and. ra(2) .eq. 0) then
         do 30 j = regl1, regh1
            do 30 k = regl2, regh2
               do 30 i = regl0, regh0
 30               dest(i,j,k) = -src(bbl0+(i-regl0),
     @                               bbh1-(j-regl1), bbl2+(k-regl2))
      else if (ra(0) .eq. 0 .and. ra(1) .eq. 1 .and. ra(2) .eq. 1) then
         do 40 k = regl2, regh2
            do 40 j = regl1, regh1
               do 40 i = regl0, regh0
 40               dest(i,j,k) = -src(bbl0+(i-regl0),
     @                               bbh1-(j-regl1), bbh2-(k-regl2))
      else if (ra(0) .eq. 1 .and. ra(1) .eq. 0 .and. ra(2) .eq. 0) then
         do 50 i = regl0, regh0
            do 50 k = regl2, regh2
               do 50 j = regl1, regh1
 50               dest(i,j,k) = -src(bbh0-(i-regl0),
     @                               bbl1+(j-regl1), bbl2+(k-regl2))
      else if (ra(0) .eq. 1 .and. ra(1) .eq. 0 .and. ra(2) .eq. 1) then
         do 60 k = regl2, regh2
            do 60 i = regl0, regh0
               do 60 j = regl1, regh1
 60               dest(i,j,k) = -src(bbh0-(i-regl0),
     @                               bbl1+(j-regl1), bbh2-(k-regl2))
      else if (ra(0) .eq. 1 .and. ra(1) .eq. 1 .and. ra(2) .eq. 0) then
         do 70 j = regl1, regh1
            do 70 i = regl0, regh0
               do 70 k = regl2, regh2
 70               dest(i,j,k) = -src(bbh0-(i-regl0),
     @                               bbh1-(j-regl1), bbl2+(k-regl2))
      else
         do 80 k = regl2, regh2
            do 80 j = regl1, regh1
               do 80 i = regl0, regh0
 80               dest(i,j,k) = -src(bbh0-(i-regl0),
     @                               bbh1-(j-regl1), bbh2-(k-regl2))
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
      integer i, j, k
      if (idir .eq. 0) then
         if (regl0 .lt. bbh0) then
            do 10 i = regl0, regh0
               do 10 k = regl2, regh2
                  do 10 j = regl1, regh1
 10                  dest(i,j,k) = 2.d0 * dest(regh0,j,k) -
     @                             src(bbh0-(i-regl0),j,k)
         else
            do 20 i = regh0, regl0, -1
               do 20 k = regl2, regh2
                  do 20 j = regl1, regh1
 20                  dest(i,j,k) = 2.d0 * dest(regh0,j,k) -
     @                             src(bbh0-(i-regl0),j,k)
         endif
      else if (idir .eq. 1) then
         if (regl1 .lt. bbh1) then
            do 30 j = regl1, regh1
               do 30 k = regl2, regh2
                  do 30 i = regl0, regh0
 30                  dest(i,j,k) = 2.d0 * dest(i,regh1,k) -
     @                             src(i,bbh1-(j-regl1),k)
         else
            do 40 j = regh1, regl1, -1
               do 40 k = regl2, regh2
                  do 40 i = regl0, regh0
 40                  dest(i,j,k) = 2.d0 * dest(i,regh1,k) -
     @                             src(i,bbh1-(j-regl1),k)
         endif
      else
         if (regl2 .lt. bbh2) then
            do 50 k = regl2, regh2
               do 50 j = regl1, regh1
                  do 50 i = regl0, regh0
 50                  dest(i,j,k) = 2.d0 * dest(i,j,regh2) -
     @                             src(i,j,bbh2-(k-regl2))
         else
            do 60 k = regh2, regl2, -1
               do 60 j = regl1, regh1
                  do 60 i = regl0, regh0
 60                  dest(i,j,k) = 2.d0 * dest(i,j,regh2) -
     @                             src(i,j,bbh2-(k-regl2))
         endif
      endif
      end

c Works for CELL-based velocity data.
c This routine is called when the inflow face velocity data has already
c been altered by a call to fill_borders.
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
      integer i, j, k
      if (idir .eq. 0) then
         if (regl0 .lt. bbh0) then
            do 10 i = regl0, regh0
               do 10 k = regl2, regh2
                  do 10 j = regl1, regh1
 10                  dest(i,j,k) = src(bbl0,j,k) + src(bbl0+1,j,k) -
     @                             src(bbh0-(i-regl0),j,k)
         else
            do 20 i = regh0, regl0, -1
               do 20 k = regl2, regh2
                  do 20 j = regl1, regh1
 20                  dest(i,j,k) = src(bbh0-1,j,k) + src(bbh0,j,k) -
     @                             src(bbl0+(regh0-i),j,k)
         endif
      else if (idir .eq. 1) then
         if (regl1 .lt. bbh1) then
            do 30 j = regl1, regh1
               do 30 k = regl2, regh2
                  do 30 i = regl0, regh0
 30                  dest(i,j,k) = src(i,bbl1,k) + src(i,bbl1+1,k) -
     @                             src(i,bbh1-(j-regl1),k)
         else
            do 40 j = regh1, regl1, -1
               do 40 k = regl2, regh2
                  do 40 i = regl0, regh0
 40                  dest(i,j,k) = src(i,bbh1-1,k) + src(i,bbh1,k) -
     @                             src(i,bbl1+(regh1-j),k)
         endif
      else
         if (regl2 .lt. bbh2) then
            do 50 k = regl2, regh2
               do 50 j = regl1, regh1
                  do 50 i = regl0, regh0
 50                  dest(i,j,k) = src(i,j,bbl2) + src(i,j,bbl2+1) -
     @                             src(i,j,bbh2-(k-regl2))
         else
            do 60 k = regh2, regl2, -1
               do 60 j = regl1, regh1
                  do 60 i = regl0, regh0
 60                  dest(i,j,k) = src(i,j,bbh2-1) + src(i,j,bbh2) -
     @                             src(i,j,bbl2+(regh2-k))
         endif
      endif
      end

c CELL-based data only.
      subroutine acint2(dest,
     @ dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ dims(bb), ir, jr, kr)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer dims(bb)
      integer ir, jr, kr
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      real*8 xoff, yoff, zoff, sx, sy, sz
      integer ic, jc, kc, i, j, k
      do 10 k = regl2, regh2
         kc = k/kr
         zoff = (mod(k,kr) + 0.5) / kr - 0.5
         do 10 j = regl1, regh1
            jc = j/jr
            yoff = (mod(j,jr) + 0.5) / jr - 0.5
            do 10 i = regl0, regh0
               ic = i/ir
               xoff = (mod(i,ir) + 0.5) / ir - 0.5
               sz = 0.5 * (src(ic,jc,kc+1) - src(ic,jc,kc-1))
               sy = 0.5 * (src(ic,jc+1,kc) - src(ic,jc-1,kc))
               sx = 0.5 * (src(ic+1,jc,kc) - src(ic-1,jc,kc))
 10         dest(i,j,k) = src(ic,jc,kc) + xoff*sx + yoff*sy + zoff*sz
      end

c NODE-based data only.
      subroutine anint2(dest,
     @ dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ dims(bb), ir, jr, kr)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer dims(bb)
      integer ir, jr, kr
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      real*8 p, q
      integer ic, jc, kc, j, k, m
      do 10 kc = bbl2, bbh2
         do 10 jc = bbl1, bbh1
            do 10 ic = bbl0, bbh0
 10            dest(ir*ic,jr*jc,kr*kc) = src(ic,jc,kc)
      if (kr .eq. 2) then
         do 20 kc = bbl2, bbh2-1
            do 20 jc = bbl1, bbh1
               do 20 ic = bbl0, bbh0
 20               dest(ir*ic,jr*jc,2*kc+1) = 0.5d0 * (src(ic,jc,kc) +
     @                                                src(ic,jc,kc+1))
      else if (kr .eq. 4) then
         do 120 kc = bbl2, bbh2-1
            do 120 jc = bbl1, bbh1
cdir$ ivdep
               do 120 ic = bbl0, bbh0
                  dest(ir*ic,jr*jc,4*kc+2) = 0.5d0 *
     @               (src(ic,jc,kc) + src(ic,jc,kc+1))
                  dest(ir*ic,jr*jc,4*kc+1) = 0.5d0 *
     @               (dest(ir*ic,jr*jc,4*kc) + dest(ir*ic,jr*jc,4*kc+2))
 120              dest(ir*ic,jr*jc,4*kc+3) = 0.5d0 *
     @               (dest(ir*ic,jr*jc,4*kc+2) + dest(ir*ic,jr*jc,4*kc+4))
      else if (kr .ne. 1) then
         do 220 m = 1, kr-1
            q = dble(m)/kr
            p = 1.0 - q
            do 220 kc = bbl2, bbh2-1
               do 220 jc = bbl1, bbh1
                  do 220 ic = bbl0, bbh0
 220                 dest(ir*ic,jr*jc,kr*kc+m) =
     @                  p * src(ic,jc,kc) + q * src(ic,jc,kc+1)
      endif
      if (jr .eq. 2) then
         do 30 jc = bbl1, bbh1-1
            do 30 ic = bbl0, bbh0
cdir$ ivdep
               do 30 k = regl2, regh2
 30               dest(ir*ic,2*jc+1,k) = 0.5d0 *
     @               (dest(ir*ic,2*jc,k) + dest(ir*ic,2*jc+2,k))
      else if (jr .eq. 4) then
         do 130 jc = bbl1, bbh1-1
            do 130 ic = bbl0, bbh0
cdir$ ivdep
               do 130 k = regl2, regh2
                  dest(ir*ic,4*jc+2,k) = 0.5d0 *
     @               (dest(ir*ic,4*jc,k) + dest(ir*ic,4*jc+4,k))
                  dest(ir*ic,4*jc+1,k) = 0.5d0 *
     @               (dest(ir*ic,4*jc,k) + dest(ir*ic,4*jc+2,k))
 130              dest(ir*ic,4*jc+3,k) = 0.5d0 *
     @               (dest(ir*ic,4*jc+2,k) + dest(ir*ic,4*jc+4,k))
      else if (jr .ne. 1) then
         do 230 m = 1, jr-1
            q = dble(m)/jr
            p = 1.0 - q
            do 230 jc = bbl1, bbh1-1
               do 230 ic = bbl0, bbh0
cdir$ ivdep
                  do 230 k = regl2, regh2
 230                 dest(ir*ic,jr*jc+m,k) =
     @                  p * dest(ir*ic,jr*jc,k) +
     @                  q * dest(ir*ic,jr*(jc+1),k)
      endif
      if (ir .eq. 2) then
         if (regh2 - regl2 .gt. regh1 - regl1) then
            do 40 ic = bbl0, bbh0-1
               do 40 j = regl1, regh1
cdir$ ivdep
                  do 40 k = regl2, regh2
 40                  dest(2*ic+1,j,k) = 0.5d0 *
     @                  (dest(2*ic,j,k) + dest(2*ic+2,j,k))
         else
            do 50 ic = bbl0, bbh0-1
               do 50 k = regl2, regh2
cdir$ ivdep
                  do 50 j = regl1, regh1
 50                  dest(2*ic+1,j,k) = 0.5d0 *
     @                  (dest(2*ic,j,k) + dest(2*ic+2,j,k))
         endif
      else if (ir .eq. 4) then
         if (regh2 - regl2 .gt. regh1 - regl1) then
            do 140 ic = bbl0, bbh0-1
               do 140 j = regl1, regh1
cdir$ ivdep
                  do 140 k = regl2, regh2
                     dest(4*ic+2,j,k) = 0.5d0 *
     @                  (dest(4*ic,j,k) + dest(4*ic+4,j,k))
                     dest(4*ic+1,j,k) = 0.5d0 *
     @                  (dest(4*ic,j,k) + dest(4*ic+2,j,k))
 140                 dest(4*ic+3,j,k) = 0.5d0 *
     @                  (dest(4*ic+2,j,k) + dest(4*ic+4,j,k))
         else
            do 150 ic = bbl0, bbh0-1
               do 150 k = regl2, regh2
cdir$ ivdep
                  do 150 j = regl1, regh1
                     dest(4*ic+2,j,k) = 0.5d0 *
     @                  (dest(4*ic,j,k) + dest(4*ic+4,j,k))
                     dest(4*ic+1,j,k) = 0.5d0 *
     @                  (dest(4*ic,j,k) + dest(4*ic+2,j,k))
 150                 dest(4*ic+3,j,k) = 0.5d0 *
     @                  (dest(4*ic+2,j,k) + dest(4*ic+4,j,k))
         endif
      else if (ir .ne. 1) then
         if (regh2 - regl2 .gt. regh1 - regl1) then
            do 240 m = 1, ir-1
               q = dble(m)/ir
               p = 1.0 - q
               do 240 ic = bbl0, bbh0-1
                  do 240 k = regl2, regh2
cdir$ ivdep
                     do 240 j = regl1, regh1
 240                    dest(ir*ic+m,j,k) =
     @                     p * dest(ir*ic,j,k) +
     @                     q * dest(ir*(ic+1),j,k)
         else
            do 250 m = 1, ir-1
               q = dble(m)/ir
               p = 1.0 - q
               do 250 ic = bbl0, bbh0-1
                  do 250 k = regl2, regh2
cdir$ ivdep
                     do 250 j = regl1, regh1
 250                    dest(ir*ic+m,j,k) =
     @                     p * dest(ir*ic,j,k) +
     @                     q * dest(ir*(ic+1),j,k)
         endif
      endif
      end

c CELL-based data only.
      subroutine acrst1(dest,
     @ dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ ir, jr, kr)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer ir, jr, kr
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      real*8 fac
      integer i, j, k, l, m, n
      do 10 k = regl2, regh2
         do 10 j = regl1, regh1
            do 10 i = regl0, regh0
 10            dest(i,j,k) = 0.d0
      do 20 l = 0, kr-1
         do 20 n = 0, jr-1
            do 20 m = 0, ir-1
               do 20 k = regl2, regh2
                  do 20 j = regl1, regh1
                     do 20 i = regl0, regh0
 20                     dest(i,j,k) = dest(i,j,k) +
     @                    src(i*ir+m, j*jr+n, k*kr+l)
      fac = 1.d0 / (ir*jr*kr)
      do 30 k = regl2, regh2
         do 30 j = regl1, regh1
            do 30 i = regl0, regh0
 30            dest(i,j,k) = dest(i,j,k) * fac
      end

c NODE-based data only.
      subroutine anrst1(dest,
     @ dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ ir, jr, kr)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer ir, jr, kr
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      integer i, j, k
      do 10 k = regl2, regh2
         do 10 j = regl1, regh1
            do 10 i = regl0, regh0
 10            dest(i,j,k) = src(i*ir, j*jr, k*kr)
      end

#define RESTR2 0.125d0 * src(2*i,2*j,2*k) + \
       0.0625d0   * (src(2*i,2*j,2*k-1) + src(2*i,2*j,2*k+1) + \
                     src(2*i,2*j-1,2*k) + src(2*i,2*j+1,2*k) + \
                     src(2*i-1,2*j,2*k) + src(2*i+1,2*j,2*k)) + \
       0.03125d0  * (src(2*i,2*j-1,2*k-1) + src(2*i,2*j-1,2*k+1) + \
                     src(2*i,2*j+1,2*k-1) + src(2*i,2*j+1,2*k+1) + \
                     src(2*i-1,2*j,2*k-1) + src(2*i-1,2*j,2*k+1) + \
                     src(2*i+1,2*j,2*k-1) + src(2*i+1,2*j,2*k+1) + \
                     src(2*i-1,2*j-1,2*k) + src(2*i-1,2*j+1,2*k) + \
                     src(2*i+1,2*j-1,2*k) + src(2*i+1,2*j+1,2*k)) + \
       0.015625d0 * (src(2*i-1,2*j-1,2*k-1) + src(2*i-1,2*j-1,2*k+1) + \
                     src(2*i-1,2*j+1,2*k-1) + src(2*i-1,2*j+1,2*k+1) + \
                     src(2*i+1,2*j-1,2*k-1) + src(2*i+1,2*j-1,2*k+1) + \
                     src(2*i+1,2*j+1,2*k-1) + src(2*i+1,2*j+1,2*k+1))

#define RR0(r) (src(4*i-r,4*j-r,4*k) + src(4*i-r,4*j+r,4*k) + \
               src(4*i+r,4*j-r,4*k) + src(4*i+r,4*j+r,4*k) + \
               src(4*i-r,4*j,4*k-r) + src(4*i-r,4*j,4*k+r) + \
               src(4*i+r,4*j,4*k-r) + src(4*i+r,4*j,4*k+r) + \
               src(4*i,4*j-r,4*k-r) + src(4*i,4*j-r,4*k+r) + \
               src(4*i,4*j+r,4*k-r) + src(4*i,4*j+r,4*k+r))
#define RS0(r,s) (src(4*i-r,4*j-s,4*k) + src(4*i-r,4*j+s,4*k) + \
                 src(4*i+r,4*j-s,4*k) + src(4*i+r,4*j+s,4*k) + \
                 src(4*i-r,4*j,4*k-s) + src(4*i-r,4*j,4*k+s) + \
                 src(4*i+r,4*j,4*k-s) + src(4*i+r,4*j,4*k+s) + \
                 src(4*i,4*j-r,4*k-s) + src(4*i,4*j-r,4*k+s) + \
                 src(4*i,4*j+r,4*k-s) + src(4*i,4*j+r,4*k+s))
#define REST4A  0.015625d0 * src(4*i,4*j,4*k) + \
       0.01171875d0 * (src(4*i,4*j,4*k-1) + src(4*i,4*j,4*k+1) + \
                       src(4*i,4*j-1,4*k) + src(4*i,4*j+1,4*k) + \
                       src(4*i-1,4*j,4*k) + src(4*i+1,4*j,4*k)) + \
       0.0078125d0  * (src(4*i,4*j,4*k-2) + src(4*i,4*j,4*k+2) + \
                       src(4*i,4*j-2,4*k) + src(4*i,4*j+2,4*k) + \
                       src(4*i-2,4*j,4*k) + src(4*i+2,4*j,4*k)) + \
       0.00390625d0 * (src(4*i,4*j,4*k-3) + src(4*i,4*j,4*k+3) + \
                       src(4*i,4*j-3,4*k) + src(4*i,4*j+3,4*k) + \
                       src(4*i-3,4*j,4*k) + src(4*i+3,4*j,4*k))
#define REST4B  0.0087890625d0 * RR0(1) + \
       0.00390625d0 * RR0(2) + 0.0009765625d0 * RR0(3)
#define REST4C  0.001953125d0 * (RS0(2,3) + RS0(3,2))
#define REST4D  0.005859375d0 * (RS0(1,2) + RS0(2,1))
#define REST4E  0.0029296875d0 * (RS0(1,3) + RS0(3,1))

c NODE-based data only.
      subroutine anrst2(dest,
     @ dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ ir, jr, kr, integ)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer ir, jr, kr, integ
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      real*8 fac0, fac1, fac2, fac
      integer i, j, k, l, m, n
      if (ir .eq. 2 .and. jr .eq. 2 .and. kr .eq. 2) then
         if (regh2 - regl2 .gt. regh1 - regl1 .and.
     @       regh2 - regl2 .gt. regh0 - regl0) then
            do 10 j = regl1, regh1
               do 10 i = regl0, regh0
                  do 10 k = regl2, regh2
 10                  dest(i,j,k) = RESTR2
         else if (regh1 - regl1 .gt. regh0 - regl0) then
            do 20 k = regl2, regh2
               do 20 i = regl0, regh0
                  do 20 j = regl1, regh1
 20                  dest(i,j,k) = RESTR2
         else
            do 30 k = regl2, regh2
               do 30 j = regl1, regh1
                  do 30 i = regl0, regh0
 30                  dest(i,j,k) = RESTR2
         endif
      else if (ir .eq. 4 .and. jr .eq. 4 .and. kr .eq. 4) then
         if (regh2 - regl2 .gt. regh1 - regl1 .and.
     @       regh2 - regl2 .gt. regh0 - regl0) then
            do 110 j = regl1, regh1
               do 110 i = regl0, regh0
                  do 101 k = regl2, regh2
 101                 dest(i,j,k) = REST4A
                  do 102 k = regl2, regh2
 102                 dest(i,j,k) = dest(i,j,k) + REST4B
                  do 103 k = regl2, regh2
 103                 dest(i,j,k) = dest(i,j,k) + REST4C
                  do 104 k = regl2, regh2
 104                 dest(i,j,k) = dest(i,j,k) + REST4D
                  do 110 k = regl2, regh2
 110                 dest(i,j,k) = dest(i,j,k) + REST4E
            fac0 = 1.d0 / (ir*ir * jr*jr * kr*kr)
            do 120 l = 1, 3
               fac2 = (4-l) * fac0
               do 120 n = 1, 3
                  fac1 = (4-n) * fac2
                  do 120 m = 1, 3 
                     fac = (4-m) * fac1
                     do 120 j = regl1, regh1
                        do 120 i = regl0, regh0
                           do 120 k = regl2, regh2
 120                          dest(i,j,k) = dest(i,j,k) +
     @                  fac * (src(4*i-m,4*j-n,4*k-l)+
     @                         src(4*i-m,4*j-n,4*k+l)+
     @                         src(4*i-m,4*j+n,4*k-l)+
     @                         src(4*i-m,4*j+n,4*k+l)+
     @                         src(4*i+m,4*j-n,4*k-l)+
     @                         src(4*i+m,4*j-n,4*k+l)+
     @                         src(4*i+m,4*j+n,4*k-l)+
     @                         src(4*i+m,4*j+n,4*k+l))
         else if (regh1 - regl1 .gt. regh0 - regl0) then
            do 130 k = regl2, regh2
               do 130 i = regl0, regh0
                  do 121 j = regl1, regh1
 121                 dest(i,j,k) = REST4A
                  do 122 j = regl1, regh1
 122                 dest(i,j,k) = dest(i,j,k) + REST4B
                  do 123 j = regl1, regh1
 123                 dest(i,j,k) = dest(i,j,k) + REST4C
                  do 124 j = regl1, regh1
 124                 dest(i,j,k) = dest(i,j,k) + REST4D
                  do 130 j = regl1, regh1
 130                 dest(i,j,k) = dest(i,j,k) + REST4E
            fac0 = 1.d0 / (ir*ir * jr*jr * kr*kr)
            do 140 l = 1, 3
               fac2 = (4-l) * fac0
               do 140 n = 1, 3
                  fac1 = (4-n) * fac2
                  do 140 m = 1, 3
                     fac = (4-m) * fac1
                     do 140 k = regl2, regh2
                        do 140 i = regl0, regh0
                           do 140 j = regl1, regh1
 140                          dest(i,j,k) = dest(i,j,k) +
     @                  fac * (src(4*i-m,4*j-n,4*k-l)+
     @                         src(4*i-m,4*j-n,4*k+l)+
     @                         src(4*i-m,4*j+n,4*k-l)+
     @                         src(4*i-m,4*j+n,4*k+l)+
     @                         src(4*i+m,4*j-n,4*k-l)+
     @                         src(4*i+m,4*j-n,4*k+l)+
     @                         src(4*i+m,4*j+n,4*k-l)+
     @                         src(4*i+m,4*j+n,4*k+l))
         else
            do 150 k = regl2, regh2
               do 150 j = regl1, regh1
                  do 141 i = regl0, regh0
 141                 dest(i,j,k) = REST4A
                  do 142 i = regl0, regh0
 142                 dest(i,j,k) = dest(i,j,k) + REST4B
                  do 143 i = regl0, regh0
 143                 dest(i,j,k) = dest(i,j,k) + REST4C
                  do 144 i = regl0, regh0
 144                 dest(i,j,k) = dest(i,j,k) + REST4D
                  do 150 i = regl0, regh0
 150                 dest(i,j,k) = dest(i,j,k) + REST4E
            fac0 = 1.d0 / (ir*ir * jr*jr * kr*kr)
            do 160 l = 1, 3
               fac2 = (4-l) * fac0
               do 160 n = 1, 3
                  fac1 = (4-n) * fac2
                  do 160 m = 1, 3
                     fac = (4-m) * fac1
                     do 160 k = regl2, regh2
                        do 160 j = regl1, regh1
                           do 160 i = regl0, regh0
 160                          dest(i,j,k) = dest(i,j,k) +
     @                  fac * (src(4*i-m,4*j-n,4*k-l)+
     @                         src(4*i-m,4*j-n,4*k+l)+
     @                         src(4*i-m,4*j+n,4*k-l)+
     @                         src(4*i-m,4*j+n,4*k+l)+
     @                         src(4*i+m,4*j-n,4*k-l)+
     @                         src(4*i+m,4*j-n,4*k+l)+
     @                         src(4*i+m,4*j+n,4*k-l)+
     @                         src(4*i+m,4*j+n,4*k+l))
         endif
      else
         if (regh2 - regl2 .gt. regh1 - regl1 .and.
     @       regh2 - regl2 .gt. regh0 - regl0) then
            do 210 j = regl1, regh1
               do 210 i = regl0, regh0
                  do 210 k = regl2, regh2
 210                 dest(i,j,k) = 0.d0
            fac0 = 1.d0 / (ir*ir * jr*jr * kr*kr)
            do 220 l = 0, kr-1
               fac2 = (kr-l) * fac0
               if (l .eq. 0) fac2 = 0.5d0 * fac2
               do 220 n = 0, jr-1
                  fac1 = (jr-n) * fac2
                  if (n .eq. 0) fac1 = 0.5d0 * fac1
                  do 220 m = 0, ir-1
                     fac = (ir-m) * fac1
                     if (m .eq. 0) fac = 0.5d0 * fac
                     do 220 j = regl1, regh1
                        do 220 i = regl0, regh0
                           do 220 k = regl2, regh2
 220                          dest(i,j,k) = dest(i,j,k) +
     @                  fac * (src(i*ir-m,j*jr-n,k*kr-l)+
     @                         src(i*ir-m,j*jr-n,k*kr+l)+
     @                         src(i*ir-m,j*jr+n,k*kr-l)+
     @                         src(i*ir-m,j*jr+n,k*kr+l)+
     @                         src(i*ir+m,j*jr-n,k*kr-l)+
     @                         src(i*ir+m,j*jr-n,k*kr+l)+
     @                         src(i*ir+m,j*jr+n,k*kr-l)+
     @                         src(i*ir+m,j*jr+n,k*kr+l))
         else if (regh1 - regl1 .gt. regh0 - regl0) then
            do 230 k = regl2, regh2
               do 230 i = regl0, regh0
                  do 230 j = regl1, regh1
 230                 dest(i,j,k) = 0.d0
            fac0 = 1.d0 / (ir*ir * jr*jr * kr*kr)
            do 240 l = 0, kr-1
               fac2 = (kr-l) * fac0
               if (l .eq. 0) fac2 = 0.5d0 * fac2
               do 240 n = 0, jr-1
                  fac1 = (jr-n) * fac2
                  if (n .eq. 0) fac1 = 0.5d0 * fac1
                  do 240 m = 0, ir-1
                     fac = (ir-m) * fac1
                     if (m .eq. 0) fac = 0.5d0 * fac
                     do 240 k = regl2, regh2
                        do 240 i = regl0, regh0
                           do 240 j = regl1, regh1
 240                          dest(i,j,k) = dest(i,j,k) +
     @                  fac * (src(i*ir-m,j*jr-n,k*kr-l)+
     @                         src(i*ir-m,j*jr-n,k*kr+l)+
     @                         src(i*ir-m,j*jr+n,k*kr-l)+
     @                         src(i*ir-m,j*jr+n,k*kr+l)+
     @                         src(i*ir+m,j*jr-n,k*kr-l)+
     @                         src(i*ir+m,j*jr-n,k*kr+l)+
     @                         src(i*ir+m,j*jr+n,k*kr-l)+
     @                         src(i*ir+m,j*jr+n,k*kr+l))
         else
            do 250 k = regl2, regh2
               do 250 j = regl1, regh1
                  do 250 i = regl0, regh0
 250                 dest(i,j,k) = 0.d0
            fac0 = 1.d0 / (ir*ir * jr*jr * kr*kr)
            do 260 l = 0, kr-1
               fac2 = (kr-l) * fac0
               if (l .eq. 0) fac2 = 0.5d0 * fac2
               do 260 n = 0, jr-1
                  fac1 = (jr-n) * fac2
                  if (n .eq. 0) fac1 = 0.5d0 * fac1
                  do 260 m = 0, ir-1
                     fac = (ir-m) * fac1
                     if (m .eq. 0) fac = 0.5d0 * fac
                     do 260 k = regl2, regh2
                        do 260 j = regl1, regh1
                           do 260 i = regl0, regh0
 260                          dest(i,j,k) = dest(i,j,k) +
     @                  fac * (src(i*ir-m,j*jr-n,k*kr-l)+
     @                         src(i*ir-m,j*jr-n,k*kr+l)+
     @                         src(i*ir-m,j*jr+n,k*kr-l)+
     @                         src(i*ir-m,j*jr+n,k*kr+l)+
     @                         src(i*ir+m,j*jr-n,k*kr-l)+
     @                         src(i*ir+m,j*jr-n,k*kr+l)+
     @                         src(i*ir+m,j*jr+n,k*kr-l)+
     @                         src(i*ir+m,j*jr+n,k*kr+l))
         endif
      endif
      if (integ .eq. 1) then
         fac = ir * jr * kr
         if (regh2 - regl2 .gt. regh1 - regl1 .and.
     @       regh2 - regl2 .gt. regh0 - regl0) then
            do 310 j = regl1, regh1
               do 310 i = regl0, regh0
                  do 310 k = regl2, regh2
 310                 dest(i,j,k) = fac * dest(i,j,k)
         else if (regh1 - regl1 .gt. regh0 - regl0) then
            do 320 k = regl2, regh2
               do 320 i = regl0, regh0
                  do 320 j = regl1, regh1
 320                 dest(i,j,k) = fac * dest(i,j,k)
         else
            do 330 k = regl2, regh2
               do 330 j = regl1, regh1
                  do 330 i = regl0, regh0
 330                 dest(i,j,k) = fac * dest(i,j,k)
         endif
      endif
      end

c NODE-based data only.
c Fills coarse region defined by reg.
c Handles coarse-fine face, with orientation determined by idim and idir
      subroutine anfr2(dest,
     @ dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ ir, jr, kr, idim, idir, integ)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer ir, jr, kr, idim, idir, integ
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      real*8 fac0, fac1, fac2, fac
      integer i, j, k, l, m, n
      if (idim .eq. 0) then
         if (integ .eq. 0) then
            fac = (0.5d0 + 0.5d0 / ir)
            fac0 = 1.d0 / (ir*ir * jr*jr * kr*kr)
         else
            fac = 1.d0
            fac0 = 1.d0 / (ir * jr * kr)
         endif
         i = regl0
         do 10 k = regl2, regh2
            do 10 j = regl1, regh1
 10            dest(i,j,k) = fac * src(i*ir,j*jr,k*kr)
         do 20 l = 0, kr-1
            fac2 = (kr-l) * fac0
            if (l .eq. 0) fac2 = 0.5d0 * fac2
            do 20 n = 0, jr-1
               fac1 = (jr-n) * fac2
               if (n .eq. 0) fac1 = 0.5d0 * fac1
               do 20 m = idir, idir*(ir-1), idir
                  fac = (ir-abs(m)) * fac1
                  do 20 k = regl2, regh2
                     do 20 j = regl1, regh1
 20                     dest(i,j,k) = dest(i,j,k) +
     @                       fac * (src(i*ir+m,j*jr-n,k*kr-l) +
     @                              src(i*ir+m,j*jr-n,k*kr+l) +
     @                              src(i*ir+m,j*jr+n,k*kr-l) +
     @                              src(i*ir+m,j*jr+n,k*kr+l))
      else if (idim .eq. 1) then
         if (integ .eq. 0) then
            fac = (0.5d0 + 0.5d0 / jr)
            fac0 = 1.d0 / (ir*ir * jr*jr * kr*kr)
         else
            fac = 1.d0
            fac0 = 1.d0 / (ir * jr * kr)
         endif
         j = regl1
         do 30 k = regl2, regh2
            do 30 i = regl0, regh0
 30            dest(i,j,k) = fac * src(i*ir,j*jr,k*kr)
         fac0 = 1.d0 / (ir*ir * jr*jr * kr*kr)
         do 40 l = 0, kr-1
            fac2 = (kr-l) * fac0
            if (l .eq. 0) fac2 = 0.5d0 * fac2
            do 40 n = idir, idir*(jr-1), idir
               fac1 = (jr-abs(n)) * fac2
               do 40 m = 0, ir-1
                  fac = (ir-m) * fac1
                  if (m .eq. 0) fac = 0.5d0 * fac
                  do 40 k = regl2, regh2
                     do 40 i = regl0, regh0
 40                     dest(i,j,k) = dest(i,j,k) +
     @                       fac * (src(i*ir-m,j*jr+n,k*kr-l) +
     @                              src(i*ir-m,j*jr+n,k*kr+l) +
     @                              src(i*ir+m,j*jr+n,k*kr-l) +
     @                              src(i*ir+m,j*jr+n,k*kr+l))
      else
         if (integ .eq. 0) then
            fac = (0.5d0 + 0.5d0 / kr)
            fac0 = 1.d0 / (ir*ir * jr*jr * kr*kr)
         else
            fac = 1.d0
            fac0 = 1.d0 / (ir * jr * kr)
         endif
         k = regl2
         do 50 j = regl1, regh1
            do 50 i = regl0, regh0
 50            dest(i,j,k) = fac * src(i*ir,j*jr,k*kr)
         fac0 = 1.d0 / (ir*ir * jr*jr * kr*kr)
         do 60 l = idir, idir*(kr-1), idir
            fac2 = (kr-abs(l)) * fac0
            do 60 n = 0, jr-1
               fac1 = (jr-n) * fac2
               if (n .eq. 0) fac1 = 0.5d0 * fac1
               do 60 m = 0, ir-1
                  fac = (ir-m) * fac1
                  if (m .eq. 0) fac = 0.5d0 * fac
                  do 60 j = regl1, regh1
                     do 60 i = regl0, regh0
 60                     dest(i,j,k) = dest(i,j,k) +
     @                       fac * (src(i*ir-m,j*jr-n,k*kr+l) +
     @                              src(i*ir-m,j*jr+n,k*kr+l) +
     @                              src(i*ir+m,j*jr-n,k*kr+l) +
     @                              src(i*ir+m,j*jr+n,k*kr+l))
      endif
      end

c NODE-based data only.
c Fills coarse region defined by reg.
c Handles any edge geometry except all-coarse or all-fine.
      subroutine aner2(dest,
     @ dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ ir, jr, kr, ivect, ga, integ)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer ir, jr, kr, ivect(0:2), ga(0:1,0:1,0:1), integ
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      real*8 cube, center, cfac, fac0, fac1, fac2, fac
      integer i, j, k, ii, ji, ki, idir, jdir, kdir, l, m, n
      cube = ir * jr * kr
      if (ivect(0) .eq. 0) then
         j = regl1
         k = regl2
         do 10 i = regl0, regh0
 10         dest(i,j,k) = 0.d0
c center gets center plus two edges
         if (integ .eq. 0) then
            center = ir / cube
            fac0 = 1.d0 / (cube**2)
            cfac = 0.25d0 * cube * fac0 * ir * (jr-1) * (kr-1)
         else
            center = 1.d0
            fac0 = 1.d0 / cube
         endif
c quadrants
c each quadrant is two octants and a face
         do 30 ki = 0, 1
            kdir = 2 * ki - 1
            do 30 ji = 0, 1
               jdir = 2 * ji - 1
               if (ga(0,ji,ki) .eq. 1) then
                  do 20 l = kdir, kdir*(kr-1), kdir
                     fac2 = (kr-abs(l)) * fac0
                     do 20 n = jdir, jdir*(jr-1), jdir
                        fac1 = (jr-abs(n)) * fac2
                        do 20 m = 0, ir-1
                           fac = (ir-m) * fac1
                           if (m .eq. 0) fac = 0.5d0 * fac
                           do 20 i = regl0, regh0
 20                           dest(i,j,k) = dest(i,j,k) +
     @                          fac * (src(i*ir-m,j*jr+n,k*kr+l) +
     @                                 src(i*ir+m,j*jr+n,k*kr+l))
               else if (integ .eq. 0) then
                  center = center + cfac
               endif
 30            continue
c faces
c each face is two faces and an edge
         fac2 = kr * fac0
         cfac = 0.5d0 * cube * fac0 * ir * (jr-1)
         do 50 ji = 0, 1
            jdir = 2 * ji - 1
            if (ga(0,ji,0) + ga(0,ji,1) .eq. 2) then
               do 40 n = jdir, jdir*(jr-1), jdir
                  fac1 = (jr-abs(n)) * fac2
                  do 40 m = 0, ir-1
                     fac = (ir-m) * fac1
                     if (m .eq. 0) fac = 0.5d0 * fac
                     do 40 i = regl0, regh0
 40                     dest(i,j,k) = dest(i,j,k) +
     @                    fac * (src(i*ir-m,j*jr+n,k*kr) +
     @                           src(i*ir+m,j*jr+n,k*kr))
            else if (integ .eq. 0) then
               center = center + cfac
            endif
 50         continue
         fac2 = jr * fac0
         cfac = 0.5d0 * cube * fac0 * ir * (kr-1)
         do 70 ki = 0, 1
            kdir = 2 * ki - 1
            if (ga(0,0,ki) + ga(0,1,ki) .eq. 2) then
               do 60 l = kdir, kdir*(kr-1), kdir
                  fac1 = (kr-abs(l)) * fac2
                  do 60 m = 0, ir-1
                     fac = (ir-m) * fac1
                     if (m .eq. 0) fac = 0.5d0 * fac
                     do 60 i = regl0, regh0
 60                     dest(i,j,k) = dest(i,j,k) +
     @                    fac * (src(i*ir-m,j*jr,k*kr+l) +
     @                           src(i*ir+m,j*jr,k*kr+l))
            else if (integ .eq. 0) then
               center = center + cfac
            endif
 70         continue
c center
         do 80 i = regl0, regh0
 80         dest(i,j,k) = dest(i,j,k) +
     @        center * src(i*ir,j*jr,k*kr)
      else if (ivect(1) .eq. 0) then
         i = regl0
         k = regl2
         do 110 j = regl1, regh1
 110        dest(i,j,k) = 0.d0
c center gets center plus two edges
         if (integ .eq. 0) then
            center = jr / cube
            fac0 = 1.d0 / (cube**2)
            cfac = 0.25d0 * cube * fac0 * jr * (ir-1) * (kr-1)
         else
            center = 1.d0
            fac0 = 1.d0 / cube
         endif
c quadrants
c each quadrant is two octants and a face
         do 130 ki = 0, 1
            kdir = 2 * ki - 1
            do 130 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,0,ki) .eq. 1) then
                  do 120 l = kdir, kdir*(kr-1), kdir
                     fac2 = (kr-abs(l)) * fac0
                     do 120 n = 0, jr-1
                        fac1 = (jr-n) * fac2
                        if (n .eq. 0) fac1 = 0.5d0 * fac1
                        do 120 m = idir, idir*(ir-1), idir
                           fac = (ir-abs(m)) * fac1
                           do 120 j = regl1, regh1
 120                          dest(i,j,k) = dest(i,j,k) +
     @                          fac * (src(i*ir+m,j*jr-n,k*kr+l) +
     @                                 src(i*ir+m,j*jr+n,k*kr+l))
               else if (integ .eq. 0) then
                  center = center + cfac
               endif
 130           continue
c faces
c each face is two faces and an edge
         fac2 = kr * fac0
         cfac = 0.5d0 * cube * fac0 * jr * (ir-1)
         do 150 ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,0,0) + ga(ii,0,1) .eq. 2) then
               do 140 n = 0, jr-1
                  fac1 = (jr-n) * fac2
                  if (n .eq. 0) fac1 = 0.5d0 * fac1
                  do 140 m = idir, idir*(ir-1), idir
                     fac = (ir-abs(m)) * fac1
                     do 140 j = regl1, regh1
 140                    dest(i,j,k) = dest(i,j,k) +
     @                    fac * (src(i*ir+m,j*jr-n,k*kr) +
     @                           src(i*ir+m,j*jr+n,k*kr))
            else if (integ .eq. 0) then
               center = center + cfac
            endif
 150        continue
         fac2 = ir * fac0
         cfac = 0.5d0 * cube * fac0 * jr * (kr-1)
         do 170 ki = 0, 1
            kdir = 2 * ki - 1
            if (ga(0,0,ki) + ga(1,0,ki) .eq. 2) then
               do 160 l = kdir, kdir*(kr-1), kdir
                  fac1 = (kr-abs(l)) * fac2
                  do 160 n = 0, jr-1
                     fac = (jr-n) * fac1
                     if (n .eq. 0) fac = 0.5d0 * fac
                     do 160 j = regl1, regh1
 160                    dest(i,j,k) = dest(i,j,k) +
     @                    fac * (src(i*ir,j*jr-n,k*kr+l) +
     @                           src(i*ir,j*jr+n,k*kr+l))
            else if (integ .eq. 0) then
               center = center + cfac
            endif
 170        continue
c center
         do 180 j = regl1, regh1
 180        dest(i,j,k) = dest(i,j,k) +
     @        center * src(i*ir,j*jr,k*kr)
      else
         i = regl0
         j = regl1
         do 210 k = regl2, regh2
 210        dest(i,j,k) = 0.d0
c center gets center plus two edges
         if (integ .eq. 0) then
            center = kr / cube
            fac0 = 1.d0 / (cube**2)
            cfac = 0.25d0 * cube * fac0 * kr * (ir-1) * (jr-1)
         else
            center = 1.d0
            fac0 = 1.d0 / cube
         endif
c quadrants
c each quadrant is two octants and a face
         do 230 ji = 0, 1
            jdir = 2 * ji - 1
            do 230 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,0) .eq. 1) then
                  do 220 l = 0, kr-1
                     fac2 = (kr-l) * fac0
                     if (l .eq. 0) fac2 = 0.5d0 * fac2
                     do 220 n = jdir, jdir*(jr-1), jdir
                        fac1 = (jr-abs(n)) * fac2
                        do 220 m = idir, idir*(ir-1), idir
                           fac = (ir-abs(m)) * fac1
                           do 220 k = regl2, regh2
 220                          dest(i,j,k) = dest(i,j,k) +
     @                          fac * (src(i*ir+m,j*jr+n,k*kr-l) +
     @                                 src(i*ir+m,j*jr+n,k*kr+l))
               else if (integ .eq. 0) then
                  center = center + cfac
               endif
 230           continue
c faces
c each face is two faces and an edge
         fac2 = jr * fac0
         cfac = 0.5d0 * cube * fac0 * kr * (ir-1)
         do 250 ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,0,0) + ga(ii,1,0) .eq. 2) then
               do 240 l = 0, kr-1
                  fac1 = (kr-l) * fac2
                  if (l .eq. 0) fac1 = 0.5d0 * fac1
                  do 240 m = idir, idir*(ir-1), idir
                     fac = (ir-abs(m)) * fac1
                     do 240 k = regl2, regh2
 240                    dest(i,j,k) = dest(i,j,k) +
     @                    fac * (src(i*ir+m,j*jr,k*kr-l) +
     @                           src(i*ir+m,j*jr,k*kr+l))
            else if (integ .eq. 0) then
               center = center + cfac
            endif
 250        continue
         fac2 = ir * fac0
         cfac = 0.5d0 * cube * fac0 * kr * (jr-1)
         do 270 ji = 0, 1
            jdir = 2 * ji - 1
            if (ga(0,ji,0) + ga(1,ji,0) .eq. 2) then
               do 260 l = 0, kr-1
                  fac1 = (kr-l) * fac2
                  if (l .eq. 0) fac1 = 0.5d0 * fac1
                  do 260 n = jdir, jdir*(jr-1), jdir
                     fac = (jr-abs(n)) * fac1
                     do 260 k = regl2, regh2
 260                    dest(i,j,k) = dest(i,j,k) +
     @                    fac * (src(i*ir,j*jr+n,k*kr-l) +
     @                           src(i*ir,j*jr+n,k*kr+l))
            else if (integ .eq. 0) then
               center = center + cfac
            endif
 270        continue
c center
         do 280 k = regl2, regh2
 280        dest(i,j,k) = dest(i,j,k) +
     @        center * src(i*ir,j*jr,k*kr)
      endif
      end

c NODE-based data only.
c Fills coarse region defined by reg.
c Handles any corner geometry except all-coarse.
      subroutine ancr2(dest,
     @ dims(dest),
     @ dims(reg),
     @ src,  dims(src),
     @ ir, jr, kr, ga, integ)
      integer dims(dest)
      integer dims(reg)
      integer dims(src)
      integer ir, jr, kr, ga(0:1,0:1,0:1), integ
      real*8 dest(dimdec(dest))
      real*8 src(dimdec(src))
      real*8 cube, center, cfac, fac0, fac1, fac2, fac
      integer i, j, k, ii, ji, ki, idir, jdir, kdir, l, m, n
      i = regl0
      j = regl1
      k = regl2
      dest(i,j,k) = 0.d0
      cube = ir * jr * kr
      if (integ .eq. 0) then
         center = 1.d0 / cube
         fac0 = 1.d0 / (cube**2)
         cfac = 0.125d0 * cube * fac0 * (ir-1) * (jr-1) * (kr-1)
      else
         center = 1.d0
         fac0 = 1.d0 / cube
      endif
c octants
      do 20 ki = 0, 1
         kdir = 2 * ki - 1
         do 20 ji = 0, 1
            jdir = 2 * ji - 1
            do 20 ii = 0, 1
               idir = 2 * ii - 1
               if (ga(ii,ji,ki) .eq. 1) then
                  do 10 l = kdir, kdir*(kr-1), kdir
                     fac2 = (kr-abs(l)) * fac0
                     do 10 n = jdir, jdir*(jr-1), jdir
                        fac1 = (jr-abs(n)) * fac2
                        do 10 m = idir, idir*(ir-1), idir
                           fac = (ir-abs(m)) * fac1
 10                        dest(i,j,k) = dest(i,j,k) +
     @                       fac * src(i*ir+m,j*jr+n,k*kr+l)
               else if (integ .eq. 0) then
                  center = center + cfac
               endif
 20            continue
c faces
      fac2 = kr * fac0
      cfac = 0.25d0 * cube * fac0 * (ir-1) * (jr-1)
      do 40 ji = 0, 1
         jdir = 2 * ji - 1
         do 40 ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,ji,0) + ga(ii,ji,1) .eq. 2) then
               do 30 n = jdir, jdir*(jr-1), jdir
                  fac1 = (jr-abs(n)) * fac2
                  do 30 m = idir, idir*(ir-1), idir
                     fac = (ir-abs(m)) * fac1
 30                  dest(i,j,k) = dest(i,j,k) +
     @                 fac * src(i*ir+m,j*jr+n,k*kr)
            else if (integ .eq. 0) then
               center = center + cfac
            endif
 40         continue
      fac2 = jr * fac0
      cfac = 0.25d0 * cube * fac0 * (ir-1) * (kr-1)
      do 60 ki = 0, 1
         kdir = 2 * ki - 1
         do 60 ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,0,ki) + ga(ii,1,ki) .eq. 2) then
               do 50 l = kdir, kdir*(kr-1), kdir
                  fac1 = (kr-abs(l)) * fac2
                  do 50 m = idir, idir*(ir-1), idir
                     fac = (ir-abs(m)) * fac1
 50                  dest(i,j,k) = dest(i,j,k) +
     @                 fac * src(i*ir+m,j*jr,k*kr+l)
            else if (integ .eq. 0) then
               center = center + cfac
            endif
 60         continue
      fac2 = ir * fac0
      cfac = 0.25d0 * cube * fac0 * (jr-1) * (kr-1)
      do 80 ki = 0, 1
         kdir = 2 * ki - 1
         do 80 ji = 0, 1
            jdir = 2 * ji - 1
            if (ga(0,ji,ki) + ga(1,ji,ki) .eq. 2) then
               do 70 l = kdir, kdir*(kr-1), kdir
                  fac1 = (kr-abs(l)) * fac2
                  do 70 n = jdir, jdir*(jr-1), jdir
                     fac = (jr-abs(n)) * fac1
 70                  dest(i,j,k) = dest(i,j,k) +
     @                 fac * src(i*ir,j*jr+n,k*kr+l)
            else if (integ .eq. 0) then
               center = center + cfac
            endif
 80         continue
c edges
      fac1 = jr * kr * fac0
      cfac = 0.5d0 * cube * fac0 * (ir-1)
      do 100 ii = 0, 1
         idir = 2 * ii - 1
         if (ga(ii,0,0) + ga(ii,0,1) +
     @       ga(ii,1,0) + ga(ii,1,1) .eq. 4) then
            do 90 m = idir, idir*(ir-1), idir
               fac = (ir-abs(m)) * fac1
 90            dest(i,j,k) = dest(i,j,k) +
     @           fac * src(i*ir+m,j*jr,k*kr)
         else if (integ .eq. 0) then
            center = center + cfac
         endif
 100     continue
      fac1 = ir * kr * fac0
      cfac = 0.5d0 * cube * fac0 * (jr-1)
      do 120 ji = 0, 1
         jdir = 2 * ji - 1
         if (ga(0,ji,0) + ga(0,ji,1) +
     @       ga(1,ji,0) + ga(1,ji,1) .eq. 4) then
            do 110 n = jdir, jdir*(jr-1), jdir
               fac = (jr-abs(n)) * fac1
 110           dest(i,j,k) = dest(i,j,k) +
     @           fac * src(i*ir,j*jr+n,k*kr)
         else if (integ .eq. 0) then
            center = center + cfac
         endif
 120     continue
      fac1 = ir * jr * fac0
      cfac = 0.5d0 * cube * fac0 * (kr-1)
      do 140 ki = 0, 1
         kdir = 2 * ki - 1
         if (ga(0,0,ki) + ga(0,1,ki) +
     @       ga(1,0,ki) + ga(1,1,ki) .eq. 4) then
            do 130 l = kdir, kdir*(kr-1), kdir
               fac = (kr-abs(l)) * fac1
 130           dest(i,j,k) = dest(i,j,k) +
     @           fac * src(i*ir,j*jr,k*kr+l)
         else if (integ .eq. 0) then
            center = center + cfac
         endif
 140     continue
c center
      dest(i,j,k) = dest(i,j,k) +
     @  center * src(i*ir,j*jr,k*kr)
      end
