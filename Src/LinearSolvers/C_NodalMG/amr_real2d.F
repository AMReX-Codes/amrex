c-----------------------------------------------------------------------
c Works for NODE-based data.
      subroutine fcpy2(
     & da,    dal0, dah0, dal1, dah1,
     & db,    dbl0, dbh0, dbl1, dbh1,
     &        regl0, regh0, regl1, regh1, w, ncomp)
      integer dal0, dah0, dal1, dah1
      integer dbl0, dbh0, dbl1, dbh1
      integer regl0, regh0, regl1, regh1
      integer w, ncomp
      double precision da(dal0:dah0,dal1:dah1, ncomp)
      double precision db(dbl0:dbh0,dbl1:dbh1, ncomp)
      integer i, j, k
      do k = 1, ncomp
         do j = 1, w
            do i = regl0 - w, regh0 + w
               da(i,regl1+j,k) = db(i,regl1+j,k)
               db(i,regl1-j,k) = da(i,regl1-j,k)
            end do
         end do
      end do
      end
c-----------------------------------------------------------------------
c Works for CELL-based data.
      subroutine iprodc(
     & v0,    v0l0, v0h0, v0l1, v0h1,
     & v1,    v1l0, v1h0, v1l1, v1h1,
     &        regl0, regh0, regl1, regh1, 
     & sum)
      integer v0l0, v0h0, v0l1, v0h1
      integer v1l0, v1h0, v1l1, v1h1
      integer regl0, regh0, regl1, regh1
      double precision v0(v0l0:v0h0,v0l1:v0h1)
      double precision v1(v1l0:v1h0,v1l1:v1h1)
      double precision sum
      integer i, j
      do j = regl1, regh1
         do i = regl0, regh0
            sum = sum + v0(i,j) * v1(i,j)
         end do
      end do
      end
c-----------------------------------------------------------------------
c Works for NODE-based data.
      subroutine iprodn(
     & v0,    v0l0, v0h0, v0l1, v0h1,
     & v1,    v1l0, v1h0, v1l1, v1h1,
     &        regl0, regh0, regl1, regh1, 
     & sum)
      integer v0l0, v0h0, v0l1, v0h1
      integer v1l0, v1h0, v1l1, v1h1
      integer regl0, regh0, regl1, regh1
      double precision v0(v0l0:v0h0,v0l1:v0h1)
      double precision v1(v1l0:v1h0,v1l1:v1h1)
      double precision sum, sum0
      integer i, j
      sum0 = 0.5d0 * (v0(regl0,regl1) * v1(regl0,regl1) +
     &                v0(regl0,regh1) * v1(regl0,regh1) +
     &                v0(regh0,regl1) * v1(regh0,regl1) +
     &                v0(regh0,regh1) * v1(regh0,regh1))
      do i = regl0 + 1, regh0 - 1
         sum0 = sum0 + (v0(i,regl1) * v1(i,regl1) +
     &                  v0(i,regh1) * v1(i,regh1))
      end do
      do j = regl1 + 1, regh1 - 1
         sum0 = sum0 + (v0(regl0,j) * v1(regl0,j) +
     &                  v0(regh0,j) * v1(regh0,j))
      end do
      sum = sum + 0.5d0 * sum0
      do j = regl1 + 1, regh1 - 1
         do i = regl0 + 1, regh0 - 1
            sum = sum + v0(i,j) * v1(i,j)
         end do
      end do
      end
c-----------------------------------------------------------------------
c Works for CELL- or NODE-based data.
      subroutine bref(
     & dest,  destl0, desth0, destl1, desth1,
     &        regl0, regh0, regl1, regh1,
     & src,   srcl0, srch0, srcl1, srch1,
     &        bbl0, bbh0, bbl1, bbh1, 
     & idir, ncomp)
      integer ncomp
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer srcl0, srch0, srcl1, srch1
      integer bbl0, bbh0, bbl1, bbh1
      integer idir
      double precision dest(destl0:desth0,destl1:desth1,ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,ncomp)
      integer i, j, nc
      do nc = 1, ncomp
      if (idir .eq. 0) then
         do i = regl0, regh0
            do j = regl1, regh1
               dest(i,j,nc) = src(bbh0-(i-regl0),j,nc)
            end do
         end do
      else
         do j = regl1, regh1
            do i = regl0, regh0
               dest(i,j,nc) = src(i,bbh1-(j-regl1),nc)
            end do
         end do
      end if
      end do
      end
c-----------------------------------------------------------------------
c Works for CELL- or NODE-based data.
      subroutine brefm(
     & dest,  destl0, desth0, destl1, desth1,
     &        regl0, regh0, regl1, regh1,
     & src,   srcl0, srch0, srcl1, srch1,
     & bbl0,  bbh0, bbl1, bbh1, ra, ncomp)
      integer ncomp
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer srcl0, srch0, srcl1, srch1
      integer bbl0, bbh0, bbl1, bbh1
      integer ra(0:1)
      double precision dest(destl0:desth0,destl1:desth1,ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,ncomp)
      integer i, j, nc
      do nc = 1, ncomp
      if (ra(0) .eq. 0 .and. ra(1) .eq. 0) then
         do j = regl1, regh1
            do i = regl0, regh0
                  dest(i,j,nc) = src(bbl0+(i-regl0),bbl1+(j-regl1),nc)
            end do
         end do
      else if (ra(0) .eq. 0 .and. ra(1) .eq. 1) then
         do j = regl1, regh1
            do i = regl0, regh0
               dest(i,j,nc) = src(bbl0+(i-regl0),bbh1-(j-regl1),nc)
            end do
         end do
      else if (ra(0) .eq. 1 .and. ra(1) .eq. 0) then
         do i = regl0, regh0
            do j = regl1, regh1
               dest(i,j,nc) = src(bbh0-(i-regl0),bbl1+(j-regl1),nc)
            end do
         end do
      else
         do j = regl1, regh1
            do i = regl0, regh0
               dest(i,j,nc) = src(bbh0-(i-regl0),bbh1-(j-regl1),nc)
            end do
         end do
      end if
      end do
      end
c-----------------------------------------------------------------------
c Works for CELL- or NODE-based data.
      subroutine bneg(
     & dest,  destl0, desth0, destl1, desth1,
     &        regl0, regh0, regl1, regh1,
     & src,   srcl0, srch0, srcl1, srch1,
     &        bbl0, bbh0, bbl1, bbh1, idir, ncomp)
      integer ncomp
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer srcl0, srch0, srcl1, srch1
      integer bbl0, bbh0, bbl1, bbh1
      integer idir
      double precision dest(destl0:desth0,destl1:desth1,ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,ncomp)
      integer i, j, nc
      do nc = 1, ncomp
      if (idir .eq. 0) then
         do i = regl0, regh0
            do j = regl1, regh1
               dest(i,j,nc) = -src(bbh0-(i-regl0),j,nc)
            end do
         end do
      else
         do j = regl1, regh1
            do i = regl0, regh0
               dest(i,j,nc) = -src(i,bbh1-(j-regl1),nc)
            end do
         end do
      end if
      end do
      end
c-----------------------------------------------------------------------
c Works for CELL- or NODE-based data.
      subroutine bnegm(
     & dest,  destl0, desth0, destl1, desth1,
     &        regl0, regh0, regl1, regh1,
     & src,   srcl0, srch0, srcl1, srch1,
     &        bbl0, bbh0, bbl1, bbh1, ra, ncomp)
      integer ncomp
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer srcl0, srch0, srcl1, srch1
      integer bbl0, bbh0, bbl1, bbh1
      integer ra(0:1)
      double precision dest(destl0:desth0,destl1:desth1, ncomp)
      double precision src(srcl0:srch0,srcl1:srch1, ncomp)
      integer i, j, nc
      do nc = 1, ncomp
      if (ra(0) .eq. 0 .and. ra(1) .eq. 0) then
         do j = regl1, regh1
            do i = regl0, regh0
               dest(i,j,nc) = -src(bbl0+(i-regl0),bbl1+(j-regl1),nc)
            end do
         end do
      else if (ra(0) .eq. 0 .and. ra(1) .eq. 1) then
         do j = regl1, regh1
            do i = regl0, regh0
               dest(i,j,nc) = -src(bbl0+(i-regl0),bbh1-(j-regl1),nc)
            end do
         end do
      else if (ra(0) .eq. 1 .and. ra(1) .eq. 0) then
         do i = regl0, regh0
            do j = regl1, regh1
               dest(i,j,nc) = -src(bbh0-(i-regl0),bbl1+(j-regl1),nc)
            end do
         end do
      else
         do i = regl0, regh0
            do j = regl1, regh1
               dest(i,j,nc) = -src(bbh0-(i-regl0),bbh1-(j-regl1),nc)
            end do
         end do
      end if
      end do
      end
c-----------------------------------------------------------------------
c Works for CELL-based velocity data.
c This routine assumes that the inflow face velocity data has not yet
c been altered.  Running fill_borders should call this routine on every
c inflow face, so that binfil can be run for subsequent fills
      subroutine binflo(
     & dest,  destl0, desth0, destl1, desth1,
     &        regl0, regh0, regl1, regh1,
     & src,   srcl0, srch0, srcl1, srch1,
     &        bbl0, bbh0, bbl1, bbh1,
     & idir, ncomp)
      integer ncomp
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer srcl0, srch0, srcl1, srch1
      integer bbl0, bbh0, bbl1, bbh1
      integer idir
      double precision dest(destl0:desth0,destl1:desth1,ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,ncomp)
      integer i, j,nc
      do nc = 1, ncomp
      if (idir .eq. 0) then
         if (regl0 .lt. bbh0) then
            do i = regl0, regh0
               do j = regl1, regh1
                  dest(i,j,nc) = 2.d0 * dest(regh0,j,nc) 
     &			- src(bbh0-(i-regl0),j,nc)
               end do
            end do
         else
            do i = regh0, regl0, -1
               do j = regl1, regh1
                  dest(i,j,nc) = 2.d0 * dest(regl0,j,nc) 
     &			- src(bbh0-(i-regl0),j,nc)
               end do
            end do
         end if
      else
         if (regl1 .lt. bbh1) then
            do j = regl1, regh1
               do i = regl0, regh0
                  dest(i,j,nc) = 2.d0 * dest(i,regh1,nc) 
     &			- src(i,bbh1-(j-regl1),nc)
               end do
            end do
         else
            do j = regh1, regl1, -1
               do i = regl0, regh0
                  dest(i,j,nc) = 2.d0 * dest(i,regl1,nc) 
     &			- src(i,bbh1-(j-regl1),nc)
               end do
            end do
         end if
      end if
      end do
      end
c-----------------------------------------------------------------------
c Works for CELL-based velocity data.
c This routine is called when the inflow face velocity data has already
c been altered by a call to fill_borders.  The box bb must have been
c extended to one cell past the boundary by the boundary::box routine.

      subroutine binfil(
     & dest,  destl0, desth0, destl1, desth1,
     &        regl0, regh0, regl1, regh1,
     & src,   srcl0, srch0, srcl1, srch1,
     &        bbl0, bbh0, bbl1, bbh1,
     & idir, ncomp)
      integer ncomp
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer srcl0, srch0, srcl1, srch1
      integer bbl0, bbh0, bbl1, bbh1
      integer idir
      double precision dest(destl0:desth0,destl1:desth1,ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,ncomp)
      integer i, j, nc
      do nc = 1, ncomp
      if (idir .eq. 0) then
         if (regl0 .lt. bbh0) then
            do i = regl0, regh0
               do j = regl1, regh1
                  dest(i,j,nc) = src(bbl0,j,nc) + src(bbl0+1,j,nc) -
     &                        src(bbh0-(i-regl0),j,nc)
               end do
            end do
         else
            do i = regh0, regl0, -1
               do j = regl1, regh1
                  dest(i,j,nc) = src(bbh0-1,j,nc) + src(bbh0,j,nc) -
     &                        src(bbl0+(regh0-i),j,nc)
               end do
            end do
         end if
      else
         if (regl1 .lt. bbh1) then
            do j = regl1, regh1
               do i = regl0, regh0
                  dest(i,j,nc) = src(i,bbl1,nc) + src(i,bbl1+1,nc) -
     &                        src(i,bbh1-(j-regl1),nc)
               end do
            end do
         else
            do j = regh1, regl1, -1
               do i = regl0, regh0
                  dest(i,j,nc) = src(i,bbh1-1,nc) + src(i,bbh1,nc) -
     &                        src(i,bbl1+(regh1-j),nc)
               end do
            end do
         end if
      end if
      end do
      end
c-----------------------------------------------------------------------
c     Interpolation...
c-----------------------------------------------------------------------
c CELL-based data only.
      subroutine acint2(
     & dest,  destl0, desth0, destl1, desth1,
     &        regl0, regh0, regl1, regh1,
     & src,   srcl0, srch0, srcl1, srch1,
     &        bbl0, bbh0, bbl1, bbh1, 
     & ir, jr, ncomp, integ, i1, i2)
      integer ncomp
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer srcl0, srch0, srcl1, srch1
      integer bbl0, bbh0, bbl1, bbh1
      integer ir, jr
	integer i1, i2, integ
      double precision dest(destl0:desth0,destl1:desth1,ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,ncomp)
      double precision xoff, yoff, sx, sy
      integer i, j, ic, jc, nc
      do nc = 1, ncomp
      do j = regl1, regh1
         jc = j/jr
         yoff = (mod(j,jr) + 0.5) / jr - 0.5
         do i = regl0, regh0
            ic = i/ir
            xoff = (mod(i,ir) + 0.5) / ir - 0.5
            sy = 0.5 * (src(ic,jc+1,nc) - src(ic,jc-1,nc))
            sx = 0.5 * (src(ic+1,jc,nc) - src(ic-1,jc,nc))
            dest(i,j,nc) = src(ic,jc,nc) + xoff * sx + yoff * sy
         end do
      end do
      end do
      end
c-----------------------------------------------------------------------
c NODE-based data only.
      subroutine anint2(
     & dest,  destl0, desth0, destl1, desth1,
     &        regl0, regh0, regl1, regh1,
     & src,   srcl0, srch0, srcl1, srch1,
     &        bbl0, bbh0, bbl1, bbh1, 
     & ir, jr, ncomp, i1, i2, integ)
      integer ncomp
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer srcl0, srch0, srcl1, srch1
      integer bbl0, bbh0, bbl1, bbh1
      integer ir, jr,nc
	integer i1, i2, integ
      double precision dest(destl0:desth0,destl1:desth1, ncomp)
      double precision src(srcl0:srch0,srcl1:srch1, ncomp)
      double precision p, q
      integer ic, jc, j, m
      do nc = 1, ncomp
         do jc = bbl1, bbh1
            do ic = bbl0, bbh0
               dest(ir*ic,jr*jc,nc) = src(ic,jc,nc)
            end do
         end do
         do m = 1, jr-1
            q = dble(m)/jr
            p = 1.0d0 - q
            do jc = bbl1, bbh1-1
               do ic = bbl0, bbh0
                  dest(ir*ic,jr*jc+m,nc) = p * src(ic,jc,nc) +
     &                                  q * src(ic,jc+1,nc)
               end do
            end do
         end do
         do ic = bbl0, bbh0-1
cdir$ ivdep
            do j = regl1, regh1
               dest(2*ic+1,j,nc) = 0.5d0 * (dest(2*ic,j,nc) +
     &                                   dest(2*ic+2,j,nc))
            end do
         end do
         do m = 1, ir-1
            q = dble(m)/ir
            p = 1.0d0 - q
            do ic = bbl0, bbh0-1
cdir$ ivdep
               do j = regl1, regh1
                  dest(ir*ic+m,j,nc) = p * dest(ir*ic,j,nc) +
     &                              q * dest(ir*(ic+1),j,nc)
               end do
            end do
         end do
      end do
      end
c-----------------------------------------------------------------------
c     Restrictions....
c-----------------------------------------------------------------------
c CELL-based data only.
      subroutine acrst1(
     & dest,  destl0, desth0, destl1, desth1,
     &        regl0, regh0, regl1, regh1,
     & src,   srcl0, srch0, srcl1, srch1,
     & ir, jr, ncomp, integ, i1, i2)
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer srcl0, srch0, srcl1, srch1
      integer ir, jr, ncomp, integ
	integer i1, i2
      double precision dest(destl0:desth0,destl1:desth1, ncomp)
      double precision src(srcl0:srch0,srcl1:srch1, ncomp)
      double precision fac
      integer i, j, m, n, nc
	do nc = 1, ncomp
         do j = regl1, regh1
            do i = regl0, regh0
               dest(i,j,nc) = 0.d0
            end do
         end do
         do n = 0, jr-1
            do m = 0, ir-1
               do j = regl1, regh1
                  do i = regl0, regh0
                     dest(i,j,nc) = dest(i,j,nc) + src(i*ir+m,j*jr+n,nc)
                  end do
               end do
            end do
         end do
         if (integ .eq. 0) then
            fac = 1.d0 / (ir*jr)
            do j = regl1, regh1
               do i = regl0, regh0
                  dest(i,j,nc) = dest(i,j,nc) * fac
               end do
            end do
         end if
	end do
      end
c-----------------------------------------------------------------------
c NODE-based data only.
      subroutine anrst1(
     & dest,  destl0, desth0, destl1, desth1,
     &        regl0, regh0, regl1, regh1,
     & src,   srcl0, srch0, srcl1, srch1,
     & ir, jr, ncomp, integ, i1, i2)
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer srcl0, srch0, srcl1, srch1
      integer ir, jr, ncomp
	integer i1, i2, integ
      double precision dest(destl0:desth0,destl1:desth1, ncomp)
      double precision src(srcl0:srch0,srcl1:srch1, ncomp)
      integer i, j,nc
	do nc = 1, ncomp
         do j = regl1, regh1
            do i = regl0, regh0
               dest(i,j,nc) = src(i*ir,j*jr,nc)
            end do
         end do
	end do
      end
c-----------------------------------------------------------------------
c NODE-based data only.
c Fills coarse region defined by reg, which must be smaller than
c fine region by at least one coarse cell on all sides.
      subroutine anrst2(
     & dest,  destl0, desth0, destl1, desth1,
     &        regl0, regh0, regl1, regh1,
     & src,   srcl0, srch0, srcl1, srch1,
     & ir, jr, ncomp, integ, i1, i2)
      integer ncomp
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer srcl0, srch0, srcl1, srch1
      integer ir, jr, integ, i1, i2
      double precision dest(destl0:desth0,destl1:desth1,ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,ncomp)
      double precision fac0, fac1, fac
      integer i, j, m, n, nc
      do nc = 1, ncomp
         do j = regl1, regh1
            do i = regl0, regh0
               dest(i,j,nc) = 0.d0
            end do
         end do
         fac0 = 1.d0 / (ir*ir * jr*jr)
         do n = 0, jr-1
            fac1 = (jr-n) * fac0
            if (n .eq. 0) fac1 = 0.5d0 * fac1
            do m = 0, ir-1
               fac = (ir-m) * fac1
               if (m .eq. 0) fac = 0.5d0 * fac
               do j = regl1, regh1
                  do i = regl0, regh0
                     dest(i,j,nc) = dest(i,j,nc) +
     &                           fac * (src(i*ir-m,j*jr-n,nc)+
     &                                  src(i*ir-m,j*jr+n,nc)+
     &                                  src(i*ir+m,j*jr-n,nc)+
     &                                  src(i*ir+m,j*jr+n,nc))
                  end do
               end do
            end do
         end do
         if (integ .eq. 1) then
            fac = ir * jr
            do j = regl1, regh1
               do i = regl0, regh0
                  dest(i,j,nc) = fac * dest(i,j,nc)
               end do
            end do
         end if
      end do
      end
c-----------------------------------------------------------------------
c NODE-based data only.
c Fills coarse region defined by reg.
      subroutine anfr2(
     & dest,  destl0, desth0, destl1, desth1,
     &        regl0, regh0, regl1, regh1,
     & src,   srcl0, srch0, srcl1, srch1,
     & ir, jr, ncomp, integ, idim, idir)
      integer ncomp
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer srcl0, srch0, srcl1, srch1
      integer ir, jr, idim, idir, integ
      double precision dest(destl0:desth0,destl1:desth1,ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,ncomp)
      double precision fac0, fac1, fac
      integer i, j, m, n, nc
      do nc = 1, ncomp
         if (idim .eq. 0) then
            if (integ .eq. 0) then
               fac = (0.5d0 + 0.5d0 / ir)
               fac0 = 1.d0 / (ir*ir * jr*jr)
            else
               fac = 1.d0
               fac0 = 1.d0 / (ir * jr)
            end if
            i = regl0
            do j = regl1, regh1
               dest(i,j,nc) = fac * src(i*ir,j*jr,nc)
            end do
            do n = 0, jr-1
               fac1 = (jr-n) * fac0
               if (n .eq. 0) fac1 = 0.5d0 * fac1
               do m = idir, idir*(ir-1), idir
                  fac = (ir-abs(m)) * fac1
                  do j = regl1, regh1
                     dest(i,j,nc) = dest(i,j,nc) +
     &			 fac * (src(i*ir+m,j*jr-n,nc)+
     &                              src(i*ir+m,j*jr+n,nc))
                  end do
               end do
            end do
         else
            if (integ .eq. 0) then
               fac = (0.5d0 + 0.5d0 / jr)
               fac0 = 1.d0 / (ir*ir * jr*jr)
            else
               fac = 1.d0
               fac0 = 1.d0 / (ir * jr)
            end if
            j = regl1
            do i = regl0, regh0
               dest(i,j,nc) = fac * src(i*ir,j*jr,nc)
            end do
            do n = idir, idir*(jr-1), idir
               fac1 = (jr-abs(n)) * fac0
               do m = 0, ir-1
                  fac = (ir-m) * fac1
                  if (m .eq. 0) fac = 0.5d0 * fac
                  do  i = regl0, regh0
                     dest(i,j,nc) = dest(i,j,nc)
     &		+ fac * (src(i*ir-m,j*jr+n,nc)
     &                        +src(i*ir+m,j*jr+n,nc))
                  end do
               end do
            end do
         end if
      end do
      end
c-----------------------------------------------------------------------
c NODE-based data only.
c Fills coarse point defined by reg.
      subroutine anor2(
     & dest,  destl0, desth0, destl1, desth1,
     &        regl0, regh0, regl1, regh1,
     & src,   srcl0, srch0, srcl1, srch1,
     & ir, jr, ncomp, integ, idir0, idir1)
      integer ncomp
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer srcl0, srch0, srcl1, srch1
      integer ir, jr, idir0, idir1, integ
      double precision dest(destl0:desth0,destl1:desth1,ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,ncomp)
      double precision fac0, fac1, fac
      integer i, j, m, n, nc
      do nc = 1, ncomp
         if (integ .eq. 0) then
            fac = (0.75d0 + 0.25d0/ir + 0.25d0/jr - 0.25d0/(ir*jr))
            fac0 = 1.d0 / (ir*ir * jr*jr)
         else
            fac = 1.d0
            fac0 = 1.d0 / (ir * jr)
         end if
         i = regl0
         j = regl1
         dest(i,j,nc) = fac * src(i*ir,j*jr,nc)
         do n = idir1, idir1*(jr-1), idir1
            fac1 = (jr-abs(n)) * fac0
            do m = idir0, idir0*(ir-1), idir0
               fac = (ir-abs(m)) * fac1
               dest(i,j,nc) = dest(i,j,nc) + fac * src(i*ir+m,j*jr+n,nc)
            end do
         end do
      end do
      end
c-----------------------------------------------------------------------
c NODE-based data only.
c Fills coarse point defined by reg.
      subroutine anir2(
     & dest,  destl0, desth0, destl1, desth1,
     &        regl0, regh0, regl1, regh1,
     & src,   srcl0, srch0, srcl1, srch1,
     & ir, jr, ncomp, integ, idir0, idir1)
      integer ncomp
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer srcl0, srch0, srcl1, srch1
      integer ir, jr, idir0, idir1, integ
      double precision dest(destl0:desth0,destl1:desth1,ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,ncomp)
      double precision fac0, fac1, fac
      integer i, j, m, n, nc
      do nc = 1, ncomp
         if (integ .eq. 0) then
            fac = (0.25d0 + 0.25d0/ir + 0.25d0/jr + 0.25d0/(ir*jr))
            fac0 = 1.d0/(ir*ir * jr*jr)
         else
            fac = 1.d0
            fac0 = 1.d0/(ir*jr)
         end if
         i = regl0
         j = regl1
         dest(i,j,nc) = fac * src(i*ir,j*jr,nc)
         do m = idir0, idir0*(ir-1), idir0
            fac1 = (ir-abs(m)) * fac0
            fac = jr * fac1
            dest(i,j,nc) = dest(i,j,nc) + fac * src(i*ir+m,j*jr,nc)
         end do
         do n = idir1, idir1*(jr-1), idir1
            fac1 = (jr-abs(n)) * fac0
            fac = ir * fac1
            dest(i,j,nc) = dest(i,j,nc) + fac * src(i*ir,j*jr+n,nc)
            do m = idir0, idir0*(ir-1), idir0
               fac = (ir-abs(m)) * fac1
               dest(i,j,nc) = dest(i,j,nc) + fac * 
     &               (src(i*ir+m,j*jr+n,nc) +
     &                src(i*ir-m,j*jr+n,nc) +
     &                src(i*ir+m,j*jr-n,nc))
            end do
         end do
      end do
      end
c-----------------------------------------------------------------------
c NODE-based data only.
c Fills coarse point defined by reg.
      subroutine andr2(
     & dest,  destl0, desth0, destl1, desth1,
     &        regl0, regh0, regl1, regh1,
     & src,   srcl0, srch0, srcl1, srch1,
     & ir, jr, ncomp, integ, idir1, i2)
      integer ncomp
      integer destl0, desth0, destl1, desth1
      integer regl0, regh0, regl1, regh1
      integer srcl0, srch0, srcl1, srch1
      integer ir, jr, idir1, integ
      double precision dest(destl0:desth0,destl1:desth1,ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,ncomp)
      double precision fac0, fac1, fac
      integer i, j, m, n, nc, i2
      do nc = 1, ncomp
         if (integ .eq. 0) then
            fac = (0.5d0 + 0.5d0 / ir + 0.5d0 / jr - 0.5d0 / (ir*jr))
            fac0 = 1.d0 / (ir*ir * jr*jr)
         else
            fac = 1.d0
            fac0 = 1.d0 / (ir * jr)
         end if
         i = regl0
         j = regl1
         dest(i,j,nc) = fac * src(i*ir,j*jr,nc)
         do n = idir1, idir1*(jr-1), idir1
            fac1 = (jr-abs(n)) * fac0
            do m = 1, ir-1
               fac = (ir-m) * fac1
               dest(i,j,nc) = dest(i,j,nc) + fac * 
     &               (src(i*ir+m,j*jr+n,nc) +
     &                src(i*ir-m,j*jr-n,nc))
            end do
         end do
      end do
      end
