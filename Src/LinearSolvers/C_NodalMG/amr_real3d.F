c-----------------------------------------------------------------------
c Works for CELL-based data.
      subroutine iprodc(
     & v0, v0l0,v0h0,v0l1,v0h1,v0l2,v0h2,
     & v1, v1l0,v1h0,v1l1,v1h1,v1l2,v1h2,
     &     regl0,regh0,regl1,regh1,regl2,regh2, 
     & sum)
      integer v0l0,v0h0,v0l1,v0h1,v0l2,v0h2
      integer v1l0,v1h0,v1l1,v1h1,v1l2,v1h2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      double precision v0(v0l0:v0h0,v0l1:v0h1,v0l2:v0h2)
      double precision v1(v1l0:v1h0,v1l1:v1h1,v1l2:v1h2)
      double precision sum
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               sum = sum + v0(i,j,k) * v1(i,j,k)
            end do
         end do
      end do
      end
c-----------------------------------------------------------------------
c Works for CELL- or NODE-based data.
      subroutine iprodn(
     & v0, v0l0,v0h0,v0l1,v0h1,v0l2,v0h2,
     & v1, v1l0,v1h0,v1l1,v1h1,v1l2,v1h2,
     &     regl0,regh0,regl1,regh1,regl2,regh2, 
     & sum)
      integer v0l0,v0h0,v0l1,v0h1,v0l2,v0h2
      integer v1l0,v1h0,v1l1,v1h1,v1l2,v1h2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      double precision v0(v0l0:v0h0,v0l1:v0h1,v0l2:v0h2)
      double precision v1(v1l0:v1h0,v1l1:v1h1,v1l2:v1h2)
      double precision sum
      integer i, j, k
      sum = sum + 0.125d0 *
     &            (v0(regl0,regl1,regl2) * v1(regl0,regl1,regl2) +
     &             v0(regl0,regl1,regh2) * v1(regl0,regl1,regh2) +
     &             v0(regl0,regh1,regl2) * v1(regl0,regh1,regl2) +
     &             v0(regl0,regh1,regh2) * v1(regl0,regh1,regh2) +
     &             v0(regh0,regl1,regl2) * v1(regh0,regl1,regl2) +
     &             v0(regh0,regl1,regh2) * v1(regh0,regl1,regh2) +
     &             v0(regh0,regh1,regl2) * v1(regh0,regh1,regl2) +
     &             v0(regh0,regh1,regh2) * v1(regh0,regh1,regh2))
      do i = regl0 + 1, regh0 - 1
         sum = sum + 0.25d0 *
     &               (v0(i,regl1,regl2) * v1(i,regl1,regl2) +
     &                v0(i,regl1,regh2) * v1(i,regl1,regh2) +
     &                v0(i,regh1,regl2) * v1(i,regh1,regl2) +
     &                v0(i,regh1,regh2) * v1(i,regh1,regh2))
      end do
      do j = regl1 + 1, regh1 - 1
          sum = sum + 0.25d0 *
     &               (v0(regl0,j,regl2) * v1(regl0,j,regl2) +
     &                v0(regl0,j,regh2) * v1(regl0,j,regh2) +
     &                v0(regh0,j,regl2) * v1(regh0,j,regl2) +
     &                v0(regh0,j,regh2) * v1(regh0,j,regh2))
      end do
      do k = regl2 + 1, regh2 - 1
          sum = sum + 0.25d0 *
     &               (v0(regl0,regl1,k) * v1(regl0,regl1,k) +
     &                v0(regl0,regh1,k) * v1(regl0,regh1,k) +
     &                v0(regh0,regl1,k) * v1(regh0,regl1,k) +
     &                v0(regh0,regh1,k) * v1(regh0,regh1,k))
      end do
      do j = regl1 + 1, regh1 - 1
         do i = regl0 + 1, regh0 - 1
             sum = sum + 0.5d0 *
     &                  (v0(i,j,regl2) * v1(i,j,regl2) +
     &                   v0(i,j,regh2) * v1(i,j,regh2))
         end do
      end do
      do k = regl2 + 1, regh2 - 1
         do i = regl0 + 1, regh0 - 1
            sum = sum + 0.5d0 *
     &                  (v0(i,regl1,k) * v1(i,regl1,k) +
     &                   v0(i,regh1,k) * v1(i,regh1,k))
         end do
      end do
      do k = regl2 + 1, regh2 - 1
         do j = regl1 + 1, regh1 - 1
             sum = sum + 0.5d0 *
     &                  (v0(regl0,j,k) * v1(regl0,j,k) +
     &                   v0(regh0,j,k) * v1(regh0,j,k))
         end do
      end do
      do  k = regl2 + 1, regh2 - 1
         do  j = regl1 + 1, regh1 - 1
            do  i = regl0 + 1, regh0 - 1
               sum = sum + v0(i,j,k) * v1(i,j,k)
            end do
         end do
      end do
      end
c-----------------------------------------------------------------------
c     These routines implement boundary conditions
c-----------------------------------------------------------------------
c Works for CELL- or NODE-based data.
      subroutine bref(
     & dest, destl0,desth0,destl1,desth1,destl2,desth2,
     &       regl0,regh0,regl1,regh1,regl2,regh2,
     & src,  srcl0,srch0,srcl1,srch1,srcl2,srch2,
     &       bbl0,bbh0,bbl1,bbh1,bbl2,bbh2, 
     & idir, ncomp)
      integer destl0,desth0,destl1,desth1,destl2,desth2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer bbl0,bbh0,bbl1,bbh1,bbl2,bbh2
      integer idir, ncomp
      double precision dest(destl0:desth0,destl1:desth1,destl2:desth2,
     &   ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2, ncomp)
      integer i, j, k, nc
      do nc = 1, ncomp
      if (idir .eq. 0) then
         do i = regl0, regh0
            do k = regl2, regh2
               do j = regl1, regh1
                  dest(i,j,k,nc) = src(bbh0-(i-regl0),j,k,nc)
               end do
            end do
         end do
      else if (idir .eq. 1) then
         do j = regl1, regh1
            do k = regl2, regh2
               do i = regl0, regh0
                  dest(i,j,k,nc) = src(i,bbh1-(j-regl1),k, nc)
               end do
            end do
         end do
      else
         do k = regl2, regh2
            do j = regl1, regh1
               do i = regl0, regh0
                  dest(i,j,k,nc) = src(i,j,bbh2-(k-regl2),nc)
               end do
            end do
         end do
      end if
      end do
      end
c-----------------------------------------------------------------------
c Works for CELL- or NODE-based data.
      subroutine brefm(
     & dest, destl0,desth0,destl1,desth1,destl2,desth2,
     &       regl0,regh0,regl1,regh1,regl2,regh2,
     & src,  srcl0,srch0,srcl1,srch1,srcl2,srch2,
     &       bbl0,bbh0,bbl1,bbh1,bbl2,bbh2, 
     & ra, ncomp)
      integer ncomp
      integer destl0,desth0,destl1,desth1,destl2,desth2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer bbl0,bbh0,bbl1,bbh1,bbl2,bbh2
      integer ra(0:2)
      double precision dest(destl0:desth0,destl1:desth1,destl2:desth2,
     &   ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2, ncomp)
      integer i, j, k, nc
      do nc = 1, ncomp
      if (ra(0) .eq. 0 .and. ra(1) .eq. 0 .and. ra(2) .eq. 0) then
            do k = regl2, regh2
               do j = regl1, regh1
                  do i = regl0, regh0
                     dest(i,j,k,nc) = src(bbl0+(i-regl0),
     &                                bbl1+(j-regl1), bbl2+(k-regl2),nc)
                  end do
               end do
            end do
      else if (ra(0) .eq. 0 .and. ra(1) .eq. 0 .and. ra(2) .eq. 1) then
         do k = regl2, regh2
            do j = regl1, regh1
               do i = regl0, regh0
                  dest(i,j,k,nc) = src(bbl0+(i-regl0),
     &                              bbl1+(j-regl1), bbh2-(k-regl2),nc)
               end do
            end do
         end do
      else if (ra(0) .eq. 0 .and. ra(1) .eq. 1 .and. ra(2) .eq. 0) then
         do j = regl1, regh1
            do k = regl2, regh2
               do i = regl0, regh0
                  dest(i,j,k,nc) = src(bbl0+(i-regl0),
     &                              bbh1-(j-regl1), bbl2+(k-regl2),nc)
               end do
            end do
         end do
      else if (ra(0) .eq. 0 .and. ra(1) .eq. 1 .and. ra(2) .eq. 1) then
         do k = regl2, regh2
            do j = regl1, regh1
               do i = regl0, regh0
                  dest(i,j,k,nc) = src(bbl0+(i-regl0),
     &                              bbh1-(j-regl1), bbh2-(k-regl2),nc)
               end do
            end do
         end do
      else if (ra(0) .eq. 1 .and. ra(1) .eq. 0 .and. ra(2) .eq. 0) then
         do i = regl0, regh0
            do k = regl2, regh2
               do j = regl1, regh1
                  dest(i,j,k,nc) = src(bbh0-(i-regl0),
     &                              bbl1+(j-regl1), bbl2+(k-regl2),nc)
               end do
            end do
         end do
      else if (ra(0) .eq. 1 .and. ra(1) .eq. 0 .and. ra(2) .eq. 1) then
         do k = regl2, regh2
            do i = regl0, regh0
               do j = regl1, regh1
                  dest(i,j,k,nc) = src(bbh0-(i-regl0),
     &                              bbl1+(j-regl1), bbh2-(k-regl2),nc)
               end do
            end do
         end do
      else if (ra(0) .eq. 1 .and. ra(1) .eq. 1 .and. ra(2) .eq. 0) then
         do j = regl1, regh1
            do i = regl0, regh0
               do k = regl2, regh2
                  dest(i,j,k,nc) = src(bbh0-(i-regl0),
     &                              bbh1-(j-regl1), bbl2+(k-regl2),nc)
               end do
            end do
         end do
      else
         do k = regl2, regh2
            do j = regl1, regh1
               do i = regl0, regh0
                  dest(i,j,k,nc) = src(bbh0-(i-regl0),
     &                              bbh1-(j-regl1), bbh2-(k-regl2),nc)
               end do
            end do
         end do
      end if
      end do
      end
c-----------------------------------------------------------------------
c Works for CELL- or NODE-based data.
      subroutine bneg(
     & dest, destl0,desth0,destl1,desth1,destl2,desth2,
     &       regl0,regh0,regl1,regh1,regl2,regh2,
     & src,  srcl0,srch0,srcl1,srch1,srcl2,srch2,
     &       bbl0,bbh0,bbl1,bbh1,bbl2,bbh2, 
     & idir, ncomp)
      integer destl0,desth0,destl1,desth1,destl2,desth2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer bbl0,bbh0,bbl1,bbh1,bbl2,bbh2
      integer idir,ncomp
      double precision dest(destl0:desth0,destl1:desth1,destl2:desth2,
     &   ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2, ncomp)
      integer i, j, k, nc
      do nc = 1, ncomp
      if (idir .eq. 0) then
         do i = regl0, regh0
            do k = regl2, regh2
               do j = regl1, regh1
                  dest(i,j,k,nc) = -src(bbh0-(i-regl0),j,k,nc)
               end do
            end do
         end do
      else if (idir .eq. 1) then
         do j = regl1, regh1
            do k = regl2, regh2
               do i = regl0, regh0
                  dest(i,j,k,nc) = -src(i,bbh1-(j-regl1),k,nc)
               end do
            end do
         end do
      else
         do k = regl2, regh2
            do j = regl1, regh1
               do i = regl0, regh0
                  dest(i,j,k,nc) = -src(i,j,bbh2-(k-regl2),nc)
               end do
            end do
         end do
      end if
      end do
      end
c-----------------------------------------------------------------------
c Works for CELL- or NODE-based data.
      subroutine bnegm(
     & dest, destl0,desth0,destl1,desth1,destl2,desth2,
     &       regl0,regh0,regl1,regh1,regl2,regh2,
     & src,  srcl0,srch0,srcl1,srch1,srcl2,srch2,
     &       bbl0,bbh0,bbl1,bbh1,bbl2,bbh2, 
     & ra, ncomp)
      integer ncomp
      integer destl0,desth0,destl1,desth1,destl2,desth2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer bbl0,bbh0,bbl1,bbh1,bbl2,bbh2
      integer ra(0:2)
      double precision dest(destl0:desth0,destl1:desth1,destl2:desth2, 
     &   ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2, ncomp)
      integer i, j, k, nc
      do nc = 1, ncomp
      if (ra(0) .eq. 0 .and. ra(1) .eq. 0 .and. ra(2) .eq. 0) then
         do k = regl2, regh2
            do j = regl1, regh1
               do i = regl0, regh0
                  dest(i,j,k,nc) = -src(bbl0+(i-regl0),
     &                            bbl1+(j-regl1), bbl2+(k-regl2),nc)
               end do
            end do
         end do
      else if (ra(0) .eq. 0 .and. ra(1) .eq. 0 .and. ra(2) .eq. 1) then
         do k = regl2, regh2
            do j = regl1, regh1
               do i = regl0, regh0
                  dest(i,j,k,nc) = -src(bbl0+(i-regl0),
     &                               bbl1+(j-regl1), bbh2-(k-regl2),nc)
               end do
            end do
         end do
      else if (ra(0) .eq. 0 .and. ra(1) .eq. 1 .and. ra(2) .eq. 0) then
         do j = regl1, regh1
            do k = regl2, regh2
               do i = regl0, regh0
                  dest(i,j,k,nc) = -src(bbl0+(i-regl0),
     &                               bbh1-(j-regl1), bbl2+(k-regl2),nc)
               end do
            end do
         end do
      else if (ra(0) .eq. 0 .and. ra(1) .eq. 1 .and. ra(2) .eq. 1) then
         do k = regl2, regh2
            do j = regl1, regh1
               do i = regl0, regh0
                  dest(i,j,k,nc) = -src(bbl0+(i-regl0),
     &                               bbh1-(j-regl1), bbh2-(k-regl2),nc)
               end do
            end do
         end do
      else if (ra(0) .eq. 1 .and. ra(1) .eq. 0 .and. ra(2) .eq. 0) then
         do i = regl0, regh0
            do k = regl2, regh2
               do j = regl1, regh1
                  dest(i,j,k,nc) = -src(bbh0-(i-regl0),
     &                               bbl1+(j-regl1), bbl2+(k-regl2),nc)
               end do
            end do
         end do
      else if (ra(0) .eq. 1 .and. ra(1) .eq. 0 .and. ra(2) .eq. 1) then
         do k = regl2, regh2
            do i = regl0, regh0
               do j = regl1, regh1
                  dest(i,j,k,nc) = -src(bbh0-(i-regl0),
     &                               bbl1+(j-regl1), bbh2-(k-regl2),nc)
               end do
            end do
         end do
      else if (ra(0) .eq. 1 .and. ra(1) .eq. 1 .and. ra(2) .eq. 0) then
         do j = regl1, regh1
            do i = regl0, regh0
               do k = regl2, regh2
                  dest(i,j,k,nc) = -src(bbh0-(i-regl0),
     &                               bbh1-(j-regl1), bbl2+(k-regl2),nc)
               end do
            end do
         end do
      else
         do k = regl2, regh2
            do j = regl1, regh1
               do i = regl0, regh0
                  dest(i,j,k,nc) = -src(bbh0-(i-regl0),
     &                               bbh1-(j-regl1), bbh2-(k-regl2),nc)
               end do
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
     & dest, destl0,desth0,destl1,desth1,destl2,desth2,
     &       regl0,regh0,regl1,regh1,regl2,regh2,
     & src,  srcl0,srch0,srcl1,srch1,srcl2,srch2,
     &       bbl0,bbh0,bbl1,bbh1,bbl2,bbh2,
     & idir, ncomp)
      integer ncomp
      integer destl0,desth0,destl1,desth1,destl2,desth2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer bbl0,bbh0,bbl1,bbh1,bbl2,bbh2
      integer idir
      double precision dest(destl0:desth0,destl1:desth1,destl2:desth2,
     &   ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2, ncomp)
      integer i, j, k, nc
      do nc = 1, ncomp
      if (idir .eq. 0) then
         if (regl0 .lt. bbh0) then
            do i = regl0, regh0
               do k = regl2, regh2
                  do j = regl1, regh1
                     dest(i,j,k,nc) = 2.d0 * dest(regh0,j,k,nc) -
     &                             src(bbh0-(i-regl0),j,k,nc)
                  end do
               end do
            end do
         else
            do i = regh0, regl0, -1
               do k = regl2, regh2
                  do j = regl1, regh1
                     dest(i,j,k,nc) = 2.d0 * dest(regh0,j,k, nc) -
     &                             src(bbh0-(i-regl0),j,k, nc)
                  end do
               end do
            end do
         end if
      else if (idir .eq. 1) then
         if (regl1 .lt. bbh1) then
            do j = regl1, regh1
               do k = regl2, regh2
                  do i = regl0, regh0
                     dest(i,j,k,nc) = 2.d0 * dest(i,regh1,k,nc) -
     &                             src(i,bbh1-(j-regl1),k,nc)
                  end do
               end do
            end do
         else
            do j = regh1, regl1, -1
               do k = regl2, regh2
                  do i = regl0, regh0
                     dest(i,j,k,nc) = 2.d0 * dest(i,regh1,k,nc) -
     &                             src(i,bbh1-(j-regl1),k,nc)
                  end do
               end do
            end do
         end if
      else
         if (regl2 .lt. bbh2) then
            do k = regl2, regh2
               do j = regl1, regh1
                  do i = regl0, regh0
                     dest(i,j,k,nc) = 2.d0 * dest(i,j,regh2,nc) -
     &                             src(i,j,bbh2-(k-regl2),nc)
                  end do
               end do
            end do
         else
            do k = regh2, regl2, -1
               do j = regl1, regh1
                  do i = regl0, regh0
                     dest(i,j,k,nc) = 2.d0 * dest(i,j,regh2,nc) -
     &                             src(i,j,bbh2-(k-regl2),nc)
                  end do
               end do
            end do
         end if
      end if
      end do
      end
c-----------------------------------------------------------------------
c Works for CELL-based velocity data.
c This routine is called when the inflow face velocity data has already
c been altered by a call to fill_borders.
      subroutine binfil(
     & dest, destl0,desth0,destl1,desth1,destl2,desth2,
     &       regl0,regh0,regl1,regh1,regl2,regh2,
     & src,  srcl0,srch0,srcl1,srch1,srcl2,srch2,
     &       bbl0,bbh0,bbl1,bbh1,bbl2,bbh2,
     & idir, ncomp)
      integer ncomp
      integer destl0,desth0,destl1,desth1,destl2,desth2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer bbl0,bbh0,bbl1,bbh1,bbl2,bbh2
      integer idir
      double precision dest(destl0:desth0,destl1:desth1,destl2:desth2,
     &   ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2, ncomp)
      integer i, j, k, nc
      do nc = 1, ncomp
      if (idir .eq. 0) then
         if (regl0 .lt. bbh0) then
            do  i = regl0, regh0
               do k = regl2, regh2
                  do j = regl1, regh1
                     dest(i,j,k,nc) = src(bbl0,j,k,nc) 
     &                           + src(bbl0+1,j,k,nc) -
     &                             src(bbh0-(i-regl0),j,k,nc)
                  end do
               end do
            end do
         else
            do i = regh0, regl0, -1
               do k = regl2, regh2
                  do j = regl1, regh1
                     dest(i,j,k,nc) = src(bbh0-1,j,k,nc) 
     &                + src(bbh0,j,k,nc) -
     &                             src(bbl0+(regh0-i),j,k,nc)
                  end do
               end do
            end do
         end if
      else if (idir .eq. 1) then
         if (regl1 .lt. bbh1) then
            do j = regl1, regh1
               do k = regl2, regh2
                  do i = regl0, regh0
                     dest(i,j,k,nc) = src(i,bbl1,k,nc)
     &                 + src(i,bbl1+1,k,nc) -
     &                             src(i,bbh1-(j-regl1),k,nc)
                  end do
               end do
            end do
         else
            do j = regh1, regl1, -1
               do k = regl2, regh2
                  do i = regl0, regh0
                     dest(i,j,k,nc) = src(i,bbh1-1,k,nc)
     &                 + src(i,bbh1,k,nc) -
     &                             src(i,bbl1+(regh1-j),k,nc)
                  end do
               end do
            end do
         end if
      else
         if (regl2 .lt. bbh2) then
            do k = regl2, regh2
               do j = regl1, regh1
                  do i = regl0, regh0
                     dest(i,j,k,nc) = src(i,j,bbl2,nc) 
     &                     + src(i,j,bbl2+1,nc) -
     &                             src(i,j,bbh2-(k-regl2),nc)
                  end do
               end do
            end do
         else
            do k = regh2, regl2, -1
               do j = regl1, regh1
                  do i = regl0, regh0
                     dest(i,j,k,nc) = src(i,j,bbh2-1,nc)
     &                      + src(i,j,bbh2,nc) -
     &                             src(i,j,bbl2+(regh2-k),nc)
                  end do
               end do
            end do
         end if
      end if
      end do
      end
c-----------------------------------------------------------------------
c     Interpolation routines
c-----------------------------------------------------------------------
c CELL-based data only.
      subroutine acint2(
     & dest, destl0,desth0,destl1,desth1,destl2,desth2,
     &       regl0,regh0,regl1,regh1,regl2,regh2,
     & src,  srcl0,srch0,srcl1,srch1,srcl2,srch2,
     &       bbl0,bbh0,bbl1,bbh1,bbl2,bbh2, 
     & ir, jr, kr, ncomp)
      integer ncomp
      integer destl0,desth0,destl1,desth1,destl2,desth2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer bbl0,bbh0,bbl1,bbh1,bbl2,bbh2
      integer ir, jr, kr
      double precision dest(destl0:desth0,destl1:desth1,destl2:desth2,
     &    ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2, ncomp)
      double precision xoff, yoff, zoff, sx, sy, sz
      integer ic, jc, kc, i, j, k, nc
      do nc = 1, ncomp
      do k = regl2, regh2
         kc = k/kr
         zoff = (mod(k,kr) + 0.5) / kr - 0.5
         do j = regl1, regh1
            jc = j/jr
            yoff = (mod(j,jr) + 0.5) / jr - 0.5
            do i = regl0, regh0
               ic = i/ir
               xoff = (mod(i,ir) + 0.5) / ir - 0.5
               sz = 0.5 * (src(ic,jc,kc+1,nc) - src(ic,jc,kc-1,nc))
               sy = 0.5 * (src(ic,jc+1,kc,nc) - src(ic,jc-1,kc,nc))
               sx = 0.5 * (src(ic+1,jc,kc,nc) - src(ic-1,jc,kc,nc))
               dest(i,j,k,nc) = src(ic,jc,kc,nc) 
     &          + xoff*sx + yoff*sy + zoff*sz
            end do
         end do
      end do
      end do
      end
c-----------------------------------------------------------------------
c NODE-based data only.
      subroutine anint2(
     & dest, destl0,desth0,destl1,desth1,destl2,desth2,
     &       regl0,regh0,regl1,regh1,regl2,regh2,
     & src,  srcl0,srch0,srcl1,srch1,srcl2,srch2,
     &       bbl0,bbh0,bbl1,bbh1,bbl2,bbh2, 
     & ir, jr, kr, ncomp)
      integer ncomp
      integer destl0,desth0,destl1,desth1,destl2,desth2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer bbl0,bbh0,bbl1,bbh1,bbl2,bbh2
      integer ir, jr, kr
      double precision dest(destl0:desth0,destl1:desth1,destl2:desth2,
     &   ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2, ncomp)
      double precision p, q
      integer ic, jc, kc, j, k, m, nc
      do nc = 1, ncomp
      do kc = bbl2, bbh2
         do jc = bbl1, bbh1
            do ic = bbl0, bbh0
               dest(ir*ic,jr*jc,kr*kc, nc) = src(ic,jc,kc, nc)
            end do
         end do
      end do
         do m = 1, kr-1
            q = dble(m)/kr
            p = 1.0d0 - q
            do kc = bbl2, bbh2-1
               do jc = bbl1, bbh1
                  do ic = bbl0, bbh0
                     dest(ir*ic,jr*jc,kr*kc+m,nc) =
     &                  p * src(ic,jc,kc,nc) + q * src(ic,jc,kc+1,nc)
                  end do
               end do
            end do
         end do
         do m = 1, jr-1
            q = dble(m)/jr
            p = 1.0d0 - q
            do jc = bbl1, bbh1-1
               do ic = bbl0, bbh0
cdir$ ivdep
                  do k = regl2, regh2
                     dest(ir*ic,jr*jc+m,k,nc) =
     &                  p * dest(ir*ic,jr*jc,k,nc) +
     &                  q * dest(ir*ic,jr*(jc+1),k,nc)
                  end do
               end do
            end do
         end do
            do  m = 1, ir-1
               q = dble(m)/ir
               p = 1.0d0 - q
               do ic = bbl0, bbh0-1
                  do k = regl2, regh2
cdir$ ivdep
                     do j = regl1, regh1
                        dest(ir*ic+m,j,k,nc) =
     &                     p * dest(ir*ic,j,k,nc) +
     &                     q * dest(ir*(ic+1),j,k,nc)
                     end do
                  end do
               end do
            end do
      end do
      end
c-----------------------------------------------------------------------
c     Restrictions
c-----------------------------------------------------------------------
c CELL-based data only.
      subroutine acrst1(
     & dest, destl0,desth0,destl1,desth1,destl2,desth2,
     &       regl0,regh0,regl1,regh1,regl2,regh2,
     & src,  srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & ir, jr, kr, ncomp, integ, i1, i2)
      integer ncomp
      integer destl0,desth0,destl1,desth1,destl2,desth2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer ir, jr, kr, integ, i1, i2
      double precision dest(destl0:desth0,destl1:desth1,destl2:desth2,
     &   ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2, ncomp)
      double precision fac
      integer i, j, k, l, m, n,nc
      do nc = 1, ncomp
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               dest(i,j,k,nc) = 0.d0
            end do
         end do
      end do
      do l = 0, kr-1
         do n = 0, jr-1
            do m = 0, ir-1
               do k = regl2, regh2
                  do j = regl1, regh1
                     do i = regl0, regh0
                        dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                    src(i*ir+m, j*jr+n, k*kr+l,nc)
                     end do
                  end do
               end do
            end do
         end do
      end do
      if (integ .eq. 0) then
         fac = 1.d0 / (ir*jr*kr)
         do k = regl2, regh2
            do j = regl1, regh1
               do i = regl0, regh0
                  dest(i,j,k,nc) = dest(i,j,k,nc) * fac
               end do
            end do
         end do
      end if
      end do
      end
c-----------------------------------------------------------------------
c NODE-based data only.
      subroutine anrst1(
     & dest, destl0,desth0,destl1,desth1,destl2,desth2,
     &       regl0,regh0,regl1,regh1,regl2,regh2,
     & src,  srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & ir, jr, kr, ncomp, integ, i1, i2)
      integer ncomp, integ
      integer destl0,desth0,destl1,desth1,destl2,desth2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer ir, jr, kr
      integer nc
      double precision dest(destl0:desth0,destl1:desth1,destl2:desth2,
     &   ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2, ncomp)
      integer i, j, k, i1, i2
      do nc = 1, ncomp
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               dest(i,j,k,nc) = src(i*ir, j*jr, k*kr, nc)
            end do
         end do
      end do
      end do
      end
c-----------------------------------------------------------------------
c NODE-based data only.
      subroutine anrst2(
     & dest, destl0,desth0,destl1,desth1,destl2,desth2,
     &       regl0,regh0,regl1,regh1,regl2,regh2,
     & src,  srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & ir, jr, kr, ncomp, integ, i1, i2)
      integer ncomp
      integer destl0,desth0,destl1,desth1,destl2,desth2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer ir, jr, kr, integ, i1, i2
      double precision dest(destl0:desth0,destl1:desth1,destl2:desth2,
     &   ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2, ncomp)
      double precision fac0, fac1, fac2, fac
      integer i, j, k, l, m, n, nc
      integer r,s
      do nc = 1, ncomp
         do k = regl2, regh2
            do j = regl1, regh1
               do i = regl0, regh0
                  dest(i,j,k,nc) = 0.d0
               end do
            end do
         end do
         fac0 = 1.d0 / (ir*ir * jr*jr * kr*kr)
         do l = 0, kr-1
            fac2 = (kr-l) * fac0
            if (l .eq. 0) fac2 = 0.5d0 * fac2
            do n = 0, jr-1
               fac1 = (jr-n) * fac2
               if (n .eq. 0) fac1 = 0.5d0 * fac1
               do m = 0, ir-1
                  fac = (ir-m) * fac1
                  if (m .eq. 0) fac = 0.5d0 * fac
                  do k = regl2, regh2
                     do j = regl1, regh1
                        do i = regl0, regh0
                           dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                  fac * (src(i*ir-m,j*jr-n,k*kr-l,nc)+
     &                         src(i*ir-m,j*jr-n,k*kr+l,nc)+
     &                         src(i*ir-m,j*jr+n,k*kr-l,nc)+
     &                         src(i*ir-m,j*jr+n,k*kr+l,nc)+
     &                         src(i*ir+m,j*jr-n,k*kr-l,nc)+
     &                         src(i*ir+m,j*jr-n,k*kr+l,nc)+
     &                         src(i*ir+m,j*jr+n,k*kr-l,nc)+
     &                         src(i*ir+m,j*jr+n,k*kr+l,nc))
                        end do
                     end do
                  end do
               end do
            end do
         end do
         if (integ .eq. 1) then
            fac = ir * jr * kr
            do k = regl2, regh2
               do j = regl1, regh1
                  do i = regl0, regh0
                     dest(i,j,k,nc) = fac * dest(i,j,k,nc)
                  end do
                end do
            end do
         end if
      end do
      end
c-----------------------------------------------------------------------
c NODE-based data only.
c Fills coarse region defined by reg.
c Handles coarse-fine face, with orientation determined by idim and idir
      subroutine anfr2(
     & dest, destl0,desth0,destl1,desth1,destl2,desth2,
     &       regl0,regh0,regl1,regh1,regl2,regh2,
     & src,  srcl0,srch0,srcl1,srch1,srcl2,srch2,
     & ir, jr, kr, ncomp, integ, idim, idir)
      integer destl0,desth0,destl1,desth1,destl2,desth2
      integer regl0,regh0,regl1,regh1,regl2,regh2
      integer srcl0,srch0,srcl1,srch1,srcl2,srch2
      integer ir, jr, kr, idim, idir, integ
      integer ncomp
      double precision dest(destl0:desth0,destl1:desth1,destl2:desth2,
     &   ncomp)
      double precision src(srcl0:srch0,srcl1:srch1,srcl2:srch2, ncomp)
      double precision fac0, fac1, fac2, fac
      integer i, j, k, l, m, n, nc
      if (idim .eq. 0) then
         if (integ .eq. 0) then
            fac = (0.5d0 + 0.5d0 / ir)
            fac0 = 1.d0 / (ir*ir * jr*jr * kr*kr)
         else
            fac = 1.d0
            fac0 = 1.d0 / (ir * jr * kr)
         end if
         i = regl0
         do nc = 1, ncomp
         do k = regl2, regh2
            do j = regl1, regh1
               dest(i,j,k,nc) = fac * src(i*ir,j*jr,k*kr,nc)
            end do
         end do
         end do
         do nc = 1, ncomp
         do l = 0, kr-1
            fac2 = (kr-l) * fac0
            if (l .eq. 0) fac2 = 0.5d0 * fac2
            do n = 0, jr-1
               fac1 = (jr-n) * fac2
               if (n .eq. 0) fac1 = 0.5d0 * fac1
               do m = idir, idir*(ir-1), idir
                  fac = (ir-abs(m)) * fac1
                  do k = regl2, regh2
                     do j = regl1, regh1
                        dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                       fac * (src(i*ir+m,j*jr-n,k*kr-l,nc) +
     &                              src(i*ir+m,j*jr-n,k*kr+l,nc) +
     &                              src(i*ir+m,j*jr+n,k*kr-l,nc) +
     &                              src(i*ir+m,j*jr+n,k*kr+l,nc))
                     end do
                  end do
               end do
            end do
         end do
         end do
      else if (idim .eq. 1) then
         if (integ .eq. 0) then
            fac = (0.5d0 + 0.5d0 / jr)
            fac0 = 1.d0 / (ir*ir * jr*jr * kr*kr)
         else
            fac = 1.d0
            fac0 = 1.d0 / (ir * jr * kr)
         end if
         j = regl1
         do nc = 1, ncomp
         do k = regl2, regh2
            do i = regl0, regh0
               dest(i,j,k,nc) = fac * src(i*ir,j*jr,k*kr,nc)
            end do
         end do
         end do
         fac0 = 1.d0 / (ir*ir * jr*jr * kr*kr)
         do nc = 1, ncomp
         do l = 0, kr-1
            fac2 = (kr-l) * fac0
            if (l .eq. 0) fac2 = 0.5d0 * fac2
            do n = idir, idir*(jr-1), idir
               fac1 = (jr-abs(n)) * fac2
               do m = 0, ir-1
                  fac = (ir-m) * fac1
                  if (m .eq. 0) fac = 0.5d0 * fac
                  do k = regl2, regh2
                     do i = regl0, regh0
                        dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                       fac * (src(i*ir-m,j*jr+n,k*kr-l,nc) +
     &                              src(i*ir-m,j*jr+n,k*kr+l,nc) +
     &                              src(i*ir+m,j*jr+n,k*kr-l,nc) +
     &                              src(i*ir+m,j*jr+n,k*kr+l,nc))
                     end do
                  end do
               end do
            end do
         end do
         end do
      else
         if (integ .eq. 0) then
            fac = (0.5d0 + 0.5d0 / kr)
            fac0 = 1.d0 / (ir*ir * jr*jr * kr*kr)
         else
            fac = 1.d0
            fac0 = 1.d0 / (ir * jr * kr)
         end if
         k = regl2
         do nc = 1, ncomp
         do j = regl1, regh1
            do i = regl0, regh0
               dest(i,j,k,nc) = fac * src(i*ir,j*jr,k*kr,nc)
            end do
         end do
         end do
         fac0 = 1.d0 / (ir*ir * jr*jr * kr*kr)
         do nc = 1, ncomp
         do l = idir, idir*(kr-1), idir
            fac2 = (kr-abs(l)) * fac0
            do n = 0, jr-1
               fac1 = (jr-n) * fac2
               if (n .eq. 0) fac1 = 0.5d0 * fac1
               do m = 0, ir-1
                  fac = (ir-m) * fac1
                  if (m .eq. 0) fac = 0.5d0 * fac
                  do j = regl1, regh1
                     do i = regl0, regh0
                        dest(i,j,k,nc) = dest(i,j,k,nc) +
     &                       fac * (src(i*ir-m,j*jr-n,k*kr+l,nc) +
     &                              src(i*ir-m,j*jr+n,k*kr+l,nc) +
     &                              src(i*ir+m,j*jr-n,k*kr+l,nc) +
     &                              src(i*ir+m,j*jr+n,k*kr+l,nc))
                     end do
                  end do
               end do
            end do
         end do
         end do
      end if
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
            dest(i,j,k,nc) = 0.d0
         end do
c center gets center plus two edges
         if (integ .eq. 0) then
            center = ir / cube
            fac0 = 1.d0 / (cube**2)
            cfac = 0.25d0 * cube * fac0 * ir * (jr-1) * (kr-1)
         else
            center = 1.d0
            fac0 = 1.d0 / cube
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
                           if (m .eq. 0) fac = 0.5d0 * fac
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
         cfac = 0.5d0 * cube * fac0 * ir * (jr-1)
         do ji = 0, 1
            jdir = 2 * ji - 1
            if (ga(0,ji,0) + ga(0,ji,1) .eq. 2) then
               do n = jdir, jdir*(jr-1), jdir
                  fac1 = (jr-abs(n)) * fac2
                  do m = 0, ir-1
                     fac = (ir-m) * fac1
                     if (m .eq. 0) fac = 0.5d0 * fac
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
         cfac = 0.5d0 * cube * fac0 * ir * (kr-1)
         do ki = 0, 1
            kdir = 2 * ki - 1
            if (ga(0,0,ki) + ga(0,1,ki) .eq. 2) then
               do l = kdir, kdir*(kr-1), kdir
                  fac1 = (kr-abs(l)) * fac2
                  do m = 0, ir-1
                     fac = (ir-m) * fac1
                     if (m .eq. 0) fac = 0.5d0 * fac
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
            dest(i,j,k,nc) = 0.d0
         end do
c center gets center plus two edges
         if (integ .eq. 0) then
            center = jr / cube
            fac0 = 1.d0 / (cube**2)
            cfac = 0.25d0 * cube * fac0 * jr * (ir-1) * (kr-1)
         else
            center = 1.d0
            fac0 = 1.d0 / cube
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
                        if (n .eq. 0) fac1 = 0.5d0 * fac1
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
         cfac = 0.5d0 * cube * fac0 * jr * (ir-1)
         do ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,0,0) + ga(ii,0,1) .eq. 2) then
               do n = 0, jr-1
                  fac1 = (jr-n) * fac2
                  if (n .eq. 0) fac1 = 0.5d0 * fac1
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
         cfac = 0.5d0 * cube * fac0 * jr * (kr-1)
         do ki = 0, 1
            kdir = 2 * ki - 1
            if (ga(0,0,ki) + ga(1,0,ki) .eq. 2) then
               do  l = kdir, kdir*(kr-1), kdir
                  fac1 = (kr-abs(l)) * fac2
                  do n = 0, jr-1
                     fac = (jr-n) * fac1
                     if (n .eq. 0) fac = 0.5d0 * fac
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
            dest(i,j,k,nc) = 0.d0
         end do
c center gets center plus two edges
         if (integ .eq. 0) then
            center = kr / cube
            fac0 = 1.d0 / (cube**2)
            cfac = 0.25d0 * cube * fac0 * kr * (ir-1) * (jr-1)
         else
            center = 1.d0
            fac0 = 1.d0 / cube
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
                     if (l .eq. 0) fac2 = 0.5d0 * fac2
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
         cfac = 0.5d0 * cube * fac0 * kr * (ir-1)
         do ii = 0, 1
            idir = 2 * ii - 1
            if (ga(ii,0,0) + ga(ii,1,0) .eq. 2) then
               do l = 0, kr-1
                  fac1 = (kr-l) * fac2
                  if (l .eq. 0) fac1 = 0.5d0 * fac1
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
         cfac = 0.5d0 * cube * fac0 * kr * (jr-1)
         do ji = 0, 1
            jdir = 2 * ji - 1
            if (ga(0,ji,0) + ga(1,ji,0) .eq. 2) then
               do l = 0, kr-1
                  fac1 = (kr-l) * fac2
                  if (l .eq. 0) fac1 = 0.5d0 * fac1
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
      dest(i,j,k,nc) = 0.d0
      cube = ir * jr * kr
      if (integ .eq. 0) then
         center = 1.d0 / cube
         fac0 = 1.d0 / (cube**2)
         cfac = 0.125d0 * cube * fac0 * (ir-1) * (jr-1) * (kr-1)
      else
         center = 1.d0
         fac0 = 1.d0 / cube
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
      cfac = 0.25d0 * cube * fac0 * (ir-1) * (jr-1)
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
      cfac = 0.25d0 * cube * fac0 * (ir-1) * (kr-1)
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
      cfac = 0.25d0 * cube * fac0 * (jr-1) * (kr-1)
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
      cfac = 0.5d0 * cube * fac0 * (ir-1)
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
      cfac = 0.5d0 * cube * fac0 * (jr-1)
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
      cfac = 0.5d0 * cube * fac0 * (kr-1)
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
