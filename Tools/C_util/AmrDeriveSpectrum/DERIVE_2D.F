
#include "CONSTANTS.H"
#include "Derived.H"

#define DIM 2

c
c The following routine does NOT have the signature of a derived function.
c
c This routine assumes domain is of form: [(0,0),(xlen,ylen)]
c
       subroutine FORT_RANDMOVE (a,DIMS(c),b,c,d,u,nspec,x,y,spec,
     &                           dx,dy,dt,dxinv,dyinv)

       implicit none

       integer DIMDEC(c), nspec, spec
       REAL_T  x, y, dx, dy, dt, dxinv, dyinv
       REAL_T  a(DIMV(c),0:nspec-1)
       REAL_T  b(DIMV(c),0:nspec-1)
       REAL_T  c(DIMV(c),0:nspec-1)
       REAL_T  d(DIMV(c),0:nspec-1)
       REAL_T  u(DIMV(c),2)

       integer i, ix, iy
       REAL_T  rn

       if (spec.lt.zero.or.spec.ge.nspec) stop 'bad species ref (df)'
       !
       ! First advect stuff ...
       !
       ix = x * dxinv
       iy = y * dyinv

       if ((ix.ge.c_l1.and.ix.le.c_h1) .and.
     &      (iy.ge.c_l2.and.iy.le.c_h2)) then

          x = x + u(ix,iy,1) * dt
          y = y + u(ix,iy,2) * dt

          if (x.lt.zero) x = -x

       end if
       !
       ! Now do the stochastic diffusion ...
       !
       ix = x * dxinv
       iy = y * dyinv

       if ((ix.ge.c_l1.and.ix.le.c_h1) .and.
     &      (iy.ge.c_l2.and.iy.le.c_h2)) then

          call blutilrand(rn)

          if (rn.lt.a(ix,iy,spec)*dt) then
             x = x - dx
          else if (rn.lt.b(ix,iy,spec)*dt) then
             x = x + dx
          else if (rn.lt.c(ix,iy,spec)*dt) then
             y = y - dy
          else if (rn.lt.d(ix,iy,spec)*dt) then
             y = y + dy
          end if

       end if

       end

       subroutine FORT_MAXDTDIFF (d,DIMS(c),nspec,x,y,spec,
     &                            dx,dy,dt,dxinv,dyinv)

       implicit none

       integer DIMDEC(c), nspec, spec
       REAL_T  x, y, dx, dy, dt, dxinv, dyinv
       REAL_T  d(DIMV(c),0:nspec-1)

       integer i, j, ilo, ihi, jlo, jhi, ix, iy
       REAL_T  rn, dmax

       if (spec.lt.zero.or.spec.ge.nspec) stop 'bad species ref (df)'

       ix   = x * dxinv
       iy   = y * dyinv
       dmax = 0

       if ((ix.ge.c_l1.and.ix.le.c_h1) .and.
     &      (iy.ge.c_l2.and.iy.le.c_h2)) then

          ilo = max(ix-2,c_l1)
          ihi = min(ix+2,c_h1)
          jlo = max(iy-2,c_l2)
          jhi = min(iy+2,c_h2)

          do i = ilo, ihi
             do j = jlo, jhi
                dmax = max(dmax,d(i,j,spec))
             end do
          end do

       end if

       if (dmax.gt.zero) dt = min(dt,.1/dmax)

       end
c
c This routine assumes domain is of form: [(0,0),(xlen,ylen)]
c
c conc is really conc.invert(1) so we can multiply instead of divide ...
c
       subroutine FORT_CHEMMOVE (rf,DIMS(c),rr,conc,nspec,nreac,
     &                           edges,edgeslen,pedges,x,y,spec,
     &                           dx,dy,dt,rxn,dxinv,dyinv)

       implicit none

       integer DIMDEC(c), nspec, nreac, spec, rxn
       integer edgeslen, edges(0:edgeslen-1), pedges(0:nspec-1)
       REAL_T  x, y, dx, dy, dt, dxinv, dyinv
       REAL_T  rf(DIMV(c),0:nreac-1)
       REAL_T  rr(DIMV(c),0:nreac-1)
       REAL_T  conc(DIMV(c),0:nspec-1)

       integer i, ix, iy, ie, nedges, rxnid, factor, tospec, nu
       REAL_T  rn, lambda, netreact

       if (spec.lt.zero.or.spec.ge.nspec) stop 'bad species ref (ch)'

       ix = x * dxinv
       iy = y * dyinv

       if ((ix.ge.c_l1.and.ix.le.c_h1) .and.
     &      (iy.ge.c_l2.and.iy.le.c_h2)) then

          call blutilrand(rn)

          if (pedges(spec).lt.zero) stop 'invalid reaction ref (ch)'

          lambda = 0
          nedges = edges(pedges(spec))

          if (nedges.le.zero) stop 'nedges must be positive (ch)'

          do ie = 0, nedges-1

             rxnid  = edges(pedges(spec)+ie*4+1)
             factor = edges(pedges(spec)+ie*4+2)
             tospec = edges(pedges(spec)+ie*4+3)
             nu     = edges(pedges(spec)+ie*4+4)

             if (factor.gt.zero) then
c                netreact = max(rf(ix,iy,rxnid)-rr(ix,iy,rxnid),zero)
                netreact = rf(ix,iy,rxnid)
                lambda   = nu*lambda+netreact*conc(ix,iy,spec)*dt/factor
             else
c                netreact = max(rr(ix,iy,rxnid)-rf(ix,iy,rxnid),zero)
                netreact = rr(ix,iy,rxnid)
                lambda   = nu*lambda-netreact*conc(ix,iy,spec)*dt/factor
             end if

             if (rn.lt.lambda) then
                rxn  = rxnid
                spec = tospec
                goto 100
             end if

          end do

       end if

100    continue

       end

       subroutine FORT_MAXLAMBDA (rf,DIMS(c),rr,conc,nspec,nreac,
     &                            edges,edgeslen,pedges)

       implicit none

       integer DIMDEC(c), nspec, nreac
       integer edgeslen, edges(0:edgeslen-1), pedges(0:nspec-1)
       REAL_T  rf(DIMV(c),0:nreac-1)
       REAL_T  rr(DIMV(c),0:nreac-1)
       REAL_T  conc(DIMV(c),0:nspec-1)

       integer i, ispec, ix, iy, ie, nedges, rxnid, factor, tospec, nu
       REAL_T  rn, lambda, lmax, netreact

       do ispec = 0, nspec-1

          lmax   = 0
          nedges = edges(pedges(ispec))

          do ix = c_l1, c_h1
             do iy = c_l2, c_h2

                lambda = 0

                do ie = 0, nedges-1
                   rxnid  = edges(pedges(ispec)+ie*4+1)
                   factor = edges(pedges(ispec)+ie*4+2)
                   tospec = edges(pedges(ispec)+ie*4+3)
                   nu     = edges(pedges(ispec)+ie*4+4)

                   if (factor.gt.zero) then
c                      netreact = max(rf(ix,iy,rxnid)-rr(ix,iy,rxnid),zero)
                      netreact = rf(ix,iy,rxnid)
                      lambda   = nu*lambda+netreact*conc(ix,iy,ispec)/factor
                   else
c                      netreact = max(rr(ix,iy,rxnid)-rf(ix,iy,rxnid),zero)
                      netreact = rr(ix,iy,rxnid)
                      lambda   = nu*lambda-netreact*conc(ix,iy,ispec)/factor
                   end if
                end do

                lmax = max(lmax,lambda)

             end do
          end do

          write(*,*) 'spec = ', ispec, 'lmax = ', lmax

       end do

       end

       subroutine FORT_MAXDTLAMBDA (rf,DIMS(c),rr,conc,nspec,nreac,
     &                              edges,edgeslen,pedges,x,y,spec,
     &                              dx,dy,dt,dxinv,dyinv)

       implicit none

       integer DIMDEC(c), nspec, nreac, spec
       integer edgeslen, edges(0:edgeslen-1), pedges(0:nspec-1)
       REAL_T  x, y, dx, dy, dt, dxinv, dyinv
       REAL_T  rf(DIMV(c),0:nreac-1)
       REAL_T  rr(DIMV(c),0:nreac-1)
       REAL_T  conc(DIMV(c),0:nspec-1)

       integer i, j, ix, iy, ie, nedges, rxnid, factor, tospec, nu
       integer ilo, ihi, jlo, jhi, ipos
       REAL_T  lambda, lmax, netreact

       if (spec.lt.zero.or.spec.ge.nspec) stop 'bad species ref (ch)'

       ix   = x * dxinv
       iy   = y * dyinv
       lmax = 0

       if (pedges(spec).lt.zero) stop 'invalid reaction ref (ch)'

       lambda = 0
       nedges = edges(pedges(spec))

       if (nedges.le.zero) stop 'nedges must be positive (ch)'

       do ie = 0, nedges-1
          rxnid  = edges(pedges(spec)+ie*4+1)
          factor = edges(pedges(spec)+ie*4+2)
          tospec = edges(pedges(spec)+ie*4+3)
          nu     = edges(pedges(spec)+ie*4+4)

          if (factor.gt.zero) then
c             netreact = max(rf(ix,iy,rxnid)-rr(ix,iy,rxnid),zero)
             netreact = rf(ix,iy,rxnid)
             lambda   = nu*lambda+netreact*conc(ix,iy,spec)/factor
          else
c             netreact = max(rr(ix,iy,rxnid)-rf(ix,iy,rxnid),zero)
             netreact = rr(ix,iy,rxnid)
             lambda   = nu*lambda-netreact*conc(ix,iy,spec)/factor
          end if
       end do

       lmax = max(lmax,lambda)

       if (lmax.gt.zero) dt = min(dt,.1/lmax)

       end

      subroutine FORT_PROBFAB (rf,DIMS(c),rr,prob,DIMS(p),
     &                         nspec,nreac,ispec,
     &                         edges,edgeslen,pedges,isrz,dx,dy)

      implicit none

      integer DIMDEC(c), DIMDEC(p), nspec, nreac, ispec, isrz
      integer edgeslen, edges(0:edgeslen-1), pedges(0:nspec-1)
      REAL_T  rf(DIMV(c),0:nreac-1)
      REAL_T  rr(DIMV(c),0:nreac-1)
      REAL_T  prob(DIMV(p)), dx, dy

      integer ix, iy, ie, nedges, rxnid, factor, nu
      REAL_T  lambda, lmax, vol

      if (ispec.lt.zero.or.ispec.gt.nspec-1) stop 'invalid ispec'

      nedges = edges(pedges(ispec))

      do iy = p_l2, p_h2
         do ix = p_l1, p_h1

            vol = dx*dy

            if (isrz.ne.0) vol = vol*half*Pi*dx*(2*p_l1+1)

            lambda = 0

            do ie = 0, nedges-1

               rxnid  = edges(pedges(ispec)+ie*4+1)
               factor = edges(pedges(ispec)+ie*4+2)
               nu     = edges(pedges(ispec)+ie*4+4)

               if (factor.gt.zero) then
                  lambda = nu*lambda + rf(ix,iy,rxnid)/factor
               else
                  lambda = nu*lambda - rr(ix,iy,rxnid)/factor
               end if

            end do

            prob(ix,iy) = lambda*vol

         end do
      end do

      end

      subroutine FORT_ACCUMPROB (prob, DIMS(p), totreact, cumprob)

      implicit none

      integer DIMDEC(p)
      REAL_T  prob(DIMV(p)), totreact, cumprob

      integer ix, iy

      do iy = p_l2, p_h2
         do ix = p_l1, p_h1

            prob(ix,iy) = cumprob+prob(ix,iy)/totreact
            cumprob     = prob(ix,iy)

         enddo
      enddo

      end

      subroutine FORT_SELECTPOS (prob,DIMS(c),rn,dx,dy,x,y)

      implicit none

      integer DIMDEC(c)
      REAL_T  prob(DIMV(c))
      REAL_T  dx,dy,x,y,rn

      integer ix, iy

      do iy = c_l2, c_h2-1
         if (rn.le.prob(c_h1,iy)) go to 100
      enddo

 100  continue

      do ix = c_l1, c_h1-1
         if (rn.le.prob(ix,iy)) go to 200
      enddo

 200  continue

      x = (ix+half)*dx
      y = (iy+half)*dy

      end
c
c     INTEGRATE (2D version)
c   
c     integral = sum { field(i,j) * dx * dy }
c   
c     Inputs / Outputs:
c   
c     field       => field of values to integrate
c     DIMS(field) => index limits of field
c     DIMS(valid) => index limits of valid data in field
c     delta       => cell size
c     integral    <= integral
c
      subroutine FORT_INTEGRATE 
     +   (field, DIMS(field), DIMS(valid), delta, integral)
      integer DIMDEC(field)
      integer DIMDEC(valid)
      REAL_T  delta(2), integral
      REAL_T  field(DIMV(field))

      integer i, j
      REAL_T  sum

      sum = 0

      do j = ARG_L2(valid), ARG_H2(valid)
         do i = ARG_L1(valid), ARG_H1(valid)
            sum = sum + field(i, j)
c            if (field (i, j) .ne. 0) sum = sum + 1
         end do
      end do

      integral = delta(1) * delta(2) * sum

      end
c
c     JPDF (2D version)
c   
c     Inputs / Outputs:
c   
c     field       => field of values to bin (1 = horizontal, 2 = vertical)
c     visible     => flags indicating visibility of cells (0.0 no, 1.0 yes)
c     DIMS(field) => index limits of field
c     DIMS(valid) => index limits of visible, and of valid data in field
c     hBins       => number of horizontal bins
c     hMax        => upper limit for horizontal bins
c     hMin        => lower limit for horizontal bins
c     vBins       => number of vertical bins
c     vMax        => upper limit for vertical bins
c     vMin        => lower limit for vertical bins
c     excessCount <= count of cells outside the bin limits
c     jointCount  <= count of cells inside the individual, two-way bins
c
      subroutine FORT_JPDF
     +     (field, visible, DIMS(field), DIMS(valid), 
     +     hBins, hMax, hMin, vBins, vMax, vMin,
     +     excessCount, jointCount)
      double precision  foobar
      integer DIMDEC(field)
      integer DIMDEC(valid)
      integer hBins, vBins
      integer excessCount
      integer jointCount
      REAL_T  hMax, hMin, vMax, vMin
      REAL_T  field(DIMV(field), 2)
c     the visible box is assumed sized to the valid region
      REAL_T  visible(DIMV(valid))
      dimension jointCount (0 : hBins * vBins - 1)

      integer bin, hBin, i, j, vBin
      REAL_T h1, h2, hVal, v1, v2, vVal

      h1 = (hBins - half) / (hMax - hMin)
      h2 = (- half) / (hMax - hMin)
      v1 = (vBins - half) / (vMax - vMin)
      v2 = (- half) / (vMax - vMin)

      do j = ARG_L2(valid), ARG_H2(valid)
         do i = ARG_L1(valid), ARG_H1(valid)
            if (visible (i, j) .eq. 1) then
               hVal = field(i, j, 1)
               if (hVal .eq. hMin) then
                  hBin = 0
               else if (hVal .eq. hMax) then
                  hBin = hBins - 1
               else
c                 if either of the Min and Max limits is an extrema, then
c                 rounding may choose an out of range bin when none is 
c                 intended.  hence the need to deal with this case, above. 
                  hBin = nint ((hVal - hMin) * h1 + (hMax - hVal) * h2)
               end if

               vVal = field(i, j, 2)
               if (vVal .eq. vMin) then
                  vBin = 0
               else if (vVal .eq. vMax) then
                  vBin = vBins - 1
               else
c                 if either of the Min and Max limits is an extrema, then
c                 rounding may choose an out of range bin when none is 
c                 intended.  hence the need to deal with this case, above. 
                  vBin = nint ((vVal - vMin) * v1 + (vMax - vVal) * v2)
               end if

               if (0 .le. hBin .and. hBin .lt. hBins .and.
     +             0 .le. vBin .and. vBin .lt. vBins) then
                  bin = hBin + hBins * vBin
                  jointCount(bin) = jointCount(bin) + 1
               else
                  excessCount = excessCount + 1
               end if
            end if
         end do
      end do

      return
      end
c
c     LIMITS (2D version)
c   
c     maximum = max { field(i,j) }
c     minimum = min { field(i,j) }
c   
c     Inputs / Outputs:
c   
c     field       => field of values to integrate
c     visible     => flags indicating visibility of cells (0.0 no, 1.0 yes)
c     DIMS(field) => index limits of field
c     DIMS(valid) => index limits of visible, and of valid data in field
c     maximum     <= max of starting value and field cells that are visible
c     minimum     <= min of starting value and field cells that are visible
c
      subroutine FORT_LIMITS
     +   (field, visible, DIMS(field), DIMS(valid), maximum, minimum)
      integer DIMDEC(field)
      integer DIMDEC(valid)
      integer ifound
      REAL_T  maximum, minimum
      REAL_T  field(DIMV(field))
c     the visible box is assumed sized to the valid region
      REAL_T  visible(DIMV(valid))

      integer i, j

      do j = ARG_L2(valid), ARG_H2(valid)
         do i = ARG_L1(valid), ARG_H1(valid)
            if (visible (i, j) .eq. 1) then
               maximum = max (maximum, field(i, j))
               minimum = min (minimum, field(i, j))
            end if
         end do
      end do

      end
c :: ----------------------------------------------------------
c :: SUMMASS
c ::             MASS = sum{ vol(i,j)*rho(i,j) }
c ::
c :: INPUTS / OUTPUTS:
c ::  rho        => density field
c ::  DIMS(rho)  => index limits of rho aray
c ::  lo,hi      => index limits of grid interior
c ::  dx         => cell size
c ::  mass      <=  total mass
c ::  r          => radius at cell center
c ::  irlo,hi    => index limits of r array
c ::  rz_flag    => == 1 if R_Z coords
c ::  tmp        => temp column array
c :: ----------------------------------------------------------
c ::
      subroutine FORT_SUMMASS(rho,DIMS(rho),DIMS(grid),dx,mass,
     &     r,irlo,irhi,rz_flag,tmp)

      integer irlo, irhi, rz_flag
      integer DIMDEC(rho)
      integer DIMDEC(grid)
      REAL_T  mass, dx(2)
      REAL_T  rho(DIMV(rho))
      REAL_T  r(irlo:irhi)
      REAL_T  tmp(DIM2(grid))
      
      integer i, j
      REAL_T  dr, dz, vol
      
      dr = dx(1)
      dz = dx(2)
      
      do j = ARG_L2(grid), ARG_H2(grid)
         tmp(j) = zero
      end do
      
      do i = ARG_L1(grid), ARG_H1(grid)
         vol = dr*dz
         if (rz_flag .eq. 1) vol = vol*two*Pi*r(i)
         do j = ARG_L2(grid), ARG_H2(grid)
            tmp(j) = tmp(j) + vol*rho(i,j)
         end do
      end do
      
      mass = zero
      do j = ARG_L2(grid), ARG_H2(grid)
         mass = mass + tmp(j)
      end do
      
      end
      
