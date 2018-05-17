#include <AMReX_REAL.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MCLO_F.H>
#include <AMReX_ArrayLim.H>

c-----------------------------------------------------------------------
      subroutine FORT_HARMONIC_AVERAGEEC (
     $     c, DIMS(c),
     $     f, DIMS(f),
     $     lo, hi, nc,
     $     cdir
     $     )
c
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cdir
      integer DIMDEC(f)
      REAL_T f(DIMV(f),nc)
      integer DIMDEC(c)
      REAL_T c(DIMV(c),nc)
c
      REAL_T factor, den
      parameter(factor=2.00)
      integer n
      integer i
      integer j
c
      if ( cdir .eq. 0 ) then
         do n = 1, nc
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  den = f(2*i,2*j,n) + f(2*i,2*j+1,n)
                  if (den .ne. 0.0) then
                    c(i,j,n) =  factor*f(2*i,2*j,n)*f(2*i,2*j+1,n)/den
                  else
                    c(i,j,n) =  0.0
                  endif
               enddo
            enddo
         enddo
      else if (cdir .eq. 1 ) then
         do n = 1, nc
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  den = f(2*i,2*j,n) + f(2*i+1,2*j,n)
                  if (den .ne. 0.0) then
                    c(i,j,n) =  factor*f(2*i,2*j,n)*f(2*i+1,2*j,n)/den
                  else
                    c(i,j,n) =  0.0
                  endif
               enddo
            enddo
         enddo
      endif
c
      end
c-----------------------------------------------------------------------
      subroutine FORT_AVERAGEEC (
     $     c, DIMS(c),
     $     f, DIMS(f),
     $     lo, hi, nc,
     $     cdir
     $     )
c
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cdir
      integer DIMDEC(f)
      REAL_T f(DIMV(f),nc)
      integer DIMDEC(c)
      REAL_T c(DIMV(c),nc)
c
      integer n
      integer i
      integer j
      REAL_T denom
      parameter(denom=0.50d0)
c
      if (cdir .eq. 0 ) then
         do n = 1, nc
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  c(i,j,n) = (f(2*i,2*j,n) + f(2*i,2*j+1,n))*denom
               enddo
            enddo
         enddo
      else if (cdir .eq. 1) then
         do n = 1, nc
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  c(i,j,n) = (f(2*i,2*j,n) + f(2*i+1,2*j,n))*denom
               enddo
            enddo
         enddo
      endif
c
      end
c-----------------------------------------------------------------------
      subroutine FORT_AVERAGECC (
     $     c, DIMS(c),
     $     f, DIMS(f),
     $     lo, hi, nc
     $     )
c
      integer nc
      integer DIMDEC(f)
      integer DIMDEC(c)
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)
c
      integer i
      integer j
      integer n
      REAL_T denom
      parameter(denom=0.25d0)
c
      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               c(i,j,n) =  (
     $              f(2*i+1,2*j+1,n) + f(2*i  ,2*j+1,n)
     $              + f(2*i+1,2*j  ,n) + f(2*i  ,2*j  ,n))*denom
            enddo
         enddo
      enddo
c
      end
c-----------------------------------------------------------------------
      subroutine FORT_APPLYBC (
     $     flagden, flagbc, maxorder,
     $     u,   DIMS(u),
     $     cdir, bct, bcl,
     $     bcval, DIMS(bcval),
     $     mask,  DIMS(mask),
     $     maskp, DIMS(maskp),
     $     maskm, DIMS(maskm),
     $     den,   DIMS(den),
     $     exttd, DIMS(exttd),
     $     tander,DIMS(tander),
     $     lo, hi, nc,
     $     h
     $     )
c
c     If the boundary is of Neumann type, set the ghost cell value to
c     that of the outermost point in the valid data (2nd order accurate)
c     and then fill the "den" array with the value "1"
c     
c     
c     If flagbc==1:
c     
c     If the boundary is of Dirichlet type, construct a polynomial
c     interpolation through the boundary location and internal points
c     (at locations x(-1:len-2) that generates the ghost cell value (at
c     location xInt).  Then fill the ghost cell with the interpolated value.
c     If flagden==1, load the "den" array with the interpolation
c     coefficient corresponding to outermost point in the valid region
c     ( the coef(0) corresponding to the location x(0) )
c      
c     Note: 
c     The bc type = LO_REFLECT_ODD is a special type of boundary condition.

c ::: the old FORT_APPLYBC put values out in ghost cells.  The new
c ::: FORT_APPLYBC makes a distinction between normal and tangential derivs.
c ::: the normal derivatives are still evaluated with ghost cells ( and
c ::: the den array for adjusting the diagonal element).  Tangential
c ::: derivatives are NOT computed from the ghost cells (note: we are
c ::: discussing here the tangential derivates which are centered on the
c ::: external cells walls of the rectangular domain.  Internal tangential
c ::: derivatives whose stencil extends outside the rectangular domain have
c ::: still ANOTHER calculational trick, one sided derivatives, which is 
c ::: implemented in the apply operator).  For these tangential derivatives,
c ::: if some of the cells in the stencil are uncovered, the derivative is
c ::: computed half from cells which are internal, plus an externally supplied
c ::: tangential derivative centered at location bcl.
c ??? exttd is the externally supplied tangential derivative which is
c ??? centered at bcl from wall.  Since the location of exttd is
c ??? indeterminate (determined by bcl), we arbitrarily index it in the
c ??? cell just outside the wall.
c ??? tander is the calculated tangential derivative which is centered at
c ??? at the wall.  It is indexed as if it were a cell centered quantity
c ??? just outside the wall.  This is anomalous since it is really edge
c ??? centered on the outer wall.
c      
c ::: other notes since previous developers didn't bother to document
c ::: cdir is mnemonic for coordinate direction, i.e. which side this is
c ::: cdir==0 -> left wall
c ::: cdir==2 -> right
c ::: cdir==1 -> bottom
c ::: cdir==3 -> top

c ::: notes about mask
c ::: since we now require other masks, there are 3.
c ::: 1) mask -- the main one covering the cells off the side of valid domain.
c ::: 2) maskp -- the plus end of transverse direction
c ::: 3) maskn -- the minus end of transverse direction

      integer maxorder
      integer nc, cdir, flagden, flagbc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(u)
      REAL_T u(DIMV(u),nc)
      integer DIMDEC(den)
      REAL_T den(DIMV(den),nc)
      integer DIMDEC(exttd)
      REAL_T exttd(DIMV(exttd),nc)
      integer DIMDEC(tander)
      REAL_T tander(DIMV(tander),nc)
      integer DIMDEC(bcval)
      REAL_T bcval(DIMV(bcval),nc)
      integer DIMDEC(mask)
      integer mask(DIMV(mask))
      integer DIMDEC(maskp)
      integer maskp(DIMV(maskp))
      integer DIMDEC(maskm)
      integer maskm(DIMV(maskm))
      integer bct(nc)
      REAL_T bcl
      REAL_T h(BL_SPACEDIM)
c
      REAL_T hx,hy,lambda
      REAL_T innder,outder
      REAL_T innloc,outloc,edgloc
      integer i
      integer j
      integer n
      logical is_dirichlet
      logical is_neumann
c
      integer lenx
      integer leny
      integer m
c
      integer Lmaxorder
      integer maxmaxorder
      parameter(maxmaxorder=4)
      REAL_T x(-1:maxmaxorder-2)
      REAL_T coef(-1:maxmaxorder-2)
      REAL_T xInt
c
      logical False, True
      parameter( False=.false.)
      parameter( True=.true. )
c
      is_dirichlet(i) = ( i .eq. LO_DIRICHLET   )
      is_neumann(i)   = ( i .eq. LO_NEUMANN )
      hx = h(1)
      hy = h(2)
c
      if ( maxorder .eq. -1 ) then
         Lmaxorder = maxmaxorder
      else
         Lmaxorder = MIN(maxorder,maxmaxorder)
      endif
      lenx = MIN(hi(1)-lo(1), Lmaxorder-2)
      leny = MIN(hi(2)-lo(2), Lmaxorder-2)
c
c     TODO:
c     In order for this to work with growing multigrid, must
c     sort xa[] because it is possible for the xb value to lay
c     within this range.
c     
c     The Left face of the grid
c
      if(cdir .eq. 0) then
         do n = 1, nc
            if (is_neumann(bct(n))) then
               do j = lo(2), hi(2)
                  u(lo(1)-1,j,n) = merge(
     $                 u(lo(1),j,n),
     $                 u(lo(1)-1,j,n),
     $                 mask(lo(1)-1,j) .gt. 0)
               enddo
               if ( flagden .eq. 1) then
                  do j = lo(2), hi(2)
                     den(lo(1),j,n) = 1.0
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: Null
c     :::::: interior part of left side
               i = lo(1)
               do j=lo(2)+1,hi(2)-1
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                  tander(i-1,j,n) = innder
               enddo
c     :::::: now end points
c     :::::: lower leftside
               j = lo(2)
               if( maskm(i,-1 + j).gt.0 ) then
                    innder = (-3*U(i,j,n) + 4*U(i,1 + j,n) - U(i,2 + j,n))/(2.*hy)
               else
                   innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
               endif
               tander(i-1,j,n) = innder
c     :::::: upper left side
                j = hi(2)
                if( maskp(i,1 + j).gt.0 ) then
                   innder = (U(i,-2 + j,n) - 4*U(i,-1 + j,n) + 3*U(i,j,n))/(2.*hy)
                else
                   innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                endif
                tander(i-1,j,n) = innder
            else if (is_dirichlet(bct(n))) then
               do m=0,lenx
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(1)
               xInt = - 0.5
               call amrex_tmg_polyInterpCoeff(xInt, x, lenx+2, coef)
               if ( flagbc .eq. 1 ) then
                  do j = lo(2), hi(2)
                     u(lo(1)-1, j, n) = merge(
     $                    bcval(lo(1)-1,j,n)*coef(-1),
     $                    u(lo(1)-1, j, n),
     $                    mask(lo(1)-1,j) .gt. 0)
                  enddo
               else
                  do j = lo(2), hi(2)
                     u(lo(1)-1, j, n) = merge(
     $                    0.0d0,
     $                    u(lo(1)-1, j, n),
     $                    mask(lo(1)-1,j) .gt. 0)
                  enddo
               endif
               do m = 0, lenx
                  do j = lo(2), hi(2)
                     u(lo(1)-1,j,n) = merge(
     $                    u(lo(1)-1,j,n)
     $                    + u(lo(1)+m, j, n)*coef(m),
     $                    u(lo(1)-1,j,n),
     $                    mask(lo(1)-1,j) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     den(lo(1),j,n) = merge(coef(0), 0.0d0,
     $                    mask(lo(1)-1,j) .gt. 0)
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: Null
c     :::::: interior part of left side
               edgloc = 0.
               i = lo(1)
               do j=lo(2)+1,hi(2)-1
                  if( mask(-1 + i,-1 + j).eq.0.and.mask(-1 + i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-U(-1 + i,-1 + j,n) + U(-1 + i,1 + j,n))/(2.*hy)
                  elseif( mask(-1 + i,j).eq.0.and.mask(-1 + i,1 + j).eq.0.and.mask(-1 + i,2 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-3*U(-1 + i,j,n) + 4*U(-1 + i,1 + j,n) - U(-1 + i,2 + j,n))/(2.*hy)
                  elseif( mask(-1 + i,-2 + j).eq.0.and.mask(-1 + i,-1 + j).eq.0.and.mask(-1 + i,j).eq.0 ) then
                     outloc = -0.5
                     outder = (U(-1 + i,-2 + j,n) - 4*U(-1 + i,-1 + j,n) + 3*U(-1 + i,j,n))/(2.*hy)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i-1,j,n)
                     else
                        outder = 0.
                     endif
                  endif
                  innloc = 0.5
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
            
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  tander(i-1,j,n) = lambda*innder+(1-lambda)*outder
               enddo
c     :::::: now end points
c     :::::: lower leftside
               j = lo(2)
               if( mask(-1 + i,-1 + j).eq.0.and.mask(-1 + i,1 + j).eq.0 ) then
                   outloc = -0.5
                   outder = (-U(-1 + i,-1 + j,n) + U(-1 + i,1 + j,n))/(2.*hy)
               elseif( mask(-1 + i,j).eq.0.and.mask(-1 + i,1 + j).eq.0.and.mask(-1 + i,2 + j).eq.0 ) then
                   outloc = -0.5
                   outder = (-3*U(-1 + i,j,n) + 4*U(-1 + i,1 + j,n) - U(-1 + i,2 + j,n))/(2.*hy)
               else
                   outloc = x(-1)
                   if( flagbc .eq. 1) then
                      outder = exttd(i-1,j,n)
                   else
                      outder = 0.
                   endif
               endif
               if( maskm(i,-1 + j).gt.0 ) then
                    innder = (-3*U(i,j,n) + 4*U(i,1 + j,n) - U(i,2 + j,n))/(2.*hy)
               else
                   innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i-1,j,n) = lambda*innder+(1-lambda)*outder
c     :::::: upper left side
                j = hi(2)
                if( mask(-1 + i,-1 + j).eq.0.and.mask(-1 + i,1 + j).eq.0 ) then
                   outloc = -0.5
                   outder = (-U(-1 + i,-1 + j,n) + U(-1 + i,1 + j,n))/(2.*hy)
                elseif( mask(-1 + i,-2 + j).eq.0.and.mask(-1 + i,-1 + j).eq.0.and.mask(-1 + i,j).eq.0 ) then
                   outloc = -0.5
                   outder = (U(-1 + i,-2 + j,n) - 4*U(-1 + i,-1 + j,n) + 3*U(-1 + i,j,n))/(2.*hy)
                else
                   outloc = x(-1)
                   if( flagbc .eq. 1) then
                      outder = exttd(i-1,j,n)
                   else
                      outder = 0.
                   endif
                endif
                if( maskp(i,1 + j).gt.0 ) then
                   innder = (U(i,-2 + j,n) - 4*U(i,-1 + j,n) + 3*U(i,j,n))/(2.*hy)
                else
                   innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                endif
                innloc = 0.5
                lambda = (edgloc-outloc)/(innloc-outloc)
                tander(i-1,j,n) = lambda*innder+(1-lambda)*outder
               
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               
               do j = lo(2), hi(2)
                  u(lo(1)-1, j, n) = merge(
     $                 -u(lo(1),j,n),
     $                 u(lo(1)-1, j, n),
     $                 mask(lo(1)-1,j) .gt. 0)
                  tander(lo(1)-1,j,n) = 0.0
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     den(lo(1),j,n) = merge(-1.0d0, 0.0d0,
     $                    mask(lo(1)-1,j) .gt. 0)
                  enddo
               endif
               
            else 
               print *,'UNKNOWN BC ON LEFT FACE IN APPLYBC'
               stop
            endif
         enddo
      endif
c     
c     The Right face of the grid
c
      if(cdir .eq. 2) then
         do n = 1, nc
            if(is_neumann(bct(n))) then
               do j = lo(2), hi(2)
                  u(hi(1)+1,j,n) = merge(
     $                 u(hi(1), j, n),
     $                 u(hi(1)+1, j, n),
     $                 mask(hi(1)+1,j) .gt. 0)
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     den(hi(1),j,n) = 1.0
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: Null
c     :::::: interior part of right side
               i = hi(1)
               do j=lo(2)+1,hi(2)-1
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                  tander(i+1,j,n) = innder
               enddo
c :::::: now end points
c :::::: lower right side
               j = lo(2)
               if( maskm(i,-1 + j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(i,1 + j,n) - U(i,2 + j,n))/(2.*hy)
               else
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
               endif
               tander(i+1,j,n) = innder
c     :::::: upper right side
               j = hi(2)
               if( maskp(i,1 + j).gt.0 ) then
                  innder = (U(i,-2 + j,n) - 4*U(i,-1 + j,n) + 3*U(i,j,n))/(2.*hy)
               else
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
               endif
               tander(i+1,j,n) = innder
            else if (is_dirichlet(bct(n))) then
               do m=0,lenx
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(1)
               xInt = - 0.5
               call amrex_tmg_polyInterpCoeff(xInt, x, lenx+2, coef)
               if ( flagbc .eq. 1 ) then
                  do j = lo(2), hi(2)
                     u(hi(1)+1,j,n) = merge(
     $                    bcval(hi(1)+1,j,n)*coef(-1),
     $                    u(hi(1)+1,j,n),
     $                    mask(hi(1)+1,j) .gt. 0)
                  enddo
               else
                  do j = lo(2), hi(2)
                     u(hi(1)+1,j,n) = merge(
     $                    0.0d0,
     $                    u(hi(1)+1,j,n),
     $                    mask(hi(1)+1,j) .gt. 0)
                  enddo
               endif
               do m = 0, lenx
                  do j = lo(2), hi(2)
                     u(hi(1)+1,j,n) = merge(
     $                    u(hi(1)+1,j,n)
     $                    + u(hi(1)-m,j,n)*coef(m),
     $                    u(hi(1)+1,j,n),
     $                    mask(hi(1)+1,j) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     den(hi(1),j,n)   = merge(coef(0), 0.0d0,
     $                    mask(hi(1)+1,j) .gt. 0)
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: Null
c     :::::: interior part of right side
               edgloc = 0.
               i = hi(1)
               do j=lo(2)+1,hi(2)-1
                  if( mask(1 + i,-1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-U(1 + i,-1 + j,n) + U(1 + i,1 + j,n))/(2.*hy)
                  elseif( mask(1 + i,j).eq.0.and.mask(1 + i,1 + j).eq.0.and.mask(1 + i,2 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-3*U(1 + i,j,n) + 4*U(1 + i,1 + j,n) - U(1 + i,2 + j,n))/(2.*hy)
                  elseif( mask(1 + i,-2 + j).eq.0.and.mask(1 + i,-1 + j).eq.0.and.mask(1 + i,j).eq.0 ) then
                     outloc = -0.5
                     outder = (U(1 + i,-2 + j,n) - 4*U(1 + i,-1 + j,n) + 3*U(1 + i,j,n))/(2.*hy)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i+1,j,n)
                     else
                        outder = 0.
                     endif
                  endif
                  innloc = 0.5
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)

                  lambda = (edgloc-outloc)/(innloc-outloc)
                  tander(i+1,j,n) = lambda*innder+(1-lambda)*outder
               enddo
c :::::: now end points
c :::::: lower right side
               j = lo(2)
               if( mask(1 + i,-1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(1 + i,-1 + j,n) + U(1 + i,1 + j,n))/(2.*hy)
               elseif( mask(1 + i,j).eq.0.and.mask(1 + i,1 + j).eq.0.and.mask(1 + i,2 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-3*U(1 + i,j,n) + 4*U(1 + i,1 + j,n) - U(1 + i,2 + j,n))/(2.*hy)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i+1,j,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskm(i,-1 + j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(i,1 + j,n) - U(i,2 + j,n))/(2.*hy)
               else
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i+1,j,n) = lambda*innder+(1-lambda)*outder
c     :::::: upper right side
               j = hi(2)
               if( mask(1 + i,-1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(1 + i,-1 + j,n) + U(1 + i,1 + j,n))/(2.*hy)
               elseif( mask(1 + i,-2 + j).eq.0.and.mask(1 + i,-1 + j).eq.0.and.mask(1 + i,j).eq.0 ) then
                  outloc = -0.5
                  outder = (U(1 + i,-2 + j,n) - 4*U(1 + i,-1 + j,n) + 3*U(1 + i,j,n))/(2.*hy)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i+1,j,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskp(i,1 + j).gt.0 ) then
                  innder = (U(i,-2 + j,n) - 4*U(i,-1 + j,n) + 3*U(i,j,n))/(2.*hy)
               else
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i+1,j,n) = lambda*innder+(1-lambda)*outder
               
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               
               do j = lo(2), hi(2)
                  u(hi(1)+1, j, n) = merge(
     $                 -u(hi(1),j,n),
     $                 u(hi(1)+1, j, n),
     $                 mask(hi(1)+1,j) .gt. 0)
c                 ::: tangential der is trivial
                  tander(hi(1)+1,j,n) = 0.0
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     den(hi(1),j,n) = merge(-1.0d0, 0.0d0,
     $                    mask(hi(1)+1,j) .gt. 0)
                  enddo
               endif
               
            else
               print *,'UNKNOWN BC ON RIGHT FACE IN APPLYBC'
               stop
            endif
         enddo
      endif
c
c     The Bottom of the Grid
c
      if(cdir .eq. 1) then
         do n = 1, nc
            if(is_neumann(bct(n))) then
               do i = lo(1),hi(1)
                  u(i,lo(2)-1,n) = merge(
     $                 u(i,lo(2),n),
     $                 u(i,lo(2)-1,n),
     $                 mask(i,lo(2)-1) .gt. 0)
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1),hi(1)
                     den(i,lo(2),n)   = 1.0
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: interior part of bottom side
               j = lo(2)
               do i=lo(1)+1,hi(1)-1
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
                  tander(i,j-1,n) = innder
               enddo
c :::::: now end points
c :::::: left bottom side
               i = lo(1)
               if( maskm(-1 + i,j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(1 + i,j,n) - U(2 + i,j,n))/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
               endif
               tander(i,j-1,n) = innder
c :::::: right bottom side
               i = hi(1)
               if( maskp(1 + i,j).gt.0 ) then
                  innder = (U(-2 + i,j,n) - 4*U(-1 + i,j,n) + 3*U(i,j,n))/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
               endif
               tander(i,j-1,n) = innder
            else if (is_dirichlet(bct(n))) then
               do m=0,leny
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(2)
               xInt = - 0.5
               call amrex_tmg_polyInterpCoeff(xInt, x, leny+2, coef)
               if ( flagbc .eq. 1 ) then
                  do i = lo(1), hi(1)
                     u(i,lo(2)-1,n) = merge(
     $                    bcval(i,lo(2)-1,n)*coef(-1),
     $                    u(i,lo(2)-1,n),
     $                    mask(i,lo(2)-1) .gt. 0)
                  enddo
               else
                  do i = lo(1), hi(1)
                     u(i,lo(2)-1,n) = merge(
     $                    0.0d0,
     $                    u(i,lo(2)-1,n),
     $                    mask(i,lo(2)-1) .gt. 0)
                  enddo
               endif
               do m = 0, leny
                  do i = lo(1), hi(1)
                     u(i, lo(2)-1, n) = merge(
     $                    u(i, lo(2)-1,n)
     $                    + u(i, lo(2)+m,n)*coef(m),
     $                    u(i, lo(2)-1, n),
     $                    mask(i, lo(2)-1) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1), hi(1)
                     den(i, lo(2),n)   = merge(coef(0), 0.0d0,
     $                    mask(i, lo(2)-1) .gt. 0)
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: interior part of bottom side
               edgloc = 0.
               j = lo(2)
               do i=lo(1)+1,hi(1)-1
                  if(  mask(-1 + i,-1 + j).eq.0.and.mask(1 + i,-1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-U(-1 + i,-1 + j,n) + U(1 + i,-1 + j,n))/(2.*hx)
                  else if(  mask(i,-1 + j).eq.0.and.mask(1 + i,-1 + j).eq.0.and.mask(2 + i,-1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-3*U(i,-1 + j,n) + 4*U(1 + i,-1 + j,n) - U(2 + i,-1 + j,n))/(2.*hx)
                  else if(  mask(-2 + i,-1 + j).eq.0.and.mask(-1 + i,-1 + j).eq.0.and.mask(i,-1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (U(-2 + i,-1 + j,n) - 4*U(-1 + i,-1 + j,n) + 3*U(i,-1 + j,n))/(2.*hx)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i,j-1,n)
                     else
                        outder = 0.
                     endif
                  endif
                  innloc = 0.5
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)

                  lambda = (edgloc-outloc)/(innloc-outloc)
                  tander(i,j-1,n) = lambda*innder+(1-lambda)*outder
               enddo
c :::::: now end points
c :::::: left bottom side
               i = lo(1)
               if(  mask(-1 + i,-1 + j).eq.0.and.mask(1 + i,-1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(-1 + i,-1 + j,n) + U(1 + i,-1 + j,n))/(2.*hx)
               else if( mask(i,-1 + j).eq.0.and.mask(1 + i,-1 + j).eq.0.and.mask(2 + i,-1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-3*U(i,-1 + j,n) + 4*U(1 + i,-1 + j,n) - U(2 + i,-1 + j,n))/(2.*hx)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i,j-1,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskm(-1 + i,j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(1 + i,j,n) - U(2 + i,j,n))/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i,j-1,n) = lambda*innder+(1-lambda)*outder
c :::::: right bottom side
               i = hi(1)
               if( mask(-1 + i,-1 + j).eq.0.and.mask(1 + i,-1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(-1 + i,-1 + j,n) + U(1 + i,-1 + j,n))/(2.*hx)
               else if( mask(-2 + i,-1 + j).eq.0.and.mask(-1 + i,-1 + j).eq.0.and.mask(i,-1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (U(-2 + i,-1 + j,n) - 4*U(-1 + i,-1 + j,n) + 3*U(i,-1 + j,n))/(2.*hx)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i,j-1,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskp(1 + i,j).gt.0 ) then
                  innder = (U(-2 + i,j,n) - 4*U(-1 + i,j,n) + 3*U(i,j,n))/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i,j-1,n) = lambda*innder+(1-lambda)*outder
               
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               
               do i = lo(1), hi(1)
                  u(i,lo(2)-1,n) = merge(
     $                 -u(i,lo(2),n),
     $                 u(i,lo(2)-1,n),
     $                 mask(i,lo(2)-1) .gt. 0)
                  tander(i,lo(2)-1,n) = 0.0
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1), hi(1)
                     den(i,lo(2),n) = merge(-1.0d0, 0.0d0,
     $                    mask(i,lo(2)-1) .gt. 0)
                  enddo
               endif
               
            else
               print *,'UNKNOWN BC ON BOTTOM FACE IN APPLYBC'
               stop
            endif
        enddo
      endif
c     
c     The top of the grid
c
      if (cdir .eq. 3) then
         do n = 1, nc
            if(is_neumann(bct(n))) then
               do i = lo(1), hi(1)
                  u(i,hi(2)+1,n) = merge(
     $                 u(i,hi(2),n),
     $                 u(i,hi(2)+1,n),
     $                 mask(i,hi(2)+1) .gt. 0)
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1), hi(1)
                     den(i,hi(2),n)   = 1.0
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: interior part of top side
               j = hi(2)
               do i=lo(1)+1,hi(1)-1
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
                  tander(i,j+1,n) = innder
               enddo
c     :::::: now end points
c     :::::: left top side
               i = lo(1)
               if( maskm(-1 + i,j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(1 + i,j,n) - U(2 + i,j,n))/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
               endif
               tander(i,j+1,n) = innder
c     :::::: right top side
               i = hi(1)
               if( maskp(1 + i,j).gt.0 ) then
                  innder = (U(-2 + i,j,n) - 4*U(-1 + i,j,n) + 3*U(i,j,n))/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
               endif
               tander(i,j+1,n) = innder
            else if (is_dirichlet(bct(n))) then
               if ( bct(n) .eq. LO_REFLECT_ODD ) leny = 0
               do m=0,leny
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(2)
               xInt = - 0.5
               call amrex_tmg_polyInterpCoeff(xInt, x, leny+2, coef)
               if ( flagbc .eq. 1 ) then
                  do i = lo(1), hi(1)
                     u(i,hi(2)+1,n) = merge(
     $                    bcval(i,hi(2)+1,n)*coef(-1),
     $                    u(i,hi(2)+1,n),
     $                    mask(i,hi(2)+1) .gt. 0)
                  enddo
               else
                  do i = lo(1), hi(1)
                     u(i,hi(2)+1,n) = merge(
     $                    0.0d0,
     $                    u(i,hi(2)+1,n),
     $                    mask(i,hi(2)+1) .gt. 0)
                  enddo
               endif
               do m = 0, leny
                  do i = lo(1), hi(1)
                     u(i, hi(2)+1,n) = merge(
     $                    u(i,hi(2)+1,n)
     $                    + u(i, hi(2)-m,n)*coef(m),
     $                    u(i,hi(2)+1,n),
     $                    mask(i,hi(2)+1) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1), hi(1)
                     den(i,hi(2),n)   = merge(coef(0), 0.0d0,
     $                    mask(i,hi(2)+1) .gt. 0)
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: interior part of top side
               edgloc = 0.
               j = hi(2)
               do i=lo(1)+1,hi(1)-1
                  if(  mask(-1 + i,1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-U(-1 + i,1 + j,n) + U(1 + i,1 + j,n))/(2.*hx)
                  else if(  mask(i,1 + j).eq.0.and.mask(1 + i,1 + j).eq.0.and.mask(2 + i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-3*U(i,1 + j,n) + 4*U(1 + i,1 + j,n) - U(2 + i,1 + j,n))/(2.*hx)
                  else if( mask(-2 + i,1 + j).eq.0.and.mask(-1 + i,1 + j).eq.0.and.mask(i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (U(-2 + i,1 + j,n) - 4*U(-1 + i,1 + j,n) + 3*U(i,1 + j,n))/(2.*hx)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i,j+1,n)
                     else
                        outder = 0.
                     endif
                  endif
                  innloc = 0.5
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
                  
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  tander(i,j+1,n) = lambda*innder+(1-lambda)*outder
               enddo
c     :::::: now end points
c     :::::: left top side
               i = lo(1)
               if( mask(-1 + i,1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(-1 + i,1 + j,n) + U(1 + i,1 + j,n))/(2.*hx)
               else if( mask(i,1 + j).eq.0.and.mask(1 + i,1 + j).eq.0.and.mask(2 + i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-3*U(i,1 + j,n) + 4*U(1 + i,1 + j,n) - U(2 + i,1 + j,n))/(2.*hx)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i,j+1,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskm(-1 + i,j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(1 + i,j,n) - U(2 + i,j,n))/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i,j+1,n) = lambda*innder+(1-lambda)*outder
c     :::::: right top side
               i = hi(1)
               if( mask(-1 + i,1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(-1 + i,1 + j,n) + U(1 + i,1 + j,n))/(2.*hx)
               else if(  mask(-2 + i,1 + j).eq.0.and.mask(-1 + i,1 + j).eq.0.and.mask(i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (U(-2 + i,1 + j,n) - 4*U(-1 + i,1 + j,n) + 3*U(i,1 + j,n))/(2.*hx)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i,j+1,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskp(1 + i,j).gt.0 ) then
                  innder = (U(-2 + i,j,n) - 4*U(-1 + i,j,n) + 3*U(i,j,n))/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i,j+1,n) = lambda*innder+(1-lambda)*outder
               
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               
               do i = lo(1), hi(1)
                  u(i,hi(2)+1,n) = merge(
     $                 -u(i,hi(2),n),
     $                 u(i,hi(2)+1,n),
     $                 mask(i,hi(2)+1) .gt. 0)
                  tander(i,hi(2)+1,n) = 0.0
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1), hi(1)
                     den(i,hi(2),n) = merge(-1.0d0, 0.0d0,
     $                    mask(i,hi(2)+1) .gt. 0)
                  enddo
               endif
               
            else
               print *,'UNKNOWN BC ON TOP FACE IN APPLYBC'
               stop
            endif
         enddo
      endif
c
      end
