#include "BC_TYPES.H"
#include "REAL.H"
#include "CONSTANTS.H"
#include "MCINTERPBNDRYDATA_F.H"
#include "ArrayLim.H"

#define SDIM 2
#define NUMDERIV 2
#define XDER   1
#define X2DER  2
      
      
c ---------------------------------------------------------------
c ::  FORT_BDINTERPXLO : Interpolation on Xlo Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  DIMS(bdry)  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  DIMS(cb)    => index limits of coarsened grid interior
c ::  nvar        => number of variables in solution, half number in bdry
c ::  ratio       => refinement ratio
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPXLO (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratio,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives,hfine)

      integer  nvar, ratio, not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(mask)
      integer  DIMDEC(crse)
      integer  DIMDEC(cb)
      REAL_T   bdry(DIMV(bdry),2*nvar)
      REAL_T   derives(DIM2(cb),NUMDERIV)      
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T   hfine
      REAL_T  hcrse

      REAL_T   xx
      integer  i, j, ic, jc, off, n
      integer  jclo, jchi

      jclo = ARG_L2(cb)
      jchi = ARG_H2(cb)
      ic   = ARG_L1(cb)-1
      i    = lo(1)-1
      hcrse = ratio*hfine

      do n = 1, nvar
c        ::::: define interp coefs
         do jc = jclo, jchi
            j = ratio*jc
            if( mask(i,j-1).eq.not_covered .and.
     &          mask(i,j+ratio).eq.not_covered ) then
               derives(jc,XDER)  = (crse(ic,jc+1,n)-crse(ic,jc-1,n))/(2*hcrse)
               derives(jc,X2DER) = (crse(ic,jc+1,n)
     &           -                  crse(ic,jc  ,n) * two
     &           +                  crse(ic,jc-1,n))/hcrse**2
            else if( mask(i,j-1).ne.not_covered .and.
     &               mask(i,j+ratio).ne.not_covered ) then
c              ::: drop order by 2
               derives(jc,XDER )  = zero
               derives(jc,X2DER)  = zero
            else if( mask(i,j+ratio).eq. not_covered ) then
c              ::: j-1 bad, so try + dir
               if( jc+2 .le. jchi+1 ) then
                  if( mask(i,j+2*ratio).eq.not_covered ) then
c                    ::: can use 3 point stencils for j,j+1,j+2
                     derives(jc,XDER) = (-crse(ic,jc+2,n)+4*crse(ic,jc+1,n)-
     &                    3*crse(ic,jc,n))/(2*hcrse)
                     derives(jc,X2DER) = (crse(ic,jc+2,n)-2*crse(ic,jc+1,n)
     &                    +crse(ic,jc,n))/hcrse**2
                  else
c                    ::: use j,j+1, drop order
                     derives(jc,XDER) = (crse(ic,jc+1,n)-crse(ic,jc,n))/hcrse
                     derives(jc,X2DER) = zero
                  endif
               else
c                 ::: use j,j+1, drop order
                  derives(jc,XDER) = (crse(ic,jc+1,n)-crse(ic,jc,n))/hcrse
                  derives(jc,X2DER) = zero
               endif
            else if( mask(i,j-1) .eq. not_covered ) then
c              ::: try - direction
               if( jc-2 .ge. jclo-1 ) then
                  if( mask(i,j-1) .eq. not_covered ) then
c                    ::: can use 3 point stencils
                     derives(jc,XDER) = (crse(ic,jc-2,n)-4*crse(ic,jc-1,n)
     &                    +3*crse(ic,jc,n))/(2*hcrse)
                     derives(jc,X2DER) = (crse(ic,jc-2,n)-2*crse(ic,jc-1,n)
     &                    +crse(ic,jc,n))/hcrse**2
                  else
                     derives(jc,XDER) = (crse(ic,jc,n)-crse(ic,jc-1,n))/hcrse
                     derives(jc,X2DER) = zero
                  endif
               else
                  derives(jc,XDER) = (crse(ic,jc,n)-crse(ic,jc-1,n))/hcrse
                  derives(jc,X2DER) = zero
               endif
            else
               write(6,*)'impossible BDINTERPYLO'
               stop
            endif
         enddo
         do off = 0, ratio - 1
            xx = (off+0.5)*hfine-hcrse/2
            do jc = jclo, jchi
               j = ratio*jc + off
               bdry(i,j,n) = crse(ic,jc,n)
     &              +derives(jc,XDER)*xx
     &              +derives(jc,X2DER)*xx**2/2
               bdry(i,j,n+nvar) =
     &              derives(jc,XDER)+derives(jc,X2DER)*xx
            enddo
         enddo
      enddo
               

      return
      end

c ---------------------------------------------------------------
c ::  FORT_BDINTERPXHI : Interpolation on Xhi Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  DIMS(bdry)  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  DIMS(cb)    => index limits of coarsened grid interior
c ::  nvar        => number of variables in solution, half number in bdry
c ::  ratio       => refinement ratio
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPXHI (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratio,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives,hfine)

      integer  nvar, ratio, not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(mask)
      integer  DIMDEC(cb)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),2*nvar)
      REAL_T   derives(DIM2(cb),NUMDERIV)      
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T  hfine
      REAL_T  hcrse

      REAL_T   xx
      integer  i, j, ic, jc, off, n
      integer  jclo, jchi

      jclo = ARG_L2(cb)
      jchi = ARG_H2(cb)
      ic   = ARG_H1(cb)+1
      i    = hi(1)+1
      hcrse = ratio*hfine
      
      do n = 1, nvar
c        ::::: define interp coefs
         do jc = jclo, jchi
            j = ratio*jc
            if( mask(i,j-1).eq.not_covered .and.
     &          mask(i,j+ratio).eq.not_covered ) then
               derives(jc,XDER)  = (crse(ic,jc+1,n)-crse(ic,jc-1,n))/(2*hcrse)
               derives(jc,X2DER) = (crse(ic,jc+1,n)
     &           -                  crse(ic,jc  ,n) * two
     &           +                  crse(ic,jc-1,n))/hcrse**2
            else if( mask(i,j-1).ne.not_covered .and.
     &               mask(i,j+ratio).ne.not_covered ) then
c              ::: drop order by 2
               derives(jc,XDER )  = zero
               derives(jc,X2DER)  = zero
            else if( mask(i,j+ratio).eq. not_covered ) then
c              ::: j-1 bad, so try + dir
               if( jc+2 .le. jchi+1 ) then
                  if( mask(i,j+2*ratio).eq.not_covered ) then
c                    ::: can use 3 point stencils for j,j+1,j+2
                     derives(jc,XDER) = (-crse(ic,jc+2,n)+4*crse(ic,jc+1,n)-
     &                    3*crse(ic,jc,n))/(2*hcrse)
                     derives(jc,X2DER) = (crse(ic,jc+2,n)-2*crse(ic,jc+1,n)
     &                    +crse(ic,jc,n))/hcrse**2
                  else
c                    ::: use j,j+1, drop order
                     derives(jc,XDER) = (crse(ic,jc+1,n)-crse(ic,jc,n))/hcrse
                     derives(jc,X2DER) = zero
                  endif
               else
c                 ::: use j,j+1, drop order
                  derives(jc,XDER) = (crse(ic,jc+1,n)-crse(ic,jc,n))/hcrse
                  derives(jc,X2DER) = zero
               endif
            else if( mask(i,j-1) .eq. not_covered ) then
c              ::: try - direction
               if( jc-2 .ge. jclo-1 ) then
                  if( mask(i,j-1) .eq. not_covered ) then
c                    ::: can use 3 point stencils
                     derives(jc,XDER) = (crse(ic,jc-2,n)-4*crse(ic,jc-1,n)
     &                    +3*crse(ic,jc,n))/(2*hcrse)
                     derives(jc,X2DER) = (crse(ic,jc-2,n)-2*crse(ic,jc-1,n)
     &                    +crse(ic,jc,n))/hcrse**2
                  else
                     derives(jc,XDER) = (crse(ic,jc,n)-crse(ic,jc-1,n))/hcrse
                     derives(jc,X2DER) = zero
                  endif
               else
                  derives(jc,XDER) = (crse(ic,jc,n)-crse(ic,jc-1,n))/hcrse
                  derives(jc,X2DER) = zero
               endif
            else
               write(6,*)'impossible BDINTERPYLO'
               stop
            endif
         enddo
         do off = 0, ratio - 1
            xx = (off+0.5)*hfine-hcrse/2
            do jc = jclo, jchi
               j = ratio*jc + off
               bdry(i,j,n) = crse(ic,jc,n)
     &              +derives(jc,XDER)*xx
     &              +derives(jc,X2DER)*xx**2/2
               bdry(i,j,n+nvar) =
     &              derives(jc,XDER)+derives(jc,X2DER)*xx
            enddo
         enddo
      enddo

      return
      end

c ---------------------------------------------------------------
c ::  FORT_BDINTERPYLO : Interpolation on Ylo Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  DIMS(bdry)  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  DIMS(cb)    => index limits of coarsened grid interior
c ::  nvar        => number of variables in solution, half number in bdry
c ::  ratio       => refinement ratio
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPYLO (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratio,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives,hfine)

      integer  nvar, ratio, not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(mask)
      integer  DIMDEC(cb)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),2*nvar)
      REAL_T   derives(DIM1(cb),NUMDERIV)
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T   xx
      integer  i, j, ic, jc, off, n
      integer  iclo, ichi
      REAL_T  hfine
      REAL_T  hcrse

      iclo = ARG_L1(cb)
      ichi = ARG_H1(cb)
      jc   = ARG_L2(cb)-1
      j    = lo(2)-1
      hcrse = ratio*hfine

      do n = 1, nvar
c        ::::: define interp coefs
         do ic = iclo, ichi
            i = ratio*ic
            if( mask(i-1,j).eq.not_covered .and. 
     &           mask(i+ratio,j).eq.not_covered ) then
c           :::: use adjoining points
               derives(ic,XDER ) = (crse(ic+1,jc,n)-crse(ic-1,jc,n))/(2*hcrse)
               derives(ic,X2DER) = (crse(ic+1,jc,n)-2*crse(ic,jc,n)+
     &              crse(ic-1,jc,n))/hcrse**2
            else if( mask(i-1,j).ne.not_covered .and.
     &               mask(i+ratio,j).ne.not_covered ) then
c              :::: must drop order by two
               derives(ic,XDER ) = zero
               derives(ic,X2DER) = zero
            else if( mask(i+ratio,j).eq.not_covered ) then
c              ::: i-1 bad, try other direction
               if( ic+2 .le. ichi+1 ) then
                  if( mask(i+2*ratio,j).eq.not_covered ) then
c                    ::: can use 3 points stencils for i,i+1,i+2
                     derives(ic,XDER ) = (-crse(ic+2,jc,n)+4*crse(ic+1,jc,n)-
     &                    3*crse(ic,jc,n))/(2*hcrse)
                     derives(ic,X2DER) = (crse(ic+2,jc,n)-2*crse(ic+1,jc,n)+
     &                    crse(ic,jc,n))/hcrse**2
                  else
c                    ::: use i,i+1
                     derives(ic,XDER ) = (crse(ic+1,jc,n)-crse(ic,jc,n))/hcrse
                     derives(ic,X2DER) = zero
                  endif
              else
c                ::: use i,i+1
                 derives(ic,XDER ) = (crse(ic+1,jc,n)-crse(ic,jc,n))/hcrse
                 derives(ic,X2DER) = zero
              endif
           else if( mask(i-1,j).eq.not_covered ) then
c             ::: i+1 bad, try other direction
              if( ic-2 .ge. iclo-1 ) then
                 if( mask(i-1,j) .eq. not_covered ) then
c                  ::: can use 3 points stencils for i,i-1,i-2
                    derives(ic,XDER ) = (+crse(ic-2,jc,n)-4*crse(ic-1,jc,n)+
     &                    3*crse(ic,jc,n))/(2*hcrse)
                    derives(ic,X2DER) = (crse(ic-2,jc,n)-2*crse(ic-1,jc,n)+
     &                    crse(ic,jc,n))/hcrse**2
                 else
c                   ::: use i,i-1
                     derives(ic,XDER ) = (crse(ic,jc,n)-crse(ic-1,jc,n))/hcrse
                     derives(ic,X2DER) = zero
                 endif
              else
c                ::: use i,i-1
                 derives(ic,XDER ) = (crse(ic,jc,n)-crse(ic-1,jc,n))/hcrse
                 derives(ic,X2DER) = zero
              endif
           else
               write(6,*)'impossible BDINTERPYLO'
              stop
           endif
         enddo
         do off = 0, ratio - 1
            xx = (off+0.5)*hfine-hcrse/2
            do ic = iclo, ichi
               i = ratio*ic + off
               bdry(i,j,n) = crse(ic,jc,n)
     &              + derives(ic,XDER)*xx
     &              + derives(ic,X2DER)*xx**2/2
               bdry(i,j,n+nvar) =
     &              + derives(ic,XDER) + derives(ic,X2DER)*xx
            enddo
         enddo
      enddo

      return
      end

c ---------------------------------------------------------------
c ::  FORT_BDINTERPYHI : Interpolation on Yhi Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  DIMS(bdry)  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  DIMS(cb)    => index limits of coarsened grid interior
c ::  nvar        => number of variables in solution, half number in bdry
c ::  ratio       => refinement ratio
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPYHI (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratio,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives,hfine)

      integer  nvar, ratio, not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(mask)
      integer  DIMDEC(cb)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),2*nvar)
      REAL_T   derives(DIM1(cb),NUMDERIV)
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T   xx
      integer  i, j, ic, jc, off, n
      integer  iclo, ichi
      REAL_T  hfine
      REAL_T  hcrse

      iclo = ARG_L1(cb)
      ichi = ARG_H1(cb)
      jc   = ARG_H2(cb)+1
      j    = hi(2)+1
      hcrse = ratio*hfine

      do n = 1, nvar
c        ::::: define interp coefs
         do ic = iclo, ichi
            i = ratio*ic
            if( mask(i-1,j).eq.not_covered .and. 
     &           mask(i+ratio,j).eq.not_covered ) then
c           :::: use adjoining points
               derives(ic,XDER ) = (crse(ic+1,jc,n)-crse(ic-1,jc,n))/(2*hcrse)
               derives(ic,X2DER) = (crse(ic+1,jc,n)-2*crse(ic,jc,n)+
     &              crse(ic-1,jc,n))/hcrse**2
            else if( mask(i-1,j).ne.not_covered .and.
     &               mask(i+ratio,j).ne.not_covered ) then
c              :::: must drop order by two
               derives(ic,XDER ) = zero
               derives(ic,X2DER) = zero
            else if( mask(i+ratio,j).eq.not_covered ) then
c              ::: i-1 bad, try other direction
               if( ic+2 .le. ichi+1 ) then
                  if( mask(i+2*ratio,j).eq.not_covered ) then
c                    ::: can use 3 points stencils for i,i+1,i+2
                     derives(ic,XDER ) = (-crse(ic+2,jc,n)+4*crse(ic+1,jc,n)-
     &                    3*crse(ic,jc,n))/(2*hcrse)
                     derives(ic,X2DER) = (crse(ic+2,jc,n)-2*crse(ic+1,jc,n)+
     &                    crse(ic,jc,n))/hcrse**2
                  else
c                    ::: use i,i+1
                     derives(ic,XDER ) = (crse(ic+1,jc,n)-crse(ic,jc,n))/hcrse
                     derives(ic,X2DER) = zero
                  endif
              else
c                ::: use i,i+1
                 derives(ic,XDER ) = (crse(ic+1,jc,n)-crse(ic,jc,n))/hcrse
                 derives(ic,X2DER) = zero
              endif
           else if( mask(i-1,j).eq.not_covered ) then
c             ::: i+1 bad, try other direction
              if( ic-2 .ge. iclo-1 ) then
                 if( mask(i-1,j) .eq. not_covered ) then
c                  ::: can use 3 points stencils for i,i-1,i-2
                    derives(ic,XDER ) = (+crse(ic-2,jc,n)-4*crse(ic-1,jc,n)+
     &                    3*crse(ic,jc,n))/(2*hcrse)
                    derives(ic,X2DER) = (crse(ic-2,jc,n)-2*crse(ic-1,jc,n)+
     &                    crse(ic,jc,n))/hcrse**2
                 else
c                   ::: use i,i-1
                     derives(ic,XDER ) = (crse(ic,jc,n)-crse(ic-1,jc,n))/hcrse
                     derives(ic,X2DER) = zero
                 endif
              else
c                ::: use i,i-1
                 derives(ic,XDER ) = (crse(ic,jc,n)-crse(ic-1,jc,n))/hcrse
                 derives(ic,X2DER) = zero
              endif
           else
               write(6,*)'impossible BDINTERPYHI'
              stop
           endif
         enddo
         do off = 0, ratio - 1
            xx = (off+0.5)*hfine-hcrse/2
            do ic = iclo, ichi
               i = ratio*ic + off
               bdry(i,j,n) = crse(ic,jc,n)
     &              + derives(ic,XDER)*xx
     &              + derives(ic,X2DER)*xx**2/2
               bdry(i,j,n+nvar) =
     &              + derives(ic,XDER) + derives(ic,X2DER)*xx
            enddo
         enddo
      enddo

      return
      end


c ------------------------------------------------------------------------
c :: FORT_BDIDERIVXLO : copy from Xlo face of fine fab (getting values
c ::      which were stuck there by the physical boundary condition routines)
c ::      to get the first nvar values in bdry.  The second nvar values 
c ::      are gotten by computing the transverse derivative.  Two formulas
c ::      are used: 1) on the interior a centered difference is used; 
c ::      2) at the end points a one-sided difference is used.  No check
c ::      is made to check that the one-sided difference does not poke out
c ::      the other side of the grid.  But if grids are that small, a lot
c ::      of things will be blowing up.

c :: Inputs/Outputs:
c ::  bdry <= fine grid bndry data strip
c ::  lo,hi => index limits of "grid". So all action takes place in cells
c ::           outside the grid
c ::  fine  => fine grid data which has phys-bc data stuffed in it
c ::  nvar  => number of components in fine.  But you write twice as many
c ::           components of data in bdry because you compute transverse
c ::           derivatives as well
c ::  hfine => delta x in transverse direction; needed for derivs
c ------------------------------------------------------------------------
      subroutine FORT_BDIDERIVXLO( bdry, DIMS(bdry),
     &     lo, hi,
     &     fine, DIMS(fine),
     &     nvar, hfine)
      integer nvar
      REAL_T  hfine
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(bdry)
      integer DIMDEC(fine)
      REAL_T  bdry(DIMV(bdry),2*nvar)
      REAL_T  fine(DIMV(fine),  nvar)

      integer i,j,n

      i = lo(1)-1
      do n=1,nvar
c ::: ::: copy boundary values over
         do j=lo(2),hi(2)
            bdry(i,j,n) = fine(i,j,n)
         enddo
c ::: ::: do stinking calculation of derivativ using centered difference
c ::: ::: this ONLY makes sense for DIRICHLET BC.  For other BC,
c ::: ::: we will not use these values.
         do j=lo(2)+1,hi(2)-1
            bdry(i,j,n+nvar) = (fine(i,j+1,n)-fine(i,j-1,n))/(2*hfine)
         enddo
         j = lo(2)
         bdry(i,j,n+nvar) = (-fine(i,j+2,n)+4*fine(i,j+1,n)-3*fine(i,j,n))/
     &        (2*hfine)
         j = hi(2)
         bdry(i,j,n+nvar) = (+fine(i,j-2,n)-4*fine(i,j-1,n)+3*fine(i,j,n))/
     &        (2*hfine)
      enddo

      return
      end
      

      subroutine FORT_BDIDERIVXHI( bdry, DIMS(bdry),
     &     lo, hi,
     &     fine, DIMS(fine),
     &     nvar, hfine)
      integer nvar
      REAL_T  hfine
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(bdry)
      integer DIMDEC(fine)
      REAL_T  bdry(DIMV(bdry),2*nvar)
      REAL_T  fine(DIMV(fine),  nvar)

      integer i,j,n

      i = hi(1)+1
      do n=1,nvar
c ::: ::: copy boundary values over
         do j=lo(2),hi(2)
            bdry(i,j,n) = fine(i,j,n)
         enddo
c ::: ::: do stinking calculation of derivativ using centered difference
c ::: ::: this ONLY makes sense for DIRICHLET BC.  For other BC,
c ::: ::: we will not use these values.
         do j=lo(2)+1,hi(2)-1
            bdry(i,j,n+nvar) = (fine(i,j+1,n)-fine(i,j-1,n))/(2*hfine)
         enddo
         j = lo(2)
         bdry(i,j,n+nvar) = (-fine(i,j+2,n)+4*fine(i,j+1,n)-3*fine(i,j,n))/
     &        (2*hfine)
         j = hi(2)
         bdry(i,j,n+nvar) = (+fine(i,j-2,n)-4*fine(i,j-1,n)+3*fine(i,j,n))/
     &        (2*hfine)
      enddo

      return
      end
      

      subroutine FORT_BDIDERIVYLO( bdry, DIMS(bdry),
     &     lo, hi,
     &     fine, DIMS(fine),
     &     nvar, hfine)
      integer nvar
      REAL_T  hfine
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(bdry)
      integer DIMDEC(fine)
      REAL_T  bdry(DIMV(bdry),2*nvar)
      REAL_T  fine(DIMV(fine),  nvar)

      integer i,j,n

      j = lo(2)-1
      do n=1,nvar
c ::: ::: copy boundary values over
         do i=lo(1),hi(1)
            bdry(i,j,n) = fine(i,j,n)
         enddo
c ::: ::: do stinking calculation of derivativ using centered difference
c ::: ::: this ONLY makes sense for DIRICHLET BC.  For other BC,
c ::: ::: we will not use these values.
         do i=lo(1)+1,hi(1)-1
            bdry(i,j,n+nvar) = (fine(i+1,j,n)-fine(i-1,j,n))/(2*hfine)
         enddo
         i = lo(1)
         bdry(i,j,n+nvar) = (-fine(i+2,j,n)+4*fine(i+1,j,n)-3*fine(i,j,n))/
     &        (2*hfine)
         i = hi(1)
         bdry(i,j,n+nvar) = (+fine(i-2,j,n)-4*fine(i-1,j,n)+3*fine(i,j,n))/
     &        (2*hfine)
      enddo

      return
      end
      


      subroutine FORT_BDIDERIVYHI( bdry, DIMS(bdry),
     &     lo, hi,
     &     fine, DIMS(fine),
     &     nvar, hfine)
      integer nvar
      REAL_T  hfine
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(bdry)
      integer DIMDEC(fine)
      REAL_T  bdry(DIMV(bdry),2*nvar)
      REAL_T  fine(DIMV(fine),  nvar)

      integer i,j,n

      j = hi(2)+1
      do n=1,nvar
c ::: ::: copy boundary values over
         do i=lo(1),hi(1)
            bdry(i,j,n) = fine(i,j,n)
         enddo
c ::: ::: do stinking calculation of derivativ using centered difference
c ::: ::: this ONLY makes sense for DIRICHLET BC.  For other BC,
c ::: ::: we will not use these values.
         do i=lo(1)+1,hi(1)-1
            bdry(i,j,n+nvar) = (fine(i+1,j,n)-fine(i-1,j,n))/(2*hfine)
         enddo
         i = lo(1)
         bdry(i,j,n+nvar) = (-fine(i+2,j,n)+4*fine(i+1,j,n)-3*fine(i,j,n))/
     &        (2*hfine)
         i = hi(1)
         bdry(i,j,n+nvar) = (+fine(i-2,j,n)-4*fine(i-1,j,n)+3*fine(i,j,n))/
     &        (2*hfine)
      enddo

      return
      end
      

