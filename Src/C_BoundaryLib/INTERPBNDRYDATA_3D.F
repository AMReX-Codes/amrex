c
c $Id: INTERPBNDRYDATA_3D.F,v 1.2 2002-10-19 18:33:57 car Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "BC_TYPES.H"
#include "REAL.H"
#include "CONSTANTS.H"
#include "INTERPBNDRYDATA_F.H"
#include "ArrayLim.H"

#define SDIM 3
#define NUMDERIV 5
#define XDER   1
#define YDER   2
#define X2DER  3
#define Y2DER  4
#define XYDER  5

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
c ::  nvar        => number of variables to interpolate
c ::  ratios(3)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPXLO (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratios,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives)

      integer  nvar, ratios(3), not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(cb)
      integer  DIMDEC(mask)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),nvar)
      REAL_T   derives(DIM23(cb),NUMDERIV)
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)

      REAL_T   xx, yy, ainterp
      integer  i, j, k, ic, jc, kc, joff, koff, n
      integer  jclo, jchi, kclo, kchi, ratioy, ratioz

      ratioy = ratios(2)
      ratioz = ratios(3)

      kclo = ARG_L3(cb)
      kchi = ARG_H3(cb)
      jclo = ARG_L2(cb)
      jchi = ARG_H2(cb)
      ic   = ARG_L1(cb)-1
      i    = lo(1)-1

      do n = 1, nvar
c        ::::: define interp coefs
         do kc = kclo, kchi
         k = ratioz*kc
         do jc = jclo, jchi
            j = ratioy*jc
            derives(jc,kc,XDER)  = half*(crse(ic,jc+1,kc,n) - crse(ic,jc-1,kc,n))
            derives(jc,kc,X2DER) = half*(crse(ic,jc+1,kc,n) - two*crse(ic,jc,kc,n)
     &                     + crse(ic,jc-1,kc,n))
            derives(jc,kc,YDER)  = half*(crse(ic,jc,kc+1,n) - crse(ic,jc,kc-1,n))
            derives(jc,kc,Y2DER) = half*(crse(ic,jc,kc+1,n) - two*crse(ic,jc,kc,n)
     &                     + crse(ic,jc,kc-1,n))
            derives(jc,kc,XYDER) = forth*(crse(ic,jc+1,kc+1,n) - crse(ic,jc-1,kc+1,n)
     &                     + crse(ic,jc-1,kc-1,n) - crse(ic,jc+1,kc-1,n))

            if (mask(i,j-1,k) .ne. not_covered) then
               derives(jc,kc,XDER)  = crse(ic,jc+1,kc,n) - crse(ic,jc,kc,n)
               derives(jc,kc,X2DER) = zero
            end if
            if (mask(i,j+ratioy,k) .ne. not_covered) then
               derives(jc,kc,XDER)  = crse(ic,jc,kc,n) - crse(ic,jc-1,kc,n)
               derives(jc,kc,X2DER) = zero
            end if
            if (mask(i,j-1,k)     .ne. not_covered .and.
     $          mask(i,j+ratioy,k) .ne. not_covered) then
               derives(jc,kc,XDER)  = zero
            end if

            if (mask(i,j,k-1) .ne. not_covered) then
               derives(jc,kc,YDER)  = crse(ic,jc,kc+1,n) - crse(ic,jc,kc,n)
               derives(jc,kc,Y2DER) = zero
            end if
            if (mask(i,j,k+ratioz) .ne. not_covered) then
               derives(jc,kc,YDER)  = crse(ic,jc,kc,n) - crse(ic,jc,kc-1,n)
               derives(jc,kc,Y2DER) = zero
            end if
            if (mask(i,j,k-1)     .ne. not_covered .and.
     $          mask(i,j,k+ratioz) .ne. not_covered) then
               derives(jc,kc,YDER)  = zero
            end if
            if ( 
     $           ( mask(i,j+ratioy,k+ratioz) .ne. not_covered ) .or.
     $           ( mask(i,j-1,k+ratioz)     .ne. not_covered ) .or.
     $           ( mask(i,j+ratioy,k-1)     .ne. not_covered ) .or.
     $           ( mask(i,j-1,k-1)         .ne. not_covered ) ) then

               derives(jc,kc,XYDER) = zero
            end if
         end do
         end do

c        ::::: interpolate to fine grid
         do koff = 0, ratioz - 1
            yy = (dble(koff - ratioz/2) + half)/ratioz
            do kc = kclo,kchi
               k = ratioz*kc + koff
               do joff = 0, ratioy - 1
                  xx = (dble(joff - ratioy/2) + half)/ratioy
                  do jc = jclo, jchi
                     j = ratioy*jc + joff
                     bdry(i,j,k,n) = crse(ic,jc,kc,n) + xx*derives(jc,kc,XDER) 
     &                     + derives(jc,kc,X2DER)*xx**2 + yy*derives(jc,kc,YDER) 
     &                     + derives(jc,kc,Y2DER)*yy**2 + xx*yy*derives(jc,kc,XYDER) 
                  end do
               end do
            end do
         end do
      end do

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
c ::  nvar        => number of variables to interpolate
c ::  ratios(3)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPXHI (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratios,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives)

      integer  nvar, ratios(3), not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(cb)
      integer  DIMDEC(mask)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),nvar)
      REAL_T   derives(DIM23(cb),NUMDERIV)
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)

      REAL_T   xx, yy, ainterp
      integer  i, j, k, ic, jc, kc, joff, koff, n
      integer  jclo, jchi, kclo, kchi, ratioy, ratioz

      ratioy = ratios(2)
      ratioz = ratios(3)

      kclo = ARG_L3(cb)
      kchi = ARG_H3(cb)
      jclo = ARG_L2(cb)
      jchi = ARG_H2(cb)
      ic   = ARG_H1(cb)+1
      i    = hi(1)+1

      do n = 1, nvar
c        ::::: define interp coefs
         do kc = kclo, kchi
         k = ratioz*kc
         do jc = jclo, jchi
            j = ratioy*jc
            derives(jc,kc,XDER)  = half*(crse(ic,jc+1,kc,n) - crse(ic,jc-1,kc,n))
            derives(jc,kc,X2DER) = half*(crse(ic,jc+1,kc,n) - two*crse(ic,jc,kc,n)
     &                     + crse(ic,jc-1,kc,n))
            derives(jc,kc,YDER)  = half*(crse(ic,jc,kc+1,n) - crse(ic,jc,kc-1,n))
            derives(jc,kc,Y2DER) = half*(crse(ic,jc,kc+1,n) - two*crse(ic,jc,kc,n)
     &                     + crse(ic,jc,kc-1,n))
            derives(jc,kc,XYDER) = forth*(crse(ic,jc+1,kc+1,n) - crse(ic,jc-1,kc+1,n)
     &                     + crse(ic,jc-1,kc-1,n) - crse(ic,jc+1,kc-1,n))

            if (mask(i,j-1,k) .ne. not_covered) then
               derives(jc,kc,XDER)  = crse(ic,jc+1,kc,n) - crse(ic,jc,kc,n)
               derives(jc,kc,X2DER) = zero
            end if
            if (mask(i,j+ratioy,k) .ne. not_covered) then
               derives(jc,kc,XDER)  = crse(ic,jc,kc,n) - crse(ic,jc-1,kc,n)
               derives(jc,kc,X2DER) = zero
            end if
            if (mask(i,j-1,k)     .ne. not_covered .and.
     $          mask(i,j+ratioy,k) .ne. not_covered) then
               derives(jc,kc,XDER)  = zero
            end if

            if (mask(i,j,k-1) .ne. not_covered) then
               derives(jc,kc,YDER)  = crse(ic,jc,kc+1,n) - crse(ic,jc,kc,n)
               derives(jc,kc,Y2DER) = zero
            end if
            if (mask(i,j,k+ratioz) .ne. not_covered) then
               derives(jc,kc,YDER)  = crse(ic,jc,kc,n) - crse(ic,jc,kc-1,n)
               derives(jc,kc,Y2DER) = zero
            end if
            if (mask(i,j,k-1)     .ne. not_covered .and.
     $          mask(i,j,k+ratioz) .ne. not_covered) then
               derives(jc,kc,YDER)  = zero
            end if
            if (
     $           ( mask(i,j+ratioy,k+ratioz) .ne. not_covered ) .or.
     $           ( mask(i,j-1,k+ratioz)     .ne. not_covered ) .or.
     $           ( mask(i,j+ratioy,k-1)     .ne. not_covered ) .or.
     $           ( mask(i,j-1,k-1)         .ne. not_covered ) ) then
               
                derives(jc,kc,XYDER) = zero
            end if
         end do
         end do

c        ::::: interpolate to fine grid
         do koff = 0, ratioz - 1
            yy = (dble(koff - ratioz/2) + half)/ratioz
            do kc = kclo,kchi
               k = ratioz*kc + koff
               do joff = 0, ratioy - 1
                  xx = (dble(joff - ratioy/2) + half)/ratioy
                  do jc = jclo, jchi
                     j = ratioy*jc + joff
                     bdry(i,j,k,n) = crse(ic,jc,kc,n) + xx*derives(jc,kc,XDER) 
     &                     + derives(jc,kc,X2DER)*xx**2 + yy*derives(jc,kc,YDER) 
     &                     + derives(jc,kc,Y2DER)*yy**2 + xx*yy*derives(jc,kc,XYDER) 
                  end do
               end do
            end do
         end do
      end do

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
c ::  nvar        => number of variables to interpolate
c ::  ratios(3)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPYLO (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratios,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives)

      integer  nvar, ratios(3), not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(cb)
      integer  DIMDEC(mask)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),nvar)
      REAL_T   derives(DIM13(cb),NUMDERIV)
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)

      REAL_T   xx, yy, ainterp
      integer  i, j, k, ic, jc, kc, ioff, koff, n
      integer  iclo, ichi, kclo, kchi, ratiox, ratioz

      ratiox = ratios(1)
      ratioz = ratios(3)

      kclo = ARG_L3(cb)
      kchi = ARG_H3(cb)
      iclo = ARG_L1(cb)
      ichi = ARG_H1(cb)
      jc   = ARG_L2(cb)-1
      j    = lo(2)-1

      do n = 1, nvar
c        ::::: define interp coefs
         do kc = kclo, kchi
         k = ratioz*kc
         do ic = iclo, ichi
            i = ratiox*ic
            derives(ic,kc,XDER)  = half*(crse(ic+1,jc,kc,n) - crse(ic-1,jc,kc,n))
            derives(ic,kc,X2DER) = half*(crse(ic+1,jc,kc,n) - two*crse(ic,jc,kc,n)
     &                     + crse(ic-1,jc,kc,n))
            derives(ic,kc,YDER)  = half*(crse(ic,jc,kc+1,n) - crse(ic,jc,kc-1,n))
            derives(ic,kc,Y2DER) = half*(crse(ic,jc,kc+1,n) - two*crse(ic,jc,kc,n)
     &                     + crse(ic,jc,kc-1,n))
            derives(ic,kc,XYDER) = forth*(crse(ic+1,jc,kc+1,n) - crse(ic-1,jc,kc+1,n)
     &                     + crse(ic-1,jc,kc-1,n) - crse(ic+1,jc,kc-1,n))

            if (mask(i-1,j,k) .ne. not_covered) then
               derives(ic,kc,XDER)  = crse(ic+1,jc,kc,n) - crse(ic,jc,kc,n)
               derives(ic,kc,X2DER) = zero
            end if
            if (mask(i+ratiox,j,k) .ne. not_covered) then
               derives(ic,kc,XDER)  = crse(ic,jc,kc,n) - crse(ic-1,jc,kc,n)
               derives(ic,kc,X2DER) = zero
            end if
            if (mask(i-1,j,k)     .ne. not_covered .and.
     $          mask(i+ratiox,j,k) .ne. not_covered) then
               derives(ic,kc,XDER)  = zero
            end if

            if (mask(i,j,k-1) .ne. not_covered) then
               derives(ic,kc,YDER)  = crse(ic,jc,kc+1,n) - crse(ic,jc,kc,n)
               derives(ic,kc,Y2DER) = zero
            end if
            if (mask(i,j,k+ratioz) .ne. not_covered) then
               derives(ic,kc,YDER)  = crse(ic,jc,kc,n) - crse(ic,jc,kc-1,n)
               derives(ic,kc,Y2DER) = zero
            end if
            if (mask(i,j,k-1)     .ne. not_covered .and.
     $          mask(i,j,k+ratioz) .ne. not_covered) then
               derives(ic,kc,YDER)  = zero
            end if
            if (
     $           ( mask(i+ratiox,j,k+ratioz) .ne. not_covered ) .or.
     $           ( mask(i-1,j,k+ratioz)     .ne. not_covered ) .or.
     $           ( mask(i+ratiox,j,k-1)     .ne. not_covered ) .or.
     $           ( mask(i-1,j,k-1)         .ne. not_covered ) ) then
               
                derives(ic,kc,XYDER) = zero
            end if
         end do
         end do

c        ::::: interpolate to fine grid
         do koff = 0, ratioz - 1
            yy = (dble(koff - ratioz/2) + half)/ratioz
            do kc = kclo,kchi
               k = ratioz*kc + koff
               do ioff = 0, ratiox - 1
                  xx = (dble(ioff - ratiox/2) + half)/ratiox
                  do ic = iclo, ichi
                     i = ratiox*ic + ioff
                     bdry(i,j,k,n) = crse(ic,jc,kc,n) + xx*derives(ic,kc,XDER) 
     &                     + derives(ic,kc,X2DER)*xx**2 + yy*derives(ic,kc,YDER) 
     &                     + derives(ic,kc,Y2DER)*yy**2 + xx*yy*derives(ic,kc,XYDER) 
                  end do
               end do
            end do
         end do
      end do

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
c ::  nvar        => number of variables to interpolate
c ::  ratios(3)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPYHI (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratios,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives)

      integer  nvar, ratios(3), not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(cb)
      integer  DIMDEC(mask)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),nvar)
      REAL_T   derives(DIM13(cb),NUMDERIV)
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)

      REAL_T   xx, yy, ainterp
      integer  i, j, k, ic, jc, kc, ioff, koff, n
      integer  iclo, ichi, kclo, kchi, ratiox, ratioz

      ratiox = ratios(1)
      ratioz = ratios(3)

      kclo = ARG_L3(cb)
      kchi = ARG_H3(cb)
      iclo = ARG_L1(cb)
      ichi = ARG_H1(cb)
      jc   = ARG_H2(cb)+1
      j    = hi(2)+1

      do n = 1, nvar
c        ::::: define interp coefs
         do kc = kclo, kchi
         k = ratioz*kc
         do ic = iclo, ichi
            i = ratiox*ic
            derives(ic,kc,XDER)  = half*(crse(ic+1,jc,kc,n) - crse(ic-1,jc,kc,n))
            derives(ic,kc,X2DER) = half*(crse(ic+1,jc,kc,n) - two*crse(ic,jc,kc,n)
     &                     + crse(ic-1,jc,kc,n))
            derives(ic,kc,YDER)  = half*(crse(ic,jc,kc+1,n) - crse(ic,jc,kc-1,n))
            derives(ic,kc,Y2DER) = half*(crse(ic,jc,kc+1,n) - two*crse(ic,jc,kc,n)
     &                     + crse(ic,jc,kc-1,n))
            derives(ic,kc,XYDER) = forth*(crse(ic+1,jc,kc+1,n) - crse(ic-1,jc,kc+1,n)
     &                     + crse(ic-1,jc,kc-1,n) - crse(ic+1,jc,kc-1,n))

            if (mask(i-1,j,k) .ne. not_covered) then
               derives(ic,kc,XDER)  = crse(ic+1,jc,kc,n) - crse(ic,jc,kc,n)
               derives(ic,kc,X2DER) = zero
            end if
            if (mask(i+ratiox,j,k) .ne. not_covered) then
               derives(ic,kc,XDER)  = crse(ic,jc,kc,n) - crse(ic-1,jc,kc,n)
               derives(ic,kc,X2DER) = zero
            end if
            if (mask(i-1,j,k)     .ne. not_covered .and.
     $          mask(i+ratiox,j,k) .ne. not_covered) then
               derives(ic,kc,XDER)  = zero
            end if

            if (mask(i,j,k-1) .ne. not_covered) then
               derives(ic,kc,YDER)  = crse(ic,jc,kc+1,n) - crse(ic,jc,kc,n)
               derives(ic,kc,Y2DER) = zero
            end if
            if (mask(i,j,k+ratioz) .ne. not_covered) then
               derives(ic,kc,YDER)  = crse(ic,jc,kc,n) - crse(ic,jc,kc-1,n)
               derives(ic,kc,Y2DER) = zero
            end if
            if (mask(i,j,k-1)     .ne. not_covered .and.
     $          mask(i,j,k+ratioz) .ne. not_covered) then
               derives(ic,kc,YDER)  = zero
            end if

            if ( 
     $           ( mask(i+ratiox,j,k+ratioz) .ne. not_covered ) .or.
     $           ( mask(i-1,j,k+ratioz)     .ne. not_covered ) .or.
     $           ( mask(i+ratiox,j,k-1)     .ne. not_covered ) .or.
     $           ( mask(i-1,j,k-1)         .ne. not_covered ) ) then

               derives(ic,kc,XYDER) = zero
            end if
         end do
         end do

c        ::::: interpolate to fine grid
         do koff = 0, ratioz - 1
            yy = (dble(koff - ratioz/2) + half)/ratioz
            do kc = kclo,kchi
               k = ratioz*kc + koff
               do ioff = 0, ratiox - 1
                  xx = (dble(ioff - ratiox/2) + half)/ratiox
                  do ic = iclo, ichi
                     i = ratiox*ic + ioff
                     bdry(i,j,k,n) = crse(ic,jc,kc,n) + xx*derives(ic,kc,XDER) 
     &                     + derives(ic,kc,X2DER)*xx**2 + yy*derives(ic,kc,YDER) 
     &                     + derives(ic,kc,Y2DER)*yy**2 + xx*yy*derives(ic,kc,XYDER) 
                  end do
               end do
            end do
         end do
      end do

      return
      end

c ---------------------------------------------------------------
c ::  FORT_BDINTERPZLO : Interpolation on Zlo Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  DIMS(bdry)  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  DIMS(cb)    => index limits of coarsened grid interior
c ::  nvar        => number of variables to interpolate
c ::  ratios(3)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPZLO (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratios,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives)

      integer  nvar, ratios(3), not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(cb)
      integer  DIMDEC(mask)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),nvar)
      REAL_T   derives(DIM12(cb),NUMDERIV)
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)

      REAL_T   xx, yy, ainterp
      integer  i, j, k, ic, jc, kc, ioff, joff, n
      integer  iclo, ichi, jclo, jchi, ratiox, ratioy

      ratiox = ratios(1)
      ratioy = ratios(2)

      jclo = ARG_L2(cb)
      jchi = ARG_H2(cb)
      iclo = ARG_L1(cb)
      ichi = ARG_H1(cb)
      kc   = ARG_L3(cb)-1
      k    = lo(3)-1

      do n = 1, nvar
c     ::::: define interp coefs
         do jc = jclo, jchi
            j = ratioy*jc
            do ic = iclo, ichi
               i = ratiox*ic
               
               derives(ic,jc,XDER)  = half*(crse(ic+1,jc,kc,n) - crse(ic-1,jc,kc,n) )
               derives(ic,jc,X2DER) = half*(crse(ic+1,jc,kc,n) - two*crse(ic,jc,kc,n) + crse(ic-1,jc,kc,n) )
               derives(ic,jc,YDER)  = half*(crse(ic,jc+1,kc,n) - crse(ic,jc-1,kc,n) )
               derives(ic,jc,Y2DER) = half*(crse(ic,jc+1,kc,n) - two*crse(ic,jc,kc,n) + crse(ic,jc-1,kc,n) )
               derives(ic,jc,XYDER) = forth*(crse(ic+1,jc+1,kc,n) - crse(ic-1,jc+1,kc,n)
     &              + crse(ic-1,jc-1,kc,n) - crse(ic+1,jc-1,kc,n))
               
               if (mask(i-1,j,k) .ne. not_covered) then
                  derives(ic,jc,XDER)  = crse(ic+1,jc,kc,n) - crse(ic,jc,kc,n)
                  derives(ic,jc,X2DER) = zero
               end if
               if (mask(i+ratiox,j,k) .ne. not_covered) then
                  derives(ic,jc,XDER)  = crse(ic,jc,kc,n) - crse(ic-1,jc,kc,n)
                  derives(ic,jc,X2DER) = zero
               end if
               if (mask(i-1,j,k)     .ne. not_covered .and. 
     $             mask(i+ratiox,j,k) .ne. not_covered) then
                  derives(ic,jc,XDER)  = zero
               end if

               if (mask(i,j-1,k) .ne. not_covered) then
                  derives(ic,jc,YDER)  = crse(ic,jc+1,kc,n) - crse(ic,jc,kc,n)
                  derives(ic,jc,Y2DER) = zero
               end if
               if (mask(i,j+ratioy,k) .ne. not_covered) then
                  derives(ic,jc,YDER)  = crse(ic,jc,kc,n) - crse(ic,jc-1,kc,n)
                  derives(ic,jc,Y2DER) = zero
               end if
               if (mask(i,j-1,k)     .ne. not_covered .and. 
     $             mask(i,j+ratioy,k) .ne. not_covered) then
                  derives(ic,jc,YDER)  = zero
               end if

               if (
     $              ( mask(i+ratiox,j+ratioy,k) .ne. not_covered ) .or.
     $              ( mask(i-1,j+ratioy,k)     .ne. not_covered ) .or.
     $              ( mask(i+ratiox,j-1,k)     .ne. not_covered ) .or.
     $              ( mask(i-1,j-1,k)         .ne. not_covered ) ) then
                  
                  derives(ic,jc,XYDER) = zero
               end if
            end do
         end do
         
c     ::::: interpolate to fine grid
         do joff = 0, ratioy - 1
            yy = (dble(joff - ratioy/2) + half)/ratioy
            do jc = jclo,jchi
               j = ratioy*jc + joff
               do ioff = 0, ratiox - 1
                  xx = (dble(ioff - ratiox/2) + half)/ratiox
                  do ic = iclo, ichi
                     i = ratiox*ic + ioff
                     bdry(i,j,k,n) = crse(ic,jc,kc,n) + xx*derives(ic,jc,XDER) 
     &                    + derives(ic,jc,X2DER)*xx**2 + yy*derives(ic,jc,YDER) 
     &                    + derives(ic,jc,Y2DER)*yy**2 + xx*yy*derives(ic,jc,XYDER) 
                  end do
               end do
            end do
         end do
      end do
      
      return
      end
      
c ---------------------------------------------------------------
c ::  FORT_BDINTERPZHI : Interpolation on Zhi Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  DIMS(bdry)  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  DIMS(cb)    => index limits of coarsened grid interior
c ::  nvar        => number of variables to interpolate
c ::  ratios(3)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array for derivatives
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPZHI (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratios,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives)

      integer  nvar, ratios(3), not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(cb)
      integer  DIMDEC(mask)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),nvar)
      REAL_T   derives(DIM12(cb),NUMDERIV)
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)

      REAL_T   xx, yy, ainterp
      integer  i, j, k, ic, jc, kc, ioff, joff, n
      integer  iclo, ichi, jclo, jchi, ratiox, ratioy

      ratiox = ratios(1)
      ratioy = ratios(2)

      jclo = ARG_L2(cb)
      jchi = ARG_H2(cb)
      iclo = ARG_L1(cb)
      ichi = ARG_H1(cb)
      kc   = ARG_H3(cb)+1
      k    = hi(3)+1

      do n = 1, nvar
c        ::::: define interp coefs
         do jc = jclo, jchi
         j = ratioy*jc
         do ic = iclo, ichi
            i = ratiox*ic
            derives(ic,jc,XDER)  = half*(crse(ic+1,jc,kc,n) - crse(ic-1,jc,kc,n))
            derives(ic,jc,X2DER) = half*(crse(ic+1,jc,kc,n) - two*crse(ic,jc,kc,n)
     &                     + crse(ic-1,jc,kc,n))
            derives(ic,jc,YDER)  = half*(crse(ic,jc+1,kc,n) - crse(ic,jc-1,kc,n))
            derives(ic,jc,Y2DER) = half*(crse(ic,jc+1,kc,n) - two*crse(ic,jc,kc,n)
     &                     + crse(ic,jc-1,kc,n))
            derives(ic,jc,XYDER) = forth*(crse(ic+1,jc+1,kc,n) - crse(ic-1,jc+1,kc,n)
     &                     + crse(ic-1,jc-1,kc,n) - crse(ic+1,jc-1,kc,n))

            if (mask(i-1,j,k) .ne. not_covered) then
               derives(ic,jc,XDER)  = crse(ic+1,jc,kc,n) - crse(ic,jc,kc,n)
               derives(ic,jc,X2DER) = zero
            end if
            if (mask(i+ratiox,j,k) .ne. not_covered) then
               derives(ic,jc,XDER)  = crse(ic,jc,kc,n) - crse(ic-1,jc,kc,n)
               derives(ic,jc,X2DER) = zero
            end if
            if (mask(i-1,j,k)     .ne. not_covered .and. 
     $          mask(i+ratiox,j,k) .ne. not_covered) then
               derives(ic,jc,XDER)  = zero
            end if

            if (mask(i,j-1,k) .ne. not_covered) then
               derives(ic,jc,YDER)  = crse(ic,jc+1,kc,n) - crse(ic,jc,kc,n)
               derives(ic,jc,Y2DER) = zero
            end if
            if (mask(i,j+ratioy,k) .ne. not_covered) then
               derives(ic,jc,YDER)  = crse(ic,jc,kc,n) - crse(ic,jc-1,kc,n)
               derives(ic,jc,Y2DER) = zero
            end if
            if (mask(i,j-1,k)     .ne. not_covered .and. 
     $          mask(i,j+ratioy,k) .ne. not_covered) then
               derives(ic,jc,YDER)  = zero
            end if
            if (
     $           ( mask(i+ratiox,j+ratioy,k) .ne. not_covered ) .or.
     $           ( mask(i-1,j+ratioy,k)     .ne. not_covered ) .or.
     $           ( mask(i+ratiox,j-1,k)     .ne. not_covered ) .or.
     $           ( mask(i-1,j-1,k)         .ne. not_covered ) ) then
               
                derives(ic,jc,XYDER) = zero
            end if
         end do
         end do

c        ::::: interpolate to fine grid
         do joff = 0, ratioy - 1
            yy = (dble(joff - ratioy/2) + half)/ratioy
            do jc = jclo,jchi
               j = ratioy*jc + joff
               do ioff = 0, ratiox - 1
                  xx = (dble(ioff - ratiox/2) + half)/ratiox
                  do ic = iclo, ichi
                     i = ratiox*ic + ioff
                     bdry(i,j,k,n) = crse(ic,jc,kc,n) + xx*derives(ic,jc,XDER) 
     &                     + derives(ic,jc,X2DER)*xx**2 + yy*derives(ic,jc,YDER) 
     &                     + derives(ic,jc,Y2DER)*yy**2 + xx*yy*derives(ic,jc,XYDER) 
                  end do
               end do
            end do
         end do
      end do

      return
      end

#undef NUMDERIV
#undef XDER
#undef YDER
#undef X2DER
#undef Y2DER
#undef XYDER

