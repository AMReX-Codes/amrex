c
c $Id: INTERPBNDRYDATA_2D.F,v 1.2 2001-01-19 22:57:48 marc Exp $
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

#define USE_REORG_VERSION
      
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
c ::  nvar        => number of variables to interpolate
c ::  ratios(2)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPXLO (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratios,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives)

      integer  nvar, ratios(2), not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(mask)
      integer  DIMDEC(crse)
      integer  DIMDEC(cb)
      REAL_T   bdry(DIMV(bdry),nvar)
      REAL_T   derives(DIM2(cb),NUMDERIV)      
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)

      REAL_T   xx, ainterp
      integer  i, j, ic, jc, off, n
      integer  jclo, jchi, ratioy

#ifndef USE_REORG_VERSION
      
      ratioy = ratios(2)

      jclo = ARG_L2(cb)
      jchi = ARG_H2(cb)
      ic   = ARG_L1(cb)-1
      i    = lo(1)-1

      do n = 1, nvar
c        ::::: define interp coefs
         do jc = jclo, jchi
            j = ratioy*jc
            derives(jc,XDER)  = half * (crse(ic,jc+1,n) - crse(ic,jc-1,n))
            derives(jc,X2DER) = crse(ic,jc+1,n)
     &           -              crse(ic,jc  ,n) * two
     &           +              crse(ic,jc-1,n)

            if (mask(i,j-1) .ne. not_covered) then
               derives(jc,XDER)  = crse(ic,jc+1,n) - crse(ic,jc  ,n)
               derives(jc,X2DER) = zero
c     Can we do better by sliding stencil up?
               if (mask(i,j+2*ratioy) .eq. not_covered) then
                  derives(jc,XDER)  = half*(-three*crse(ic,jc,n)
     &           +                           four*crse(ic,jc+1,n)
     &           -                           crse(ic,jc+2,n))
                  derives(jc,X2DER) = crse(ic,jc+2,n)
     &           -                    crse(ic,jc+1,n) * two
     &           +                    crse(ic  ,jc,n)
               end if
            end if
            if (mask(i,j+ratioy) .ne. not_covered) then
               derives(jc,XDER)  = crse(ic,jc  ,n) - crse(ic,jc-1,n)
               derives(jc,X2DER) = zero
c     Can we do better by sliding stencil down?
               if (mask(i,j+ratioy) .eq. not_covered) then
                  derives(jc,XDER)  = half*(three*crse(ic,jc,n)
     &           -                          four*crse(ic,jc-1,n)
     &           +                          crse(ic,jc-2,n))
                  derives(jc,X2DER) = crse(ic,jc,n)
     &           -                    crse(ic,jc-1,n) * two
     &           +                    crse(ic,jc-2,n)
               end if
            end if
            if (mask(i,j-1)     .ne. not_covered .and. 
     $          mask(i,j+ratioy) .ne. not_covered) then
               derives(jc,XDER)  = zero
            end if
         end do

c        ::::: interpolate to fine grid
         do off = 0, ratioy - 1
            xx = (dfloat(off - ratioy/2) + half)/ratioy
            do jc = jclo, jchi
               j = ratioy*jc + off
               bdry(i,j,n) = crse(ic,jc,n)
     &              +        derives(jc,XDER) *xx
     &              +        half*derives(jc,X2DER)*xx**2
            end do
         end do
      end do

#else
      integer Norder, NN, m
      parameter (Norder = 3)
      REAL_T x(Norder), y(Norder), c(Norder), xInt
      ratioy = ratios(2)

      jclo = ARG_L2(cb)
      jchi = ARG_H2(cb)
      ic   = ARG_L1(cb)-1
      i    = lo(1)-1
      
      do n=1,nvar
         do jc=jclo,jchi
            j = ratioy*jc
            
            NN = 1
            y(NN) = crse(ic,jc,n)
            x(NN) = zero

            if (mask(i,j-1).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic,jc-1,n)
               x(NN) = -one
            else if (mask(i,j+2*ratioy).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic,jc+2,n)
               x(NN) = two
            endif
            
            if (mask(i,j+ratioy).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic,jc+1,n)
               x(NN) = one
            else if (mask(i,j-ratioy-1).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic,jc-2,n)
               x(NN) = -two
            endif
               
            if ( (mask(i,j-1).ne.not_covered).and.
     &           (mask(i,j+ratioy).ne.not_covered) ) NN = 1
            
            do off = 0,ratioy-1
               xInt = (dfloat(off - ratioy/2) + half)/ratioy
               call polyInterpCoeff(xInt, x, NN, c)
               bdry(i,j+off,n) = zero
               do m=1,NN
                  bdry(i,j+off,n) = bdry(i,j+off,n) + c(m)*y(m)
               end do
            end do
         end do
      end do
      
#endif
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
c ::  ratios(2)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPXHI (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratios,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives)

      integer  nvar, ratios(2), not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(mask)
      integer  DIMDEC(cb)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),nvar)
      REAL_T   derives(DIM2(cb),NUMDERIV)      
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)

      REAL_T   xx, ainterp
      integer  i, j, ic, jc, off, n
      integer  jclo, jchi, ratioy

#ifndef USE_REORG_VERSION
      
      ratioy = ratios(2)

      jclo = ARG_L2(cb)
      jchi = ARG_H2(cb)
      ic   = ARG_H1(cb)+1
      i    = hi(1)+1
      
      do n = 1, nvar
c        ::::: define interp coefs
         do jc = jclo, jchi
            j = ratioy*jc
            derives(jc,XDER)  = half * (crse(ic,jc+1,n) - crse(ic,jc-1,n))
            derives(jc,X2DER) = crse(ic,jc+1,n)
     &           -              crse(ic,jc  ,n) * two
     &           +              crse(ic,jc-1,n)

            if (mask(i,j-1) .ne. not_covered) then
               derives(jc,XDER)  = crse(ic,jc+1,n) - crse(ic,jc  ,n)
               derives(jc,X2DER) = zero
c     Can we do better by sliding stencil up?
               if (mask(i,j+2*ratioy) .eq. not_covered) then
                  derives(jc,XDER)  = half*(-three*crse(ic,jc,n)
     &           +                           four*crse(ic,jc+1,n)
     &           -                           crse(ic,jc+2,n))
                  derives(jc,X2DER) = crse(ic,jc+2,n)
     &           -                    crse(ic,jc+1,n) * two
     &           +                    crse(ic  ,jc,n)
               end if
            end if
            if (mask(i,j+ratioy) .ne. not_covered) then
               derives(jc,XDER)  = crse(ic,jc  ,n) - crse(ic,jc-1,n)
               derives(jc,X2DER) = zero
c     Can we do better by sliding stencil down?
               if (mask(i,j+ratioy) .eq. not_covered) then
                  derives(jc,XDER)  = half*(three*crse(ic,jc,n)
     &           -                          four*crse(ic,jc-1,n)
     &           +                          crse(ic,jc-2,n))
                  derives(jc,X2DER) = crse(ic,jc,n)
     &           -                    crse(ic,jc-1,n) * two
     &           +                    crse(ic,jc-2,n)
               end if
            end if
            if (mask(i,j-1) .ne. not_covered .and.
     $          mask(i,j+ratioy) .ne. not_covered) then
               derives(jc,XDER)  = zero
               derives(jc,X2DER)  = zero
            end if
         end do

c        ::::: interpolate to fine grid
         do off = 0, ratioy - 1
            xx = (dfloat(off - ratioy/2) + half)/ratioy
            do jc = jclo, jchi
               j = ratioy*jc + off
               bdry(i,j,n) = crse(ic,jc,n)
     &              +        derives(jc,XDER) *xx
     &              +        half*derives(jc,X2DER)*xx**2
            end do
         end do
      end do

#else
      integer Norder, NN, m
      parameter (Norder = 3)
      REAL_T x(Norder), y(Norder), c(Norder), xInt
      ratioy = ratios(2)

      jclo = ARG_L2(cb)
      jchi = ARG_H2(cb)
      ic   = ARG_H1(cb)+1
      i    = hi(1)+1
      
      do n=1,nvar
         do jc=jclo,jchi
            j = ratioy*jc
            
            NN = 1
            y(NN) = crse(ic,jc,n)
            x(NN) = zero

            if (mask(i,j-1).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic,jc-1,n)
               x(NN) = -one
            else if (mask(i,j+2*ratioy).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic,jc+2,n)
               x(NN) = two
            endif
            
            if (mask(i,j+ratioy).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic,jc+1,n)
               x(NN) = one
            else if (mask(i,j-ratioy-1).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic,jc-2,n)
               x(NN) = -two
            endif
               
            if ( (mask(i,j-1).ne.not_covered).and.
     &           (mask(i,j+ratioy).ne.not_covered) ) NN = 1
            
            do off = 0,ratioy-1
               xInt = (dfloat(off - ratioy/2) + half)/ratioy
               call polyInterpCoeff(xInt, x, NN, c)
               bdry(i,j+off,n) = zero
               do m=1,NN
                  bdry(i,j+off,n) = bdry(i,j+off,n) + c(m)*y(m)
               end do
            end do
         end do
      end do
      
#endif
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
c ::  ratios(2)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPYLO (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratios,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives)

      integer  nvar, ratios(2), not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(mask)
      integer  DIMDEC(cb)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),nvar)
      REAL_T   derives(DIM1(cb),NUMDERIV)
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T   xx, ainterp
      integer  i, j, ic, jc, off, n
      integer  iclo, ichi, ratiox

#ifndef USE_REORG_VERSION
      
      ratiox = ratios(1)

      iclo = ARG_L1(cb)
      ichi = ARG_H1(cb)
      jc   = ARG_L2(cb)-1
      j    = lo(2)-1

      do n = 1, nvar
c        ::::: define interp coefs
         do ic = iclo, ichi
            i = ratiox*ic
            derives(ic,XDER)  = half * (crse(ic+1,jc,n) - crse(ic-1,jc,n))
            derives(ic,X2DER) = crse(ic+1,jc,n)
     &           -              crse(ic  ,jc,n) * two
     &           +              crse(ic-1,jc,n)
            if (mask(i-1,j) .ne. not_covered) then
               derives(ic,XDER)  = crse(ic+1,jc,n) - crse(ic  ,jc,n)
               derives(ic,X2DER) = zero
c     Can we do better by sliding stencil to the right
               if (mask(i+2*ratiox,j) .eq. not_covered) then
                  derives(ic,XDER)  = half*(-three*crse(ic,jc,n)
     &           +                           four*crse(ic+1,jc,n)
     &           -                           crse(ic+2,jc,n))
                  derives(ic,X2DER) = crse(ic+2,jc,n)
     &           -                    crse(ic+1,jc,n) * two
     &           +                    crse(ic  ,jc,n)
               end if
            end if
            if (mask(i+ratiox,j) .ne. not_covered) then
               derives(ic,XDER)  = crse(ic  ,jc,n) - crse(ic-1,jc,n)
               derives(ic,X2DER) = zero
            end if
c     Can we do better by sliding stencil to the left
               if (mask(i-ratiox,j) .eq. not_covered) then
                  derives(ic,XDER)  = half*(three*crse(ic,jc,n)
     &           -                          four*crse(ic-1,jc,n)
     &           +                          crse(ic-2,jc,n))
                  derives(ic,X2DER) = crse(ic,jc,n)
     &           -                    crse(ic-1,jc,n) * two
     &           +                    crse(ic-2,jc,n)
               end if
            if (mask(i-1,j) .ne. not_covered .and.
     $          mask(i+ratiox,j) .ne. not_covered) then
               derives(ic,XDER)  = zero
               derives(ic,X2DER)  = zero
            end if
         end do

c        ::::: interpolate to fine grid
         do off = 0, ratiox - 1
            xx = (dfloat(off - ratiox/2) + half)/ratiox
            do ic = iclo, ichi
               i = ratiox*ic + off
               bdry(i,j,n) = crse(ic,jc,n)
     &              +        derives(ic,XDER) *xx
     &              +        half*derives(ic,X2DER)*xx**2
            end do
         end do
      end do

#else
      integer Norder, NN, m
      parameter (Norder = 3)
      REAL_T x(Norder), y(Norder), c(Norder), xInt
      ratiox = ratios(1)

      iclo = ARG_L1(cb)
      ichi = ARG_H1(cb)
      jc   = ARG_L2(cb)-1
      j    = lo(2)-1
      
      do n=1,nvar
         do ic=iclo,ichi
            i = ratiox*ic
            
            NN = 1
            y(NN) = crse(ic,jc,n)
            x(NN) = zero

            if (mask(i-1,j).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic-1,jc,n)
               x(NN) = -one
            else if (mask(i+2*ratiox,j).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic+2,jc,n)
               x(NN) = two
            endif
            
            if (mask(i+ratiox,j).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic+1,jc,n)
               x(NN) = one
            else if (mask(i-ratiox-1,j).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic-2,jc,n)
               x(NN) = -two
            endif
               
            if ( (mask(i-1,j).ne.not_covered).and.
     &           (mask(i+ratiox,j).ne.not_covered) ) NN = 1
            
            do off = 0,ratiox-1
               xInt = (dfloat(off - ratiox/2) + half)/ratiox
               call polyInterpCoeff(xInt, x, NN, c)
               bdry(i+off,j,n) = zero
               do m=1,NN
                  bdry(i+off,j,n) = bdry(i+off,j,n) + c(m)*y(m)
               end do
            end do
         end do
      end do
      
#endif
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
c ::  ratios(2)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  DIMS(mask)  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  DIMS(crse)  => index limits of crse array
c ::  derives     => crse grid tmp array
c ---------------------------------------------------------------

      subroutine FORT_BDINTERPYHI (bdry,DIMS(bdry),
     &           lo,hi,DIMS(cb),nvar,ratios,not_covered,
     &           mask,DIMS(mask),crse,DIMS(crse),derives)

      integer  nvar, ratios(2), not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(bdry)
      integer  DIMDEC(mask)
      integer  DIMDEC(cb)
      integer  DIMDEC(crse)
      REAL_T   bdry(DIMV(bdry),nvar)
      REAL_T   derives(DIM1(cb),NUMDERIV)
      integer  mask(DIMV(mask))
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T   xx, ainterp
      integer  i, j, ic, jc, off, n
      integer  iclo, ichi, ratiox

#ifndef USE_REORG_VERSION
      
      ratiox = ratios(1)

      iclo = ARG_L1(cb)
      ichi = ARG_H1(cb)
      jc   = ARG_H2(cb)+1
      j    = hi(2)+1

      do n = 1, nvar
c        ::::: define interp coefs
         do ic = iclo, ichi
            i = ratiox*ic
            derives(ic,XDER)  = half * (crse(ic+1,jc,n) - crse(ic-1,jc,n))
            derives(ic,X2DER) = crse(ic+1,jc,n)
     &           -              crse(ic  ,jc,n) * two
     &           +              crse(ic-1,jc,n)
            if (mask(i-1,j) .ne. not_covered) then
               derives(ic,XDER)  = crse(ic+1,jc,n) - crse(ic  ,jc,n)
               derives(ic,X2DER) = zero
c     Can we do better by sliding stencil to the right
               if (mask(i+2*ratiox,j) .eq. not_covered) then
                  derives(ic,XDER)  = half*(-three*crse(ic,jc,n)
     &           +                           four*crse(ic+1,jc,n)
     &           -                           crse(ic+2,jc,n))
                  derives(ic,X2DER) = crse(ic+2,jc,n)
     &           -                    crse(ic+1,jc,n) * two
     &           +                    crse(ic  ,jc,n)
               end if
            end if
            if (mask(i+ratiox,j) .ne. not_covered) then
               derives(ic,XDER)  = crse(ic  ,jc,n) - crse(ic-1,jc,n)
               derives(ic,X2DER) = zero
c     Can we do better by sliding stencil to the left
               if (mask(i-ratiox,j) .eq. not_covered) then
                  derives(ic,XDER)  = half*(three*crse(ic,jc,n)
     &           -                          four*crse(ic-1,jc,n)
     &           +                          crse(ic-2,jc,n))
                  derives(ic,X2DER) = crse(ic,jc,n)
     &           -                    crse(ic-1,jc,n) * two
     &           +                    crse(ic-2,jc,n)
               end if
            end if
            if (mask(i-1,j) .ne. not_covered .and.
     $          mask(i+ratiox,j) .ne. not_covered) then
               derives(ic,XDER)  = zero
               derives(ic,X2DER)  = zero
            end if
         end do

c        ::::: interpolate to fine grid
         do off = 0, ratiox - 1
            xx = (dfloat(off - ratiox/2) + half)/ratiox
            do ic = iclo, ichi
               i = ratiox*ic + off
               bdry(i,j,n) = crse(ic,jc,n)
     &              +        derives(ic,XDER) *xx
     &              +        half*derives(ic,X2DER)*xx**2
            end do
         end do
      end do

#else
      integer Norder, NN, m
      parameter (Norder = 3)
      REAL_T x(Norder), y(Norder), c(Norder), xInt
      ratiox = ratios(1)

      iclo = ARG_L1(cb)
      ichi = ARG_H1(cb)
      jc   = ARG_H2(cb)+1
      j    = hi(2)+1
      
      do n=1,nvar
         do ic=iclo,ichi
            i = ratiox*ic
            
            NN = 1
            y(NN) = crse(ic,jc,n)
            x(NN) = zero

            if (mask(i-1,j).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic-1,jc,n)
               x(NN) = -one
            else if (mask(i+2*ratiox,j).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic+2,jc,n)
               x(NN) = two
            endif
            
            if (mask(i+ratiox,j).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic+1,jc,n)
               x(NN) = one
            else if (mask(i-ratiox-1,j).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic-2,jc,n)
               x(NN) = -two
            endif
               
            if ( (mask(i-1,j).ne.not_covered).and.
     &           (mask(i+ratiox,j).ne.not_covered) ) NN = 1
            
            do off = 0,ratiox-1
               xInt = (dfloat(off - ratiox/2) + half)/ratiox
               call polyInterpCoeff(xInt, x, NN, c)
               bdry(i+off,j,n) = zero
               do m=1,NN
                  bdry(i+off,j,n) = bdry(i+off,j,n) + c(m)*y(m)
               end do
            end do
         end do
      end do
      
#endif
      return
      end

