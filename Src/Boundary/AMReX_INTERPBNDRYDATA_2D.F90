
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_BC_TYPES.H"
#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_INTERPBNDRYDATA_F.H"

#define SDIM 2
#define NUMDERIV 2
#define XDER   1
#define X2DER  2
      
      
! ---------------------------------------------------------------
! ::  FORT_BDINTERPXLO : Interpolation on Xlo Face
! ::       Quadratic Interpolation from crse data
! ::       in directions transverse to face of grid
! ::
! ::  Inputs/Outputs:
! ::  bdry       <=  fine grid bndry data strip
! ::  bdry_l1,bdry_l2,bdry_h1,bdry_h2  => index limits of bdry
! ::  lo,hi       => index limits of grd interior
! ::  cb_l1,cb_l2,cb_h1,cb_h2    => index limits of coarsened grid interior
! ::  nvar        => number of variables to interpolate
! ::  ratios(2)   => refinement ratios
! ::  not_covered => mask is set to this value if cell is not
! ::                 covered by another fine grid and not outside the domain.
! ::  mask        => fine grid mask bndry strip
! ::  mask_l1,mask_l2,mask_h1,mask_h2  => index limits of mask array
! ::  crse        => crse grid bndry data strip
! ::  crse_l1,crse_l2,crse_h1,crse_h2  => index limits of crse array
! ::  derives     => crse grid tmp array
! ---------------------------------------------------------------

    subroutine FORT_BDINTERPXLO (bdry,bdry_l1,bdry_l2,bdry_h1,bdry_h2, &
                 lo,hi,cb_l1,cb_l2,cb_h1,cb_h2,nvar,ratios,not_covered, &
                 mask,mask_l1,mask_l2,mask_h1,mask_h2,crse,crse_l1,crse_l2,crse_h1,crse_h2,derives,max_order)

      implicit none

      integer  nvar, ratios(2), not_covered,max_order
      integer  lo(SDIM), hi(SDIM)
      integer  bdry_l1,bdry_l2,bdry_h1,bdry_h2
      integer  mask_l1,mask_l2,mask_h1,mask_h2
      integer  crse_l1,crse_l2,crse_h1,crse_h2
      integer  cb_l1,cb_l2,cb_h1,cb_h2
      REAL_T   bdry(bdry_l1:bdry_h1,bdry_l2:bdry_h2,nvar)
      REAL_T   derives(cb_l2:cb_h2,NUMDERIV)      
      integer  mask(mask_l1:mask_h1,mask_l2:mask_h2)
      REAL_T   crse(crse_l1:crse_h1,crse_l2:crse_h2,nvar)

      integer  i, j, ic, jc, off, n
      integer  jclo, jchi, ratioy

      integer Norder, NN, m
      parameter (Norder = 3)
      REAL_T x(Norder), y(Norder), c(Norder), xInt
      ratioy = ratios(2)

      jclo = cb_l2
      jchi = cb_h2
      ic   = cb_l1-1
      i    = lo(1)-1

      if (max_order.eq.1) then
         do n = 1, nvar
            do off = 0, ratioy - 1
               do jc = jclo, jchi
                  j = ratioy*jc + off
                  bdry(i,j,n) = crse(ic,jc,n)
               end do
            end do
         end do
      else

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
               
            if ( (mask(i,j-1).ne.not_covered).and. &
                 (mask(i,j+ratioy).ne.not_covered) ) NN = 1
            
            do off = 0,ratioy-1
               xInt = (dble(off - ratioy/2) + half)/ratioy
               call polyInterpCoeff(xInt, x, NN, c)
               bdry(i,j+off,n) = zero
               do m=1,NN
                  bdry(i,j+off,n) = bdry(i,j+off,n) + c(m)*y(m)
               end do
            end do
         end do
      end do
      
      endif

    end subroutine FORT_BDINTERPXLO

! ---------------------------------------------------------------
! ::  FORT_BDINTERPXHI : Interpolation on Xhi Face
! ::       Quadratic Interpolation from crse data
! ::       in directions transverse to face of grid
! ::
! ::  Inputs/Outputs:
! ::  bdry       <=  fine grid bndry data strip
! ::  bdry_l1,bdry_l2,bdry_h1,bdry_h2  => index limits of bdry
! ::  lo,hi       => index limits of grd interior
! ::  cb_l1,cb_l2,cb_h1,cb_h2    => index limits of coarsened grid interior
! ::  nvar        => number of variables to interpolate
! ::  ratios(2)   => refinement ratios
! ::  not_covered => mask is set to this value if cell is not
! ::                 covered by another fine grid and not outside the domain.
! ::  mask        => fine grid mask bndry strip
! ::  mask_l1,mask_l2,mask_h1,mask_h2  => index limits of mask array
! ::  crse        => crse grid bndry data strip
! ::  crse_l1,crse_l2,crse_h1,crse_h2  => index limits of crse array
! ::  derives     => crse grid tmp array
! ---------------------------------------------------------------

    subroutine FORT_BDINTERPXHI (bdry,bdry_l1,bdry_l2,bdry_h1,bdry_h2, &
                 lo,hi,cb_l1,cb_l2,cb_h1,cb_h2,nvar,ratios,not_covered, &
                 mask,mask_l1,mask_l2,mask_h1,mask_h2,crse,crse_l1,crse_l2,crse_h1,crse_h2,derives,max_order)

      implicit none

      integer  nvar, ratios(2), not_covered,max_order
      integer  lo(SDIM), hi(SDIM)
      integer  bdry_l1,bdry_l2,bdry_h1,bdry_h2
      integer  mask_l1,mask_l2,mask_h1,mask_h2
      integer  cb_l1,cb_l2,cb_h1,cb_h2
      integer  crse_l1,crse_l2,crse_h1,crse_h2
      REAL_T   bdry(bdry_l1:bdry_h1,bdry_l2:bdry_h2,nvar)
      REAL_T   derives(cb_l2:cb_h2,NUMDERIV)      
      integer  mask(mask_l1:mask_h1,mask_l2:mask_h2)
      REAL_T   crse(crse_l1:crse_h1,crse_l2:crse_h2,nvar)

      integer  i, j, ic, jc, off, n
      integer  jclo, jchi, ratioy

      integer Norder, NN, m
      parameter (Norder = 3)
      REAL_T x(Norder), y(Norder), c(Norder), xInt

      ratioy = ratios(2)

      jclo = cb_l2
      jchi = cb_h2
      ic   = cb_h1+1
      i    = hi(1)+1

      if (max_order.eq.1) then
         do n = 1, nvar
            do off = 0, ratioy - 1
               do jc = jclo, jchi
                  j = ratioy*jc + off
                  bdry(i,j,n) = crse(ic,jc,n)
               end do
            end do
         end do
      else
      
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
               
            if ( (mask(i,j-1).ne.not_covered).and. &
                 (mask(i,j+ratioy).ne.not_covered) ) NN = 1
            
            do off = 0,ratioy-1
               xInt = (dble(off - ratioy/2) + half)/ratioy
               call polyInterpCoeff(xInt, x, NN, c)
               bdry(i,j+off,n) = zero
               do m=1,NN
                  bdry(i,j+off,n) = bdry(i,j+off,n) + c(m)*y(m)
               end do
            end do
         end do
      end do
      
      endif

    end subroutine FORT_BDINTERPXHI

! ---------------------------------------------------------------
! ::  FORT_BDINTERPYLO : Interpolation on Ylo Face
! ::       Quadratic Interpolation from crse data
! ::       in directions transverse to face of grid
! ::
! ::  Inputs/Outputs:
! ::  bdry       <=  fine grid bndry data strip
! ::  bdry_l1,bdry_l2,bdry_h1,bdry_h2  => index limits of bdry
! ::  lo,hi       => index limits of grd interior
! ::  cb_l1,cb_l2,cb_h1,cb_h2    => index limits of coarsened grid interior
! ::  nvar        => number of variables to interpolate
! ::  ratios(2)   => refinement ratios
! ::  not_covered => mask is set to this value if cell is not
! ::                 covered by another fine grid and not outside the domain.
! ::  mask        => fine grid mask bndry strip
! ::  mask_l1,mask_l2,mask_h1,mask_h2  => index limits of mask array
! ::  crse        => crse grid bndry data strip
! ::  crse_l1,crse_l2,crse_h1,crse_h2  => index limits of crse array
! ::  derives     => crse grid tmp array
! ---------------------------------------------------------------

    subroutine FORT_BDINTERPYLO (bdry,bdry_l1,bdry_l2,bdry_h1,bdry_h2, &
                 lo,hi,cb_l1,cb_l2,cb_h1,cb_h2,nvar,ratios,not_covered, &
                 mask,mask_l1,mask_l2,mask_h1,mask_h2,crse,crse_l1,crse_l2,crse_h1,crse_h2,derives,max_order)

      implicit none

      integer  nvar, ratios(2), not_covered,max_order
      integer  lo(SDIM), hi(SDIM)
      integer  bdry_l1,bdry_l2,bdry_h1,bdry_h2
      integer  mask_l1,mask_l2,mask_h1,mask_h2
      integer  cb_l1,cb_l2,cb_h1,cb_h2
      integer  crse_l1,crse_l2,crse_h1,crse_h2
      REAL_T   bdry(bdry_l1:bdry_h1,bdry_l2:bdry_h2,nvar)
      REAL_T   derives(cb_l1:cb_h1,NUMDERIV)
      integer  mask(mask_l1:mask_h1,mask_l2:mask_h2)
      REAL_T   crse(crse_l1:crse_h1,crse_l2:crse_h2,nvar)

      integer  i, j, ic, jc, off, n
      integer  iclo, ichi, ratiox

      integer Norder, NN, m
      parameter (Norder = 3)
      REAL_T x(Norder), y(Norder), c(Norder), xInt

      ratiox = ratios(1)

      iclo = cb_l1
      ichi = cb_h1
      jc   = cb_l2-1
      j    = lo(2)-1

      if (max_order.eq.1) then
      do n = 1, nvar
         do off = 0, ratiox - 1
            do ic = iclo, ichi
               i = ratiox*ic + off
               bdry(i,j,n) = crse(ic,jc,n)
            end do
         end do
      end do
      else
      
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
               
            if ( (mask(i-1,j).ne.not_covered).and. &
                 (mask(i+ratiox,j).ne.not_covered) ) NN = 1
            
            do off = 0,ratiox-1
               xInt = (dble(off - ratiox/2) + half)/ratiox
               call polyInterpCoeff(xInt, x, NN, c)
               bdry(i+off,j,n) = zero
               do m=1,NN
                  bdry(i+off,j,n) = bdry(i+off,j,n) + c(m)*y(m)
               end do
            end do
         end do
      end do
      
      endif

    end subroutine FORT_BDINTERPYLO

! ---------------------------------------------------------------
! ::  FORT_BDINTERPYHI : Interpolation on Yhi Face
! ::       Quadratic Interpolation from crse data
! ::       in directions transverse to face of grid
! ::
! ::  Inputs/Outputs:
! ::  bdry       <=  fine grid bndry data strip
! ::  bdry_l1,bdry_l2,bdry_h1,bdry_h2  => index limits of bdry
! ::  lo,hi       => index limits of grd interior
! ::  cb_l1,cb_l2,cb_h1,cb_h2    => index limits of coarsened grid interior
! ::  nvar        => number of variables to interpolate
! ::  ratios(2)   => refinement ratios
! ::  not_covered => mask is set to this value if cell is not
! ::                 covered by another fine grid and not outside the domain.
! ::  mask        => fine grid mask bndry strip
! ::  mask_l1,mask_l2,mask_h1,mask_h2  => index limits of mask array
! ::  crse        => crse grid bndry data strip
! ::  crse_l1,crse_l2,crse_h1,crse_h2  => index limits of crse array
! ::  derives     => crse grid tmp array
! ---------------------------------------------------------------

    subroutine FORT_BDINTERPYHI (bdry,bdry_l1,bdry_l2,bdry_h1,bdry_h2, &
                 lo,hi,cb_l1,cb_l2,cb_h1,cb_h2,nvar,ratios,not_covered, &
                 mask,mask_l1,mask_l2,mask_h1,mask_h2,crse,crse_l1,crse_l2,crse_h1,crse_h2,derives,max_order)

      implicit none

      integer  nvar, ratios(2), not_covered,max_order
      integer  lo(SDIM), hi(SDIM)
      integer  bdry_l1,bdry_l2,bdry_h1,bdry_h2
      integer  mask_l1,mask_l2,mask_h1,mask_h2
      integer  cb_l1,cb_l2,cb_h1,cb_h2
      integer  crse_l1,crse_l2,crse_h1,crse_h2
      REAL_T   bdry(bdry_l1:bdry_h1,bdry_l2:bdry_h2,nvar)
      REAL_T   derives(cb_l1:cb_h1,NUMDERIV)
      integer  mask(mask_l1:mask_h1,mask_l2:mask_h2)
      REAL_T   crse(crse_l1:crse_h1,crse_l2:crse_h2,nvar)

      integer  i, j, ic, jc, off, n
      integer  iclo, ichi, ratiox

      integer Norder, NN, m
      parameter (Norder = 3)
      REAL_T x(Norder), y(Norder), c(Norder), xInt
      
      ratiox = ratios(1)

      iclo = cb_l1
      ichi = cb_h1
      jc   = cb_h2+1
      j    = hi(2)+1

      if (max_order.eq.1) then
         do n = 1, nvar
            do off = 0, ratiox - 1
               do ic = iclo, ichi
                  i = ratiox*ic + off
                  bdry(i,j,n) = crse(ic,jc,n)
               end do
            end do
         end do
      else

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
               
            if ( (mask(i-1,j).ne.not_covered).and. &
                 (mask(i+ratiox,j).ne.not_covered) ) NN = 1
            
            do off = 0,ratiox-1
               xInt = (dble(off - ratiox/2) + half)/ratiox
               call polyInterpCoeff(xInt, x, NN, c)
               bdry(i+off,j,n) = zero
               do m=1,NN
                  bdry(i+off,j,n) = bdry(i+off,j,n) + c(m)*y(m)
               end do
            end do
         end do
      end do
      
      endif

    end subroutine FORT_BDINTERPYHI
