#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_BC_TYPES.H"
#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_INTERPBNDRYDATA_F.H"
#include "AMReX_ArrayLim.H"

#define SDIM 1
#define NUMDERIV 2
      
! ---------------------------------------------------------------
! ::  FORT_BDINTERPXLO : Interpolation on Xlo Face
! ::       Quadratic Interpolation from crse data
! ::       in directions transverse to face of grid
! ::
! ::  Inputs/Outputs:
! ::  bdry       <=  fine grid bndry data strip
! ::  bdry_l1,bdry_h1  => index limits of bdry
! ::  lo,hi       => index limits of grd interior
! ::  cb_l1,cb_h1    => index limits of coarsened grid interior
! ::  nvar        => number of variables to interpolate
! ::  ratios(1)   => refinement ratio
! ::  not_covered => mask is set to this value if cell is not
! ::                 covered by another fine grid and not outside the domain.
! ::  mask        => fine grid mask bndry strip
! ::  mask_l1,mask_h1  => index limits of mask array
! ::  crse        => crse grid bndry data strip
! ::  crse_l1,crse_h1  => index limits of crse array
! ::  derives     => crse grid tmp array
! ---------------------------------------------------------------

    subroutine FORT_BDINTERPXLO (bdry,bdry_l1,bdry_h1, &
                 lo,hi,cb_l1,cb_h1,nvar,ratios,not_covered, &
                 mask,mask_l1,mask_h1,crse,crse_l1,crse_h1,derives)

      implicit none

      integer  nvar, ratios(1), not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  bdry_l1,bdry_h1
      integer  mask_l1,mask_h1
      integer  crse_l1,crse_h1
      integer  cb_l1,cb_h1
      REAL_T   bdry(bdry_l1:bdry_h1,nvar)
      REAL_T   derives(1,NUMDERIV)
      integer  mask(mask_l1:mask_h1)
      REAL_T   crse(crse_l1:crse_h1,nvar)

      integer  i, ic, n

      ic   = ARG_L1(cb)-1
      i    = lo(1)-1
      
      do n=1,nvar
          bdry(i,n) = crse(ic,n)
      end do
      
    end subroutine FORT_BDINTERPXLO

! ---------------------------------------------------------------
! ::  FORT_BDINTERPXHI : Interpolation on Xhi Face
! ::       Quadratic Interpolation from crse data
! ::       in directions transverse to face of grid
! ::
! ::  Inputs/Outputs:
! ::  bdry       <=  fine grid bndry data strip
! ::  bdry_l1,bdry_h1  => index limits of bdry
! ::  lo,hi       => index limits of grd interior
! ::  cb_l1,cb_h1    => index limits of coarsened grid interior
! ::  nvar        => number of variables to interpolate
! ::  ratios(1)   => refinement ratio
! ::  not_covered => mask is set to this value if cell is not
! ::                 covered by another fine grid and not outside the domain.
! ::  mask        => fine grid mask bndry strip
! ::  mask_l1,mask_h1  => index limits of mask array
! ::  crse        => crse grid bndry data strip
! ::  crse_l1,crse_h1  => index limits of crse array
! ::  derives     => crse grid tmp array
! ---------------------------------------------------------------

    subroutine FORT_BDINTERPXHI (bdry,bdry_l1,bdry_h1, &
                 lo,hi,cb_l1,cb_h1,nvar,ratios,not_covered, &
                 mask,mask_l1,mask_h1,crse,crse_l1,crse_h1,derives)

      implicit none

      integer  nvar, ratios(1), not_covered
      integer  lo(SDIM), hi(SDIM)
      integer  bdry_l1,bdry_h1
      integer  mask_l1,mask_h1
      integer  cb_l1,cb_h1
      integer  crse_l1,crse_h1
      REAL_T   bdry(bdry_l1:bdry_h1,nvar)
      REAL_T   derives(1,NUMDERIV)
      integer  mask(mask_l1:mask_h1)
      REAL_T   crse(crse_l1:crse_h1,nvar)

      integer  i, ic, n

      ic   = ARG_H1(cb)+1
      i    = hi(1)+1

      do n=1,nvar
          bdry(i,n) = crse(ic,n)
      end do
      
    end subroutine FORT_BDINTERPXHI

