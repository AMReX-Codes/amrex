
module amrex_interpbndrydata_module

  use amrex_fort_module
  use amrex_constants_module

  implicit none

contains

#define SDIM 1
#define NUMDERIV 2

! ---------------------------------------------------------------
!>  AMREX_BDINTERPXLO : Interpolation on Xlo Face
!!       Quadratic Interpolation from crse data
!!       in directions transverse to face of grid
!!
!!  \param bdry       <=  fine grid bndry data strip
!!  \param bdry_l1,bdry_h1  => index limits of bdry
!!  \param lo,hi       => index limits of grd interior
!!  \param cb_l1,cb_h1    => index limits of coarsened grid interior
!!  \param nvar        => number of variables to interpolate
!!  \param ratios(1)   => refinement ratio
!!  \param not_covered => mask is set to this value if cell is not
!!                 covered by another fine grid and not outside the domain.
!!  \param mask        => fine grid mask bndry strip
!!  \param mask_l1,mask_h1  => index limits of mask array
!!  \param crse        => crse grid bndry data strip
!!  \param crse_l1,crse_h1  => index limits of crse array
!!  \param derives     => crse grid tmp array
! ---------------------------------------------------------------

    subroutine AMREX_BDINTERPXLO (bdry,bdry_l1,bdry_h1, &
                 lo,hi,cb_l1,cb_h1,nvar,ratios,not_covered, &
                 mask,mask_l1,mask_h1,crse,crse_l1,crse_h1,derives,max_order) &
                 bind(c,name='amrex_bdinterpxlo')

      implicit none

      integer  nvar, ratios(1), not_covered, max_order
      integer  lo(SDIM), hi(SDIM)
      integer  bdry_l1,bdry_h1
      integer  mask_l1,mask_h1
      integer  crse_l1,crse_h1
      integer  cb_l1,cb_h1
      real(amrex_real)   bdry(bdry_l1:bdry_h1,nvar)
      real(amrex_real)   derives(1,NUMDERIV)
      integer  mask(mask_l1:mask_h1)
      real(amrex_real)   crse(crse_l1:crse_h1,nvar)

      integer  i, ic, n

      ic   = cb_l1-1
      i    = lo(1)-1

      do n=1,nvar
          bdry(i,n) = crse(ic,n)
      end do

    end subroutine AMREX_BDINTERPXLO

! ---------------------------------------------------------------
!>  AMREX_BDINTERPXHI : Interpolation on Xhi Face
!!       Quadratic Interpolation from crse data
!!       in directions transverse to face of grid
!!
!!  \param bdry       <=  fine grid bndry data strip
!!  \param bdry_l1,bdry_h1  => index limits of bdry
!!  \param lo,hi       => index limits of grd interior
!!  \param cb_l1,cb_h1    => index limits of coarsened grid interior
!!  \param nvar        => number of variables to interpolate
!!  \param ratios(1)   => refinement ratio
!!  \param not_covered => mask is set to this value if cell is not
!!                 covered by another fine grid and not outside the domain.
!!  \param mask        => fine grid mask bndry strip
!!  \param mask_l1,mask_h1  => index limits of mask array
!!  \param crse        => crse grid bndry data strip
!!  \param crse_l1,crse_h1  => index limits of crse array
!!  \param derives     => crse grid tmp array
! ---------------------------------------------------------------

    subroutine AMREX_BDINTERPXHI (bdry,bdry_l1,bdry_h1, &
                 lo,hi,cb_l1,cb_h1,nvar,ratios,not_covered, &
                 mask,mask_l1,mask_h1,crse,crse_l1,crse_h1,derives,max_order) &
                 bind(c,name='amrex_bdinterpxhi')

      implicit none

      integer  nvar, ratios(1), not_covered, max_order
      integer  lo(SDIM), hi(SDIM)
      integer  bdry_l1,bdry_h1
      integer  mask_l1,mask_h1
      integer  cb_l1,cb_h1
      integer  crse_l1,crse_h1
      real(amrex_real)   bdry(bdry_l1:bdry_h1,nvar)
      real(amrex_real)   derives(1,NUMDERIV)
      integer  mask(mask_l1:mask_h1)
      real(amrex_real)   crse(crse_l1:crse_h1,nvar)

      integer  i, ic, n

      ic   = cb_h1+1
      i    = hi(1)+1

      do n=1,nvar
          bdry(i,n) = crse(ic,n)
      end do

    end subroutine AMREX_BDINTERPXHI

end module amrex_interpbndrydata_module
