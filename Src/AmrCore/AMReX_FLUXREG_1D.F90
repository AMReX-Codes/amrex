
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_FLUXREG_F.H"
#include "AMReX_ArrayLim.H"

#define SDIM 1

! ::: -----------------------------------------------------------
! ::: Add fine grid flux to flux register.  Flux array is a fine grid
! ::: edge based object, Register is a coarse grid edge based object.
! ::: It is assumed that the coarsened flux region contains the register
! ::: region.
! :::
! ::: INPUTS/OUTPUTS:
! ::: reg       <=> edge centered coarse grid flux register
! ::: DIMS(reg)  => index limits for reg
! ::: flx        => edge centered fine grid flux array
! ::: DIMS(flx)  => index limits for flx
! ::: numcomp    => number of components to update
! ::: dir        => direction normal to flux register
! ::: ratio(2)   => refinement ratios between coarse and fine
! ::: mult       => scalar multiplicative factor
! ::: -----------------------------------------------------------

    subroutine FORT_FRFINEADD(reg,DIMS(reg),flx,DIMS(flx), &
                              numcomp,dir,ratio,mult)

      implicit none

      integer    DIMDEC(reg)
      integer    DIMDEC(flx)
      integer    ratio(1), dir, numcomp
      REAL_T     mult
      REAL_T     reg(DIMV(reg),numcomp)
      REAL_T     flx(DIMV(flx),numcomp)

      integer    n, i, ic
      integer    ratiox

      ratiox = ratio(1)

      if (dir .eq. 0) then
         ! flux normal to X direction
         ic = ARG_L1(reg)
         i = ic*ratiox
         if (ARG_L1(reg) .ne. ARG_H1(reg)) then
            call bl_abort("FORT_FRFINEADD: bad register direction")
         end if
         if (i .lt. ARG_L1(flx) .or. i .gt. ARG_H1(flx)) then
            call bl_abort("FORT_FRFINEADD: index outside flux range")
         end if
         do n = 1, numcomp
            reg(ic,n) = reg(ic,n) + mult*flx(i,n)
         end do
      end if

    end subroutine FORT_FRFINEADD

! ::: -----------------------------------------------------------
! ::: Add fine grid flux times area to flux register.
! ::: Flux array is a fine grid edge based object, Register is a
! ::: coarse grid edge based object.
! ::: It is assumed that the coarsened flux region contains the register
! ::: region.
! :::
! ::: INPUTS/OUTPUTS:
! ::: reg       <=> edge centered coarse grid flux register
! ::: rlo,rhi    => index limits for reg
! ::: flx        => edge centered fine grid flux array
! ::: DIMS(flx)  => index limits for flx
! ::: area       => edge centered area array
! ::: DIMS(area) => index limits for area
! ::: numcomp    => number of components to update
! ::: dir        => direction normal to flux register
! ::: ratio(2)   => refinements ratio between coarse and fine
! ::: mult       => scalar multiplicative factor
! ::: -----------------------------------------------------------

    subroutine FORT_FRFAADD(reg,DIMS(reg),flx,DIMS(flx),area,DIMS(area), &
                            numcomp,dir,ratio,mult)

      implicit none

      integer    DIMDEC(reg)
      integer    DIMDEC(flx)
      integer    DIMDEC(area)
      integer    ratio(1), dir, numcomp
      REAL_T     mult
      REAL_T     reg(DIMV(reg),numcomp)
      REAL_T     flx(DIMV(flx),numcomp)
      REAL_T     area(DIMV(area))

      integer    n, i, ic
      integer    ratiox

      ratiox = ratio(1)

      if (dir .eq. 0) then
         ! flux normal to X direction
         ic = ARG_L1(reg)
         i = ic*ratiox
         if (ARG_L1(reg) .ne. ARG_H1(reg)) then
            call bl_abort("FORT_FRFAADD: bad register direction")
         end if
         if (i .lt. ARG_L1(flx) .or. i .gt. ARG_H1(flx)) then
            call bl_abort("FORT_FRFAADD: index outside flux range")
         end if
         do n = 1, numcomp
            reg(ic,n) = reg(ic,n) + mult*area(i)*flx(i,n)
         end do
      end if

    end subroutine FORT_FRFAADD

    subroutine FORT_FRREFLUX (lo, hi, s, slo, shi, f, flo, fhi, &
                              v, vlo, vhi, nc, mult, dir, isloface)

      implicit none

      integer, intent(in) :: lo(1), hi(1), slo(1), shi(1)
      integer, intent(in) :: flo(1), fhi(1), vlo(1), vhi(1)
      integer, intent(in) :: nc, dir, isloface
      REAL_T , intent(in) :: mult
      REAL_T , intent(inout) :: s(slo(1):shi(1),nc)
      REAL_T , intent(in   ) :: f(flo(1):fhi(1),nc)
      REAL_T , intent(in   ) :: v(vlo(1):vhi(1))

      integer :: i, n
      if (isloface .eq. 1) then
         do n = 1, nc
            do i = lo(1), hi(1)
               s(i,n) = s(i,n)-mult*f(i+1,n)/v(i)
            end do
         end do
      else
         do n = 1, nc
            do i = lo(1), hi(1)
               s(i,n) = s(i,n)+mult*f(i,n)/v(i)
            end do
         end do
      end if

    end subroutine FORT_FRREFLUX
