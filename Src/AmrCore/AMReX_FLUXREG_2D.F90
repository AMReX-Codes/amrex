
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_FLUXREG_F.H"
#include "AMReX_ArrayLim.H"

#define SDIM 2

! ::: -----------------------------------------------------------
! ::: Add fine grid flux to flux register.  Flux array is a fine grid
! ::: edge based object, Register is a coarse grid edge based object.
! ::: It is assumed that the coarsened flux region contains the register
! ::: region.
! :::
! ::: INPUTS/OUTPUTS:
! ::: reg       <=> edge centered coarse grid flux register
! ::: reg_l1,reg_l2,reg_h1,reg_h2  => index limits for reg
! ::: flx        => edge centered fine grid flux array
! ::: flx_l1,flx_l2,flx_h1,flx_h2  => index limits for flx
! ::: numcomp    => number of components to update
! ::: dir        => direction normal to flux register
! ::: ratio(2)   => refinement ratios between coarse and fine
! ::: mult       => scalar multiplicative factor
! ::: -----------------------------------------------------------

    subroutine FORT_FRFINEADD(reg,reg_l1,reg_l2,reg_h1,reg_h2,flx,flx_l1,flx_l2,flx_h1,flx_h2, &
                              numcomp,dir,ratio,mult)

      implicit none

      integer    reg_l1,reg_l2,reg_h1,reg_h2
      integer    flx_l1,flx_l2,flx_h1,flx_h2
      integer    ratio(2), dir, numcomp
      REAL_T     mult
      REAL_T     reg(reg_l1:reg_h1,reg_l2:reg_h2,numcomp)
      REAL_T     flx(flx_l1:flx_h1,flx_l2:flx_h2,numcomp)

      integer    n, i, j, ic, jc, off
      integer    ratiox, ratioy

      ratiox = ratio(1)
      ratioy = ratio(2)

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
            do off = 0, ratioy-1
               do jc = ARG_L2(reg), ARG_H2(reg)
                  j = ratioy*jc + off
                  reg(ic,jc,n) = reg(ic,jc,n) + mult*flx(i,j,n)
               end do
            end do
         end do
      else
         ! flux normal to Y direction
         jc = ARG_L2(reg)
         j = jc*ratioy
         if (ARG_L2(reg) .ne. ARG_H2(reg)) then
            call bl_abort("FORT_FRFINEADD: bad register direction")
         end if
         if (j .lt. ARG_L2(flx) .or. j .gt. ARG_H2(flx)) then
            call bl_abort("FORT_FRFINEADD: index outside flux range")
         end if
         do n = 1, numcomp
            do off = 0, ratiox-1
               do ic = ARG_L1(reg), ARG_H1(reg)
                  i = ratiox*ic + off
                  reg(ic,jc,n) = reg(ic,jc,n) + mult*flx(i,j,n)
               end do
            end do
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
! ::: flx_l1,flx_l2,flx_h1,flx_h2  => index limits for flx
! ::: area       => edge centered area array
! ::: area_l1,area_l2,area_h1,area_h2 => index limits for area
! ::: numcomp    => number of components to update
! ::: dir        => direction normal to flux register
! ::: ratio(2)   => refinements ratio between coarse and fine
! ::: mult       => scalar multiplicative factor
! ::: -----------------------------------------------------------

    subroutine FORT_FRFAADD(reg,reg_l1,reg_l2,reg_h1,reg_h2,flx,flx_l1,flx_l2,flx_h1,flx_h2,area,area_l1,area_l2,area_h1,area_h2, &
                            numcomp,dir,ratio,mult)

      implicit none

      integer    reg_l1,reg_l2,reg_h1,reg_h2
      integer    flx_l1,flx_l2,flx_h1,flx_h2
      integer    area_l1,area_l2,area_h1,area_h2
      integer    ratio(2), dir, numcomp
      REAL_T     mult
      REAL_T     reg(reg_l1:reg_h1,reg_l2:reg_h2,numcomp)
      REAL_T     flx(flx_l1:flx_h1,flx_l2:flx_h2,numcomp)
      REAL_T     area(area_l1:area_h1,area_l2:area_h2)

      integer    n, i, j, ic, jc, off
      integer    ratiox, ratioy

      ratiox = ratio(1)
      ratioy = ratio(2)

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
            do off = 0, ratioy-1
               do jc = ARG_L2(reg), ARG_H2(reg)
                  j = ratioy*jc + off
                  reg(ic,jc,n) = reg(ic,jc,n) + mult*area(i,j)*flx(i,j,n)
               end do
            end do
         end do
      else
         ! flux normal to Y direction
         jc = ARG_L2(reg)
         j = jc*ratioy
         if (ARG_L2(reg) .ne. ARG_H2(reg)) then
            call bl_abort("FORT_FRFAADD: bad register direction")
         end if
         if (j .lt. ARG_L2(flx) .or. j .gt. ARG_H2(flx)) then
            call bl_abort("FORT_FRFAADD: index outside flux range")
         end if
         do n = 1, numcomp
            do off = 0, ratiox-1
               do ic = ARG_L1(reg), ARG_H1(reg)
                  i = ratiox*ic + off
                  reg(ic,jc,n) = reg(ic,jc,n) + mult*area(i,j)*flx(i,j,n)
               end do
            end do
         end do
      end if

    end subroutine FORT_FRFAADD

    subroutine FORT_FRREFLUX (lo, hi, s, slo, shi, f, flo, fhi, &
                              v, vlo, vhi, nc, mult, dir, isloface)

      implicit none

      integer, intent(in) :: lo(2), hi(2), slo(2), shi(2)
      integer, intent(in) :: flo(2), fhi(2), vlo(2), vhi(2)
      integer, intent(in) :: nc, dir, isloface
      REAL_T , intent(in) :: mult
      REAL_T , intent(inout) :: s(slo(1):shi(1),slo(2):shi(2),nc)
      REAL_T , intent(in   ) :: f(flo(1):fhi(1),flo(2):fhi(2),nc)
      REAL_T , intent(in   ) :: v(vlo(1):vhi(1),vlo(2):vhi(2))
      REAL_T bufval

      integer :: i, j, n
      if (isloface .eq. 1) then
         if (dir .eq. 0) then
            do n = 1, nc
            do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               bufval = f(i+1,j,n)/v(i,j)
               s(i,j,n) = s(i,j,n)-mult*bufval
            end do
            end do
            end do
         else
            do n = 1, nc
            do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               bufval = f(i,j+1,n)/v(i,j)
               s(i,j,n) = s(i,j,n)-mult*bufval
            end do
            end do
            end do
         end if
      else
            do n = 1, nc
            do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               bufval = f(i,j,n)/v(i,j)
               s(i,j,n) = s(i,j,n)+mult*bufval
            end do
            end do
            end do
      end if

    end subroutine FORT_FRREFLUX
