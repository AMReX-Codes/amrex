
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "FLUXREG_F.H"
#include "ArrayLim.H"

#define SDIM 1

c ::: -----------------------------------------------------------
c ::: Add fine grid flux to flux register.  Flux array is a fine grid
c ::: edge based object, Register is a coarse grid edge based object.
c ::: It is assumed that the coarsened flux region contains the register
c ::: region.
c :::
c ::: INPUTS/OUTPUTS:
c ::: reg       <=> edge centered coarse grid flux register
c ::: DIMS(reg)  => index limits for reg
c ::: flx        => edge centered fine grid flux array
c ::: DIMS(flx)  => index limits for flx
c ::: numcomp    => number of components to update
c ::: dir        => direction normal to flux register
c ::: ratio(2)   => refinement ratios between coarse and fine
c ::: mult       => scalar multiplicative factor
c ::: -----------------------------------------------------------

      subroutine FORT_FRFINEADD(reg,DIMS(reg),flx,DIMS(flx),
     &                          numcomp,dir,ratio,mult)

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
c        ::::: flux normal to X direction
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

      end

c ::: -----------------------------------------------------------
c ::: Add fine grid flux times area to flux register.
c ::: Flux array is a fine grid edge based object, Register is a
c ::: coarse grid edge based object.
c ::: It is assumed that the coarsened flux region contains the register
c ::: region.
c :::
c ::: INPUTS/OUTPUTS:
c ::: reg       <=> edge centered coarse grid flux register
c ::: rlo,rhi    => index limits for reg
c ::: flx        => edge centered fine grid flux array
c ::: DIMS(flx)  => index limits for flx
c ::: area       => edge centered area array
c ::: DIMS(area) => index limits for area
c ::: numcomp    => number of components to update
c ::: dir        => direction normal to flux register
c ::: ratio(2)   => refinements ratio between coarse and fine
c ::: mult       => scalar multiplicative factor
c ::: -----------------------------------------------------------

      subroutine FORT_FRFAADD(reg,DIMS(reg),flx,DIMS(flx),area,DIMS(area),
     &                        numcomp,dir,ratio,mult)

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
c        ::::: flux normal to X direction
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

      end

      subroutine FORT_FRREFLUX (lo, hi, s, slo, shi, f, flo, fhi,
     &     v, vlo, vhi, nc, mult, dir, isloface)
      implicit none
      integer, intent(in) :: lo(1), hi(1), slo(1), shi(1)
      integer, intent(in) :: flo(1), fhi(1), vlo(1), vhi(1)
      integer, intent(in) :: nc, dir, isloface
      REAL_T , intent(in) :: mult
      REAL_T , intent(inout) :: s(slo(1):shi(1),nc)
      REAL_T , intent(in   ) :: f(flo(1):fhi(1),nc)
      REAL_T , intent(in   ) :: v(vlo(1):vhi(1))
      !
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
