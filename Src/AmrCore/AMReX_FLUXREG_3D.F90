
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_FLUXREG_F.H"
#include "AMReX_ArrayLim.H"

#define SDIM 3

! ::: -----------------------------------------------------------
! ::: Add fine grid flux to flux register.  Flux array is a fine grid
! ::: edge based object, Register is a coarse grid edge based object.      
! ::: It is assumed that the coarsened flux region contains the register
! ::: region.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: reg       <=> edge centered coarse grid flux register
! ::: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3  => index limits for reg
! ::: flx        => edge centered fine grid flux array
! ::: flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3  => index limits for flx
! ::: numcomp    => number of components to update
! ::: dir        => direction normal to flux register
! ::: ratio(3)   => refinement ratios between coarse and fine
! ::: mult       => scalar multiplicative factor      
! ::: -----------------------------------------------------------

    subroutine FORT_FRFINEADD(reg,reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3,flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
                              numcomp,dir,ratio,mult)

      implicit none

      integer    reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
      integer    flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
      integer    ratio(3), dir, numcomp
      REAL_T     mult
      REAL_T     reg(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3,numcomp)
      REAL_T     flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,numcomp)
      
      integer    n, i, j, k, ic, jc, kc, ioff, joff, koff
      integer    ratiox, ratioy, ratioz

      ratiox = ratio(1)
      ratioy = ratio(2)
      ratioz = ratio(3)

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

         do koff = 0, ratioz-1
            do n = 1, numcomp
               do kc = ARG_L3(reg), ARG_H3(reg)
                  k = ratioz*kc + koff
                  do joff = 0, ratioy-1            
                     do jc = ARG_L2(reg), ARG_H2(reg)
                        j = ratioy*jc + joff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + mult*flx(i,j,k,n)
                     end do
                  end do
               end do
            end do
         end do

      else if (dir .eq. 1) then
         ! flux normal to Y direction
         jc = ARG_L2(reg)
         j = jc*ratioy
         if (ARG_L2(reg) .ne. ARG_H2(reg)) then
            call bl_abort("FORT_FRFINEADD: bad register direction")
         end if
         if (j .lt. ARG_L2(flx) .or. j .gt. ARG_H2(flx)) then
            call bl_abort("FORT_FRFINEADD: index outside flux range")
         end if

         do koff = 0, ratioz-1
            do n = 1, numcomp
               do kc = ARG_L3(reg), ARG_H3(reg)
                  k = ratioz*kc + koff
                  do ioff = 0, ratiox-1            
                     do ic = ARG_L1(reg), ARG_H1(reg)
                        i = ratiox*ic + ioff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + mult*flx(i,j,k,n)
                     end do
                  end do
               end do
            end do
         end do

      else

         ! flux normal to Z direction

         kc = ARG_L3(reg)
         k = kc*ratioz
         if (ARG_L3(reg) .ne. ARG_H3(reg)) then
            call bl_abort("FORT_FRFINEADD: bad register direction")
         end if
         if (k .lt. ARG_L3(flx) .or. k .gt. ARG_H3(flx)) then
            call bl_abort("FORT_FRFINEADD: index outside flux range")
         end if

         do joff = 0, ratioy-1
            do n = 1, numcomp
               do jc = ARG_L2(reg), ARG_H2(reg)
                  j = ratioy*jc + joff
                  do ioff = 0, ratiox-1            
                     do ic = ARG_L1(reg), ARG_H1(reg)
                        i = ratiox*ic + ioff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + mult*flx(i,j,k,n)
                     end do
                  end do
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
! ::: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3  => index limits for reg
! ::: flx        => edge centered fine grid flux array
! ::: flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3  => index limits for flx
! ::: area       => edge centered area array
! ::: area_l1,area_l2,area_l3,area_h1,area_h2,area_h3 => index limits for area
! ::: numcomp    => number of components to update
! ::: dir        => direction normal to flux register
! ::: ratio(3)   => refinement ratios between coarse and fine
! ::: mult       => scalar multiplicative factor      
! ::: -----------------------------------------------------------

    subroutine FORT_FRFAADD(reg,reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3,flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3,area,area_l1,area_l2,area_l3,area_h1,area_h2,area_h3, &
                            numcomp,dir,ratio,mult)

      implicit none

      integer    reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
      integer    flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
      integer    area_l1,area_l2,area_l3,area_h1,area_h2,area_h3
      integer    ratio(3), dir, numcomp
      REAL_T     mult
      REAL_T     reg(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3,numcomp)
      REAL_T     flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,numcomp)
      REAL_T     area(area_l1:area_h1,area_l2:area_h2,area_l3:area_h3)
      
      integer    n, i, j, k, ic, jc, kc, ioff, joff, koff
      integer    ratiox, ratioy, ratioz

      ratiox = ratio(1)
      ratioy = ratio(2)
      ratioz = ratio(3)

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

         do koff = 0, ratioz-1
            do n = 1, numcomp
               do kc = ARG_L3(reg), ARG_H3(reg)
                  k = ratioz*kc + koff
                  do joff = 0, ratioy-1            
                     do jc = ARG_L2(reg), ARG_H2(reg)
                        j = ratioy*jc + joff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + &
                             mult*area(i,j,k)*flx(i,j,k,n)
                     end do
                  end do
               end do
            end do
         end do

      else if (dir .eq. 1) then

         ! flux normal to Y direction

         jc = ARG_L2(reg)
         j = jc*ratioy
         if (ARG_L2(reg) .ne. ARG_H2(reg)) then
            call bl_abort("FORT_FRFAADD: bad register direction")
         end if
         if (j .lt. ARG_L2(flx) .or. j .gt. ARG_H2(flx)) then
            call bl_abort("FORT_FRFAADD: index outside flux range")
         end if

         do koff = 0, ratioz-1
            do n = 1, numcomp
               do kc = ARG_L3(reg), ARG_H3(reg)
                  k = ratioz*kc + koff
                  do ioff = 0, ratiox-1            
                     do ic = ARG_L1(reg), ARG_H1(reg)
                        i = ratiox*ic + ioff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + &
                             mult*area(i,j,k)*flx(i,j,k,n)
                     end do
                  end do
               end do
            end do
         end do

      else

         ! flux normal to Z direction

         kc = ARG_L3(reg)
         k = kc*ratioz
         if (ARG_L3(reg) .ne. ARG_H3(reg)) then
            call bl_abort("FORT_FRFAADD: bad register direction")
         end if
         if (k .lt. ARG_L3(flx) .or. k .gt. ARG_H3(flx)) then
            call bl_abort("FORT_FRFAADD: index outside flux range")
         end if

         do joff = 0, ratioy-1
            do n = 1, numcomp
               do jc = ARG_L2(reg), ARG_H2(reg)
                  j = ratioy*jc + joff
                  do ioff = 0, ratiox-1            
                     do ic = ARG_L1(reg), ARG_H1(reg)
                        i = ratiox*ic + ioff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + &
                             mult*area(i,j,k)*flx(i,j,k,n)
                     end do
                  end do
               end do
            end do
         end do

      end if
      
    end subroutine FORT_FRFAADD

    subroutine FORT_FRREFLUX (lo, hi, s, slo, shi, f, flo, fhi, &
                              v, vlo, vhi, nc, mult, dir, isloface)

      implicit none

      integer, intent(in) :: lo(3), hi(3), slo(3), shi(3)
      integer, intent(in) :: flo(3), fhi(3), vlo(3), vhi(3)
      integer, intent(in) :: nc, dir, isloface
      REAL_T , intent(in) :: mult
      REAL_T , intent(inout) :: s(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),nc)
      REAL_T , intent(in   ) :: f(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3),nc)
      REAL_T , intent(in   ) :: v(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      integer :: i, j, k, n

      if (isloface .eq. 1) then
         if (dir .eq. 0) then
            do n = 1, nc
            do k = lo(3), hi(3)
            do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               s(i,j,k,n) = s(i,j,k,n)-mult*f(i+1,j,k,n)/v(i,j,k)
            end do
            end do
            end do
            end do
         else if (dir .eq. 1) then
            do n = 1, nc
            do k = lo(3), hi(3)
            do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               s(i,j,k,n) = s(i,j,k,n)-mult*f(i,j+1,k,n)/v(i,j,k)
            end do
            end do
            end do
            end do
         else
            do n = 1, nc
            do k = lo(3), hi(3)
            do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               s(i,j,k,n) = s(i,j,k,n)-mult*f(i,j,k+1,n)/v(i,j,k)
            end do
            end do
            end do
            end do
         end if
      else
            do n = 1, nc
            do k = lo(3), hi(3)
            do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               s(i,j,k,n) = s(i,j,k,n)+mult*f(i,j,k,n)/v(i,j,k)
            end do
            end do
            end do
            end do
      end if

    end subroutine FORT_FRREFLUX
