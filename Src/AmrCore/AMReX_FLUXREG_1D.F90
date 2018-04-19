
module amrex_fluxreg_module

  use amrex_fort_module
  use amrex_constants_module

  implicit none

contains

! ::: -----------------------------------------------------------
! ::: Add fine grid flux to flux register.  Flux array is a fine grid
! ::: edge based object, Register is a coarse grid edge based object.
! ::: It is assumed that the coarsened flux region contains the register
! ::: region.
! :::
! ::: INPUTS/OUTPUTS:
! ::: reg       <=> edge centered coarse grid flux register
! ::: reg_l1,reg_h1  => index limits for reg
! ::: flx        => edge centered fine grid flux array
! ::: flx_l1,flx_h1  => index limits for flx
! ::: numcomp    => number of components to update
! ::: dir        => direction normal to flux register
! ::: ratio(2)   => refinement ratios between coarse and fine
! ::: mult       => scalar multiplicative factor
! ::: -----------------------------------------------------------

    subroutine amrex_frfineadd(reg,reg_l1,reg_h1,flx,flx_l1,flx_h1, &
                              numcomp,dir,ratio,mult) bind(c,name='amrex_frfineadd')

      implicit none

      integer    reg_l1,reg_h1
      integer    flx_l1,flx_h1
      integer    ratio(1), dir, numcomp
      real(amrex_real)     mult
      real(amrex_real)     reg(reg_l1:reg_h1,numcomp)
      real(amrex_real)     flx(flx_l1:flx_h1,numcomp)

      integer    n, i, ic
      integer    ratiox

      ratiox = ratio(1)

      if (dir .eq. 0) then
         ! flux normal to X direction
         ic = reg_l1
         i = ic*ratiox
         if (reg_l1 .ne. reg_h1) then
            call bl_abort("amrex_frfineadd: bad register direction")
         end if
         if (i .lt. flx_l1 .or. i .gt. flx_h1) then
            call bl_abort("amrex_frfineadd: index outside flux range")
         end if
         do n = 1, numcomp
            reg(ic,n) = reg(ic,n) + mult*flx(i,n)
         end do
      end if

    end subroutine amrex_frfineadd

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
! ::: flx_l1,flx_h1  => index limits for flx
! ::: area       => edge centered area array
! ::: area_l1,area_h1 => index limits for area
! ::: numcomp    => number of components to update
! ::: dir        => direction normal to flux register
! ::: ratio(2)   => refinements ratio between coarse and fine
! ::: mult       => scalar multiplicative factor
! ::: -----------------------------------------------------------

    subroutine amrex_frfaadd(reg,reg_l1,reg_h1,flx,flx_l1,flx_h1,area,area_l1,area_h1, &
                            numcomp,dir,ratio,mult) bind(c,name='amrex_frfaadd')

      implicit none

      integer    reg_l1,reg_h1
      integer    flx_l1,flx_h1
      integer    area_l1,area_h1
      integer    ratio(1), dir, numcomp
      real(amrex_real)     mult
      real(amrex_real)     reg(reg_l1:reg_h1,numcomp)
      real(amrex_real)     flx(flx_l1:flx_h1,numcomp)
      real(amrex_real)     area(area_l1:area_h1)

      integer    n, i, ic
      integer    ratiox

      ratiox = ratio(1)

      if (dir .eq. 0) then
         ! flux normal to X direction
         ic = reg_l1
         i = ic*ratiox
         if (reg_l1 .ne. reg_h1) then
            call bl_abort("amrex_frfaadd: bad register direction")
         end if
         if (i .lt. flx_l1 .or. i .gt. flx_h1) then
            call bl_abort("amrex_frfaadd: index outside flux range")
         end if
         do n = 1, numcomp
            reg(ic,n) = reg(ic,n) + mult*area(i)*flx(i,n)
         end do
      end if

    end subroutine amrex_frfaadd

    subroutine amrex_frreflux (lo, hi, s, slo, shi, f, flo, fhi, &
                              v, vlo, vhi, nc, mult, dir, isloface) bind(c,name='amrex_frreflux')

      implicit none

      integer, intent(in) :: lo(1), hi(1), slo(1), shi(1)
      integer, intent(in) :: flo(1), fhi(1), vlo(1), vhi(1)
      integer, intent(in) :: nc, dir, isloface
      real(amrex_real) , intent(in) :: mult
      real(amrex_real) , intent(inout) :: s(slo(1):shi(1),nc)
      real(amrex_real) , intent(in   ) :: f(flo(1):fhi(1),nc)
      real(amrex_real) , intent(in   ) :: v(vlo(1):vhi(1))

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

    end subroutine amrex_frreflux

end module amrex_fluxreg_module
