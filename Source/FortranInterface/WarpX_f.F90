! Copyright 2019 Andrew Myers, Axel Huebl, David Grote
! Ligia Diana Amorim, Mathieu Lobet, Maxence Thevenet
! Remi Lehe, Weiqun Zhang
!
! This file is part of WarpX.
!
! License: BSD-3-Clause-LBNL


module warpx_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real, amrex_spacedim

  implicit none

contains

  subroutine warpx_build_buffer_masks (lo, hi, msk, mlo, mhi, gmsk, glo, ghi, ng) &
       bind(c, name='warpx_build_buffer_masks')
    integer, dimension(3), intent(in) :: lo, hi, mlo, mhi, glo, ghi
    integer, intent(in   ) :: gmsk(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3))
    integer, intent(inout) ::  msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    integer, intent(in) :: ng

    integer :: i,j,k

    if (amrex_spacedim .eq. 2) then

       k = lo(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (any(gmsk(i-ng:i+ng,j-ng:j+ng,k).eq.0)) then
                msk(i,j,k) = 0
             else
                msk(i,j,k) = 1
             end if
          end do
       end do

    else

       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (any(gmsk(i-ng:i+ng,j-ng:j+ng,k-ng:k+ng).eq.0)) then
                   msk(i,j,k) = 0
                else
                   msk(i,j,k) = 1
                end if
             end do
          end do
       end do

    end if

  end subroutine warpx_build_buffer_masks

end module warpx_module
