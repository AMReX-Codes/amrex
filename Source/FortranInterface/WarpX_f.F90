
module warpx_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real, amrex_spacedim

  implicit none

contains

  subroutine warpx_compute_E (lo, hi, &
       phi, phlo, phhi, &
       Ex,  Exlo, Exhi, &
       Ey,  Eylo, Eyhi, &
       Ez,  Ezlo, Ezhi, &
       dx) bind(c,name='warpx_compute_E')
    integer(c_int), intent(in) :: lo(3), hi(3), phlo(3), phhi(3), Exlo(3), Exhi(3),  &
         Eylo(3), Eyhi(3), Ezlo(3), Ezhi(3)
    real(amrex_real), intent(in)  :: dx(3)
    real(amrex_real), intent(in   ) :: phi(phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))
    real(amrex_real), intent(inout) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),Exlo(3):Exhi(3))
    real(amrex_real), intent(inout) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),Eylo(3):Eyhi(3))
    real(amrex_real), intent(inout) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),Ezlo(3):Ezhi(3))

    integer :: i, j, k
    real(amrex_real) :: dxinv(3)

    dxinv = 1.0 / dx

    do k    = lo(3), hi(3)
       do j = lo(2), hi(2)

          do i = lo(1), hi(1)-1
             Ex(i,j,k) = dxinv(1) * (phi(i,j,k) - phi(i+1,j,k))
          end do

          if (j < hi(2)) then
             do i = lo(1), hi(1)
                Ey(i,j,k) = dxinv(2) * (phi(i,j,k) - phi(i,j+1,k))
             end do
          end if

          if (k < hi(3)) then
             do i = lo(1), hi(1)
                Ez(i,j,k) = dxinv(3) * (phi(i,j,k) - phi(i,j,k+1))
             end do
          end if

       end do
    end do

  end subroutine warpx_compute_E

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
