
  subroutine amrex_eb_compute_normals ( lo, hi,           &
                                        flag, fglo, fghi, &
                                        normal, nlo, nhi, &
                                        apx, axlo, axhi,  &
                                        apy, aylo, ayhi,  &
                                        apz, azlo, azhi ) &
             bind(C, name="amrex_eb_compute_normals")

    use amrex_fort_module, only: amrex_real
    use iso_c_binding,     only: c_int

    use amrex_constants_module , only : one
    use amrex_ebcellflag_module, only : is_single_valued_cell

    implicit none

    integer         , intent(in   ) :: lo(3),hi(3),nlo(3),nhi(3),fglo(3),fghi(3)
    integer         , intent(in   ) :: axlo(3),axhi(3),aylo(3),ayhi(3),azlo(3),azhi(3)
    integer         , intent(in   ) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))
    real(amrex_real), intent(  out) :: normal(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3),3)
    real(amrex_real), intent(in   ) :: apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(amrex_real), intent(in   ) :: apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(amrex_real), intent(in   ) :: apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))

    real(amrex_real) :: axm, axp, aym, ayp, azm, azp
    real(amrex_real) :: apnorm, apnorminv, anrmx, anrmy, anrmz

    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if ( is_single_valued_cell(flag(i,j,k)) ) then
                axm = apx(i,  j  , k  )
                axp = apx(i+1,j  , k  )
                aym = apy(i,  j  , k  )
                ayp = apy(i,  j+1, k  )
                azm = apz(i,  j  , k  )
                azp = apz(i,  j  , k+1)

                apnorm = sqrt((axm-axp)**2 + (aym-ayp)**2 + (azm-azp)**2)

                apnorminv = one / apnorm
                anrmx = (axp-axm) * apnorminv   ! pointing to the wall
                anrmy = (ayp-aym) * apnorminv
                anrmz = (azp-azm) * apnorminv

                ! To fit the convention of previous mfix
                normal(i,j,k,1) = -anrmx
                normal(i,j,k,2) = -anrmy
                normal(i,j,k,3) = -anrmz
             end if

          end do
        end do
    end do

  end subroutine amrex_eb_compute_normals
