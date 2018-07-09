
subroutine vort_3d(lo, hi, d, nComps, dlo, dhi, vel, vort, delta) bind(c)
    use amrex_fort_module, only : amrex_real
    implicit none
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), nComps
    real(amrex_real), intent(inout) :: d(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3), 0:(nComps - 1))
    integer, intent(in) :: vel(3), vort(3)
    real(amrex_real), intent(in) :: delta(3)
    real(amrex_real) :: tdx(3)
    integer :: i, j, k
    real(amrex_real) :: uy, uz, vx, vz, wx, wy

    do i = 1, 3
        tdx(i) = REAL(1.0, amrex_real) / (2 * delta(i))
    enddo

    do k = lo(3),hi(3)
        do j = lo(2),hi(2)
            do i = lo(1),hi(1)
                uy = (d(i, j + 1, k, vel(1)) - d(i, j - 1, k, vel(1))) * tdx(2)
                uz = (d(i, j, k + 1, vel(1)) - d(i, j, k - 1, vel(1))) * tdx(3)
                vx = (d(i + 1, j, k, vel(2)) - d(i - 1, j, k, vel(2))) * tdx(1)
                vz = (d(i, j, k + 1, vel(2)) - d(i, j, k - 1, vel(2))) * tdx(3)
                wx = (d(i + 1, j, k, vel(3)) - d(i - 1, j, k, vel(3))) * tdx(1)
                wy = (d(i, j + 1, k, vel(3)) - d(i, j - 1, k, vel(3))) * tdx(2)

                d(i, j, k, vort(1)) = wy - vz
                d(i, j, k, vort(2)) = uz - wx
                d(i, j, k, vort(3)) = vx - uy
            end do
        end do
    end do
end


subroutine divu_3d(lo, hi, d, nComps, dlo, dhi, vel, divu, delta) bind(c)
    use amrex_fort_module, only : amrex_real
    implicit none
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), nComps
    real(amrex_real), intent(inout) :: d(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3), 0:(nComps - 1))
    integer, intent(in) :: vel(3), divu
    real(amrex_real), intent(in) :: delta(3)
    real(amrex_real) :: tdx(3)
    integer :: i, j, k
    real(amrex_real) :: ux, vy, wz

    do i = 1, 3
        tdx(i) = REAL(1.0, amrex_real) / (2 * delta(i))
    enddo

    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                ux = (d(i + 1, j, k, vel(1)) - d(i - 1, j, k, vel(1))) * tdx(1)
                vy = (d(i, j + 1, k, vel(2)) - d(i, j - 1, k, vel(2))) * tdx(2)
                wz = (d(i, j, k + 1, vel(3)) - d(i, j, k - 1, vel(3))) * tdx(3)

                d(i, j, k, divu) = ux + vy + wz
            end do
        end do
    end do
end


subroutine copy_3d(                                         &
    srclo, srchi, srcd, srcNComps, srcdlo, srcdhi,          &
    dstlo, dsthi, dstd, dstNComps, dstdlo, dstdhi,          &
    srccomp, dstcomp) bind(c)
    use amrex_fort_module, only : amrex_real
    implicit none
    integer, intent(in) ::                                  &
        srclo(3), srchi(3), srcdlo(3), srcdhi(3), srccomp
    integer, intent(in) ::                                  &
        dstlo(3), dsthi(3), dstdlo(3), dstdhi(3), dstcomp
    integer, intent(in) :: srcNComps, dstNComps
    real(amrex_real), intent(in) :: srcd(                   &
        srcdlo(1):srcdhi(1),                                &
        srcdlo(2):srcdhi(2),                                &
        srcdlo(3):srcdhi(3),                                &
        0:(srcNComps - 1)                                   &
        )
    real(amrex_real), intent(inout) :: dstd(                &
        dstdlo(1):dstdhi(1),                                &
        dstdlo(2):dstdhi(2),                                &
        dstdlo(3):dstdhi(3),                                &
        0:(dstNComps - 1)                                   &
        )
    integer :: i, j, k

    do k = dstdlo(3), dstdhi(3)
        do j = dstdlo(2), dstdhi(2)
            do i = dstdlo(1), dstdhi(1)
                dstd(i, j, k, dstcomp) = srcd(i, j, k, srccomp)
            end do
        end do
    end do
end

