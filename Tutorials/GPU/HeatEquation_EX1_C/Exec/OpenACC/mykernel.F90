
module my_kernel_module
    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    implicit none
contains

    !==================================
    !  OpenACC offloading subroutines
    !==================================
    subroutine init_phi (lo, hi, dat, dlo, dhi, dx, prob_lo) &
        bind(c,name="init_phi")
    integer(c_int), intent(in) :: lo(3), hi(3), dlo(3), dhi(3)
    real(amrex_real), intent(in) :: dx(3), prob_lo(3)
    real(amrex_real), intent(inout) :: dat(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))

    integer(c_int) :: i,j,k
    real(amrex_real) :: x,y,z,r2

    !$acc parallel deviceptr(dat)
    !$acc loop gang vector collapse(3) private(x, y, z, r2)
    do          k=lo(3), hi(3)
        do      j=lo(2), hi(2)
            do  i=lo(1), hi(1)
                z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
                y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
                x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
                r2 = ((x-0.25d0)**2 + (y-0.25d0)**2 + (z-0.25d0)**2) / 0.01d0
                dat(i,j,k) = 1.d0 + exp(-r2)
            end do
        end do
    end do
    !$acc end loop
    !$acc end parallel
    end subroutine init_phi


    subroutine compute_flux_x(lo, hi, fluxx, f_lo, f_hi, &
                                  phi, p_lo, p_hi, dxinv)    &
        bind(c,name="compute_flux_x")
    integer(c_int), intent(in)   :: lo(3), hi(3), f_lo(3), f_hi(3), p_lo(3), p_hi(3)
    real(amrex_real), intent(in) :: phi(p_lo(1):p_hi(1), p_lo(2):p_hi(2), p_lo(3):p_hi(3))
    real(amrex_real), intent(in), value :: dxinv
    real(amrex_real), intent(inout) :: fluxx(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3))
    
    integer(c_int) :: i,j,k
    !$acc parallel deviceptr(phi, fluxx)
    !$acc loop gang vector collapse(3)
    do          k=lo(3), hi(3)
        do      j=lo(2), hi(2)
            do  i=lo(1), hi(1)
                fluxx(i,j,k) = dxinv* ( phi(i,j,k)-phi(i-1,j,k))
            end do
        end do
    end do
    !$acc end loop
    !$acc end parallel

    end subroutine compute_flux_x


    subroutine compute_flux_y(lo, hi, fluxy, f_lo, f_hi, &
                                  phi, p_lo, p_hi, dyinv)    &
        bind(c,name="compute_flux_y")
    integer(c_int), intent(in)   :: lo(3), hi(3), f_lo(3), f_hi(3), p_lo(3), p_hi(3)
    real(amrex_real), intent(in) :: phi(p_lo(1):p_hi(1), p_lo(2):p_hi(2), p_lo(3):p_hi(3))
    real(amrex_real), intent(in), value :: dyinv
    real(amrex_real), intent(inout) :: fluxy(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3))
    
    integer(c_int) :: i,j,k
    !$acc parallel deviceptr(phi, fluxy)
    !$acc loop gang vector collapse(3)
    do          k=lo(3), hi(3)
        do      j=lo(2), hi(2)
            do  i=lo(1), hi(1)
                fluxy(i,j,k) = dyinv* ( phi(i,j,k)-phi(i,j-1,k))
            end do
        end do
    end do
    !$acc end loop
    !$acc end parallel

    end subroutine compute_flux_y


    subroutine compute_flux_z(lo, hi, fluxz, f_lo, f_hi, &
                                  phi, p_lo, p_hi, dzinv)    &
        bind(c,name="compute_flux_z")
    integer(c_int), intent(in)   :: lo(3), hi(3), f_lo(3), f_hi(3), p_lo(3), p_hi(3)
    real(amrex_real), intent(in) :: phi(p_lo(1):p_hi(1), p_lo(2):p_hi(2), p_lo(3):p_hi(3))
    real(amrex_real), intent(in), value :: dzinv
    real(amrex_real), intent(inout) :: fluxz(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3))
    
    integer(c_int) :: i,j,k
    !$acc parallel deviceptr(phi, fluxz)
    !$acc loop gang vector collapse(3)
    do          k=lo(3), hi(3)
        do      j=lo(2), hi(2)
            do  i=lo(1), hi(1)
                fluxz(i,j,k) = dzinv* ( phi(i,j,k)-phi(i,j,k-1))
            end do
        end do
    end do
    !$acc end loop
    !$acc end parallel

    end subroutine compute_flux_z


    subroutine update_phi(lo,hi,&
                              fluxx,fxlo,fxhi, &
                              fluxy,fylo,fyhi, &
                              fluxz,fzlo,fzhi, &
                              phi_old,polo,pohi, &
                              phi_new,pnlo,pnhi, &
                              dt,dxinv,dyinv,dzinv) &
        bind(c,name="update_phi")
    integer(c_int), intent(in)   :: lo(3), hi(3), fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3)
    integer(c_int), intent(in)   :: polo(3), pohi(3), pnlo(3), pnhi(3)

    real(amrex_real), intent(in) :: fluxx(fxlo(1):fxhi(1), fxlo(2):fxhi(2), fxlo(3):fxhi(3))
    real(amrex_real), intent(in) :: fluxy(fylo(1):fyhi(1), fylo(2):fyhi(2), fylo(3):fyhi(3))
    real(amrex_real), intent(in) :: fluxz(fzlo(1):fzhi(1), fzlo(2):fzhi(2), fzlo(3):fzhi(3))
    real(amrex_real), intent(in) :: phi_old(polo(1):pohi(1), polo(2):pohi(2), polo(3):pohi(3))
    real(amrex_real), intent(inout) :: phi_new(pnlo(1):pnhi(1), pnlo(2):pnhi(2), pnlo(3):pnhi(3))
    real(amrex_real), intent(in), value :: dt, dxinv, dyinv, dzinv
    
    integer(c_int) :: i,j,k
    !print *, "dt=", dt, "dxinv(:)", dxinv, dyinv, dzinv
    !$acc parallel deviceptr(fluxx, fluxy, fluxz, phi_old, phi_new)
    !$acc loop gang vector collapse(3)
    do          k=lo(3), hi(3)
        do      j=lo(2), hi(2)
            do  i=lo(1), hi(1)
                !print *, i,j,k, "old", phi_old(i,j,k), "new", phi_new(i,j,k)
                phi_new(i,j,k) = phi_old(i,j,k) &
                + dt * dxinv * (fluxx(i+1,j  ,k  ) - fluxx(i,j,k)) &
                + dt * dyinv * (fluxy(i  ,j+1,k  ) - fluxy(i,j,k)) &
                + dt * dzinv * (fluxz(i  ,j  ,k+1) - fluxz(i,j,k))
            end do
        end do
    end do
    !$acc end loop
    !$acc end parallel
    end subroutine update_phi

end module my_kernel_module
