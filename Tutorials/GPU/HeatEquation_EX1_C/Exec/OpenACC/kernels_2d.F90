
module my_kernel_acc_module
    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    implicit none
contains

    !==================================
    !  OpenACC offloading subroutines
    !==================================
    subroutine init_phi_acc (lo, hi, dat, dlo, dhi, dx, prob_lo) &
        bind(c,name="init_phi_acc")
    integer(c_int), intent(in) :: lo(2), hi(2), dlo(2), dhi(2)
    real(amrex_real), intent(in) :: dx(2), prob_lo(2)
    real(amrex_real), intent(inout) :: dat(dlo(1):dhi(1),dlo(2):dhi(2))

    integer(c_int) :: i,j
    real(amrex_real) :: x,y,r2

    !$acc parallel deviceptr(dat)
    !$acc loop gang vector collapse(2) private(x, y, r2)
    do      j=lo(2), hi(2)
        do  i=lo(1), hi(1)
            y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
            r2 = ((x-0.25d0)**2 + (y-0.25d0)**2) / 0.01d0
            dat(i,j) = 1.d0 + exp(-r2)
        end do
    end do
    !$acc end loop
    !$acc end parallel
    end subroutine init_phi_acc


    subroutine compute_flux_x_acc(lo, hi, fluxx, f_lo, f_hi, &
                                  phi, p_lo, p_hi, dxinv)    &
        bind(c,name="compute_flux_x_acc")
    integer(c_int), intent(in)   :: lo(2), hi(2), f_lo(2), f_hi(2), p_lo(2), p_hi(2)
    real(amrex_real), intent(in) :: phi(p_lo(1):p_hi(1), p_lo(2):p_hi(2))
    real(amrex_real), intent(in), value :: dxinv
    real(amrex_real), intent(inout) :: fluxx(f_lo(1):f_hi(1), f_lo(2):f_hi(2))
    integer(c_int) :: i,j
    !$acc parallel deviceptr(phi, fluxx)
    !$acc loop gang vector collapse(2)
        do      j=lo(2), hi(2)
            do  i=lo(1), hi(1)
                fluxx(i,j) = dxinv* ( phi(i,j)-phi(i-1,j))
            end do
        end do
    !$acc end loop
    !$acc end parallel

    end subroutine compute_flux_x_acc


    subroutine compute_flux_y_acc(lo, hi, fluxy, f_lo, f_hi, &
                                  phi, p_lo, p_hi, dyinv)    &
        bind(c,name="compute_flux_y_acc")
    integer(c_int), intent(in)   :: lo(2), hi(2), f_lo(2), f_hi(2), p_lo(2), p_hi(2)
    real(amrex_real), intent(in) :: phi(p_lo(1):p_hi(1), p_lo(2):p_hi(2))
    real(amrex_real), intent(in), value :: dyinv
    real(amrex_real), intent(inout) :: fluxy(f_lo(1):f_hi(1), f_lo(2):f_hi(2))
    
    integer(c_int) :: i,j
    !$acc parallel deviceptr(phi, fluxy)
    !$acc loop gang vector collapse(2)
    do      j=lo(2), hi(2)
        do  i=lo(1), hi(1)
            fluxy(i,j) = dyinv* ( phi(i,j) - phi(i,j-1))
        end do
    end do
    !$acc end loop
    !$acc end parallel

    end subroutine compute_flux_y_acc


    subroutine update_phi_acc(lo,hi,&
                              fluxx,fxlo,fxhi, &
                              fluxy,fylo,fyhi, &
                              phi_old,polo,pohi, &
                              phi_new,pnlo,pnhi, &
                              dt,dxinv,dyinv) &
        bind(c,name="update_phi_acc")
    integer(c_int), intent(in)      :: lo(2), hi(2), polo(2), pohi(2), pnlo(2), pnhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
    real(amrex_real), intent(in)    :: phi_old(polo(1):pohi(1),polo(2):pohi(2))
    real(amrex_real), intent(inout) :: phi_new(pnlo(1):pnhi(1),pnlo(2):pnhi(2))
    real(amrex_real), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    real(amrex_real), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2))
    real(amrex_real), intent(in), value :: dt, dxinv, dyinv

    ! local variables
    integer(c_int) :: i,j
    !$acc parallel deviceptr(fluxx, fluxy, phi_old, phi_new)
    !$acc loop gang vector collapse(2)
    do      j=lo(2), hi(2)
        do  i=lo(1), hi(1)
            phi_new(i,j) = phi_old(i,j) &
            + dt * dxinv * (fluxx(i+1,j) - fluxx(i,j)) &
            + dt * dyinv * (fluxy(i  ,j+1) - fluxy(i,j))
        end do
    end do
    !$acc end loop
    !$acc end parallel
    end subroutine update_phi_acc

end module my_kernel_acc_module
