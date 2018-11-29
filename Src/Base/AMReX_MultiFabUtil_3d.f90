subroutine amrex_compute_divergence (lo, hi, divu, dlo, dhi, u, ulo, uhi, &
     v, vlo, vhi, w, wlo, whi, dxinv) bind(c)
  use amrex_fort_module, only : amrex_real
  implicit none
  integer, dimension(3), intent(in) :: lo, hi, dlo, dhi, ulo, uhi, vlo, vhi, wlo, whi
  real(amrex_real), intent(inout) :: divu(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
  real(amrex_real), intent(in   ) ::    u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
  real(amrex_real), intent(in   ) ::    v(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
  real(amrex_real), intent(in   ) ::    w(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
  real(amrex_real), intent(in) :: dxinv(3)
  integer :: i,j,k
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           divu(i,j,k) = dxinv(1) * (u(i+1,j,k)-u(i,j,k)) &
                +        dxinv(2) * (v(i,j+1,k)-v(i,j,k)) &
                +        dxinv(3) * (w(i,j,k+1)-w(i,j,k))
        end do
     end do
  end do
end subroutine amrex_compute_divergence
