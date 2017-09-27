module diff_coef_module
  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private
  public :: compute_diff_coef
contains

  subroutine compute_diff_coef (q,qlo,qhi,lambda,mu,xi,clo,chi)
    use cns_physics_module, only : Pr, C_S, T_S, cp, use_const_visc, const_visc_mu, const_visc_ki, const_lambda
    use cns_module, only : qrho,qu,qv,qw,qp,qc,qeint,qtemp,qvar
    integer, intent(in) :: qlo(3), qhi(3), clo(3), chi(3)
    real(rt), intent(in   ) :: q     (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),qvar)
    real(rt), intent(inout) :: lambda(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(rt), intent(inout) :: mu    (clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(rt), intent(inout) :: xi    (clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))

    integer :: i,j,k
    real(rt) :: Prinv

    if (use_const_visc) then
       do       k = clo(3), chi(3)
          do    j = clo(2), chi(2)
             do i = clo(1), chi(1)
                mu(i,j,k) = const_visc_mu
                xi(i,j,k) = const_visc_ki
                lambda(i,j,k) = const_lambda
             end do
          end do
       end do
    else
       Prinv = 1.d0/Pr
    
       do       k = clo(3), chi(3)
          do    j = clo(2), chi(2)
             do i = clo(1), chi(1)
                mu(i,j,k) = C_S * sqrt(q(i,j,k,qtemp)) * q(i,j,k,qtemp) / (q(i,j,k,qtemp)+T_S)
                xi(i,j,k) = 0.d0
                lambda(i,j,k) = mu(i,j,k)*cp*Prinv
             end do
          end do
       end do
    end if
  end subroutine compute_diff_coef

end module diff_coef_module
