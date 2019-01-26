subroutine integrate_state_fcvode(lo, hi, &
                                  state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                  diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                                  a, half_dt, min_iter, max_iter)
!
    use amrex_error_module, only : amrex_abort
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NVAR, NDIAG

    implicit none

    integer         , intent(in) :: lo(3), hi(3)
    integer         , intent(in) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
    real(rt), intent(inout) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    real(rt), intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, NDIAG)
    real(rt), intent(in)    :: a, half_dt
    integer         , intent(inout) :: max_iter, min_iter

!   integer :: i, j, k
!   real(rt) :: z, rho
!   real(rt) :: T_orig, ne_orig, e_orig
!   real(rt) :: T_out , ne_out , e_out, mu, mean_rhob

!   logical  :: firstcall

    call amrex_abort("Cannot call fcvode without compiling with USE_CVODE=TRUE")

end subroutine integrate_state_fcvode
