subroutine integrate_state(lo, hi, &
                           state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                           diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                           a, half_dt, min_iter, max_iter) &
                           bind(C, name="integrate_state")

!
!   Calculates the sources to be added later on.
!
!   Parameters
!   ----------
!   lo : double array (3)
!       The low corner of the current box.
!   hi : double array (3)
!       The high corner of the current box.
!   state_* : double arrays
!       The state vars
!   diag_eos* : double arrays
!       Temp and Ne
!   a : double
!       The current a
!   half_dt : double
!       time step size, in Mpc km^-1 s ~ 10^12 yr.
!
!   Returns
!   -------
!   state : double array (dims) @todo
!       The state vars
!
   
    use amrex_fort_module, only : rt => amrex_real
    use amrex_error_module, only: amrex_abort
    use meth_params_module, only : NVAR, NDIAG, heat_cool_type

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in   ) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
    real(rt), intent(inout) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    real(rt), intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, NDIAG)
    real(rt), intent(in   ) :: a, half_dt
    integer         , intent(inout) :: min_iter, max_iter

    if (heat_cool_type .eq. 1) then
        call amrex_abort("ERROR: heat_cool_type = 1 is not in function anymore.")
    else if (heat_cool_type .eq. 3) then
        call integrate_state_vode(lo, hi, state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                          diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                                  a, half_dt, min_iter, max_iter)
    else if (heat_cool_type .eq. 5) then
        call integrate_state_fcvode(lo, hi, state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                          diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                                  a, half_dt, min_iter, max_iter)
    else if (heat_cool_type .eq. 7) then
        call integrate_state_fcvode_vec(lo, hi, state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                          diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                                  a, half_dt, min_iter, max_iter)

    end if

end subroutine integrate_state
