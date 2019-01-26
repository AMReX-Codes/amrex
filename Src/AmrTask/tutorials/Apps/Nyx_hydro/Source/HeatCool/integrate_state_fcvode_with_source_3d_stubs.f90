subroutine integrate_state_fcvode_with_source(lo, hi, &
                                state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                state_n ,sn_l1,sn_l2,sn_l3,sn_h1,sn_h2,sn_h3, &
                                diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                                hydro_src, src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                                reset_src,srcr_l1,srcr_l2,srcr_l3,srcr_h1,srcr_h2,srcr_h3, &
                                I_R, ir_l1, ir_l2, ir_l3, ir_h1, ir_h2, ir_h3, &
                                a, delta_time, min_iter, max_iter) &
                                bind(C, name="integrate_state_fcvode_with_source")
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
!   diag_eos_* : double arrays
!       Temp and Ne
!   hydro_src_* : doubles arrays
!       The source terms to be added to state (iterative approx.)
!   reset_src_* : doubles arrays
!       The source terms based on the reset correction
!   double array (3)
!       The low corner of the entire domain
!   a : double
!       The current a
!   delta_time : double
!       time step size, in Mpc km^-1 s ~ 10^12 yr.
!
!   Returns
!   -------
!   state : double array (dims) @todo
!       The state vars
!
    use amrex_fort_module, only : rt => amrex_real
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : NVAR, NDIAG

    implicit none

    integer         , intent(in) :: lo(3), hi(3)
    integer         , intent(in) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in) :: sn_l1, sn_l2, sn_l3, sn_h1, sn_h2, sn_h3
    integer         , intent(in) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
    integer         , intent(in) :: src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
    integer         , intent(in) :: srcr_l1, srcr_l2, srcr_l3, srcr_h1, srcr_h2, srcr_h3
    integer         , intent(in) :: ir_l1, ir_l2, ir_l3, ir_h1, ir_h2, ir_h3
    real(rt), intent(in   ) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    real(rt), intent(inout) ::  state_n(sn_l1:sn_h1, sn_l2:sn_h2,sn_l3:sn_h3, NVAR)
    real(rt), intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, NDIAG)
    real(rt), intent(in   ) :: hydro_src(src_l1:src_h1, src_l2:src_h2,src_l3:src_h3, NVAR)
    real(rt), intent(in   ) :: reset_src(srcr_l1:srcr_h1, srcr_l2:srcr_h2,srcr_l3:srcr_h3, 1)
    real(rt), intent(inout) :: I_R(ir_l1:ir_h1, ir_l2:ir_h2,ir_l3:ir_h3)
    real(rt), intent(in)    :: a, delta_time
    integer         , intent(inout) :: max_iter, min_iter

    call amrex_abort("Cannot call fcvode without compiling with USE_CVODE=TRUE")

end subroutine integrate_state_fcvode_with_source
