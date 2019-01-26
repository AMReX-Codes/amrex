subroutine integrate_state(lo, hi, &
                           state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                           diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                           a, half_dt, min_iter, max_iter) &
                           bind(C, name="integrate_state")

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NVAR, NDIAG

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in   ) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer, intent(in   ) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
    real(rt), intent(in   ) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    real(rt), intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, NDIAG)
    real(rt), intent(in   ) :: a, half_dt
    integer,  intent(inout) :: min_iter, max_iter

end subroutine integrate_state

subroutine integrate_state_with_source(lo, hi, &
                                state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                state_n ,sn_l1,sn_l2,sn_l3,sn_h1,sn_h2,sn_h3, &
                                diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                                hydro_src, src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                                reset_src,srcr_l1,srcr_l2,srcr_l3,srcr_h1,srcr_h2,srcr_h3, &
                                I_R, ir_l1, ir_l2, ir_l3, ir_h1, ir_h2, ir_h3, &
                                a, delta_time, min_iter, max_iter) &
                                bind(C, name="integrate_state_with_source")
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
    use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, &
                                   NDIAG, TEMP_COMP, NE_COMP, ZHI_COMP, gamma_minus_1

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

    integer :: i, j, k
    real(rt) :: asq,aendsq,ahalf,ahalf_inv,delta_rho,delta_e,delta_rhoe
    real(rt) :: z, z_end, a_end, rho, H_reion_z, He_reion_z
    real(rt) :: rho_orig, T_orig, ne_orig, e_orig
    real(rt) :: rho_out, T_out, ne_out, e_out
    real(rt) :: rho_src, rhoe_src, e_src
    real(rt) :: mu, mean_rhob, T_H, T_He
    real(rt) :: species(5)

end subroutine integrate_state_with_source
 
module adjust_heat_cool_module

  implicit none
 
  contains

    subroutine adjust_heat_cool(lo,hi, &
                                u_old,uo_l1,uo_l2,uo_l3,uo_h1,uo_h2,uo_h3, &
                                u_new,un_l1,un_l2,un_l3,un_h1,un_h2,un_h3, &
                                src_old, so_l1,so_l2,so_l3,so_h1,so_h2,so_h3, &
                                src_new, sn_l1,sn_l2,sn_l3,sn_h1,sn_h2,sn_h3, &
                                a_old, a_new, dt)

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NVAR

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: uo_l1, uo_l2, uo_l3, uo_h1, uo_h2, uo_h3
      integer          :: un_l1, un_l2, un_l3, un_h1, un_h2, un_h3
      integer          :: so_l1,so_l2,so_l3,so_h1,so_h2,so_h3
      integer          :: sn_l1,sn_l2,sn_l3,sn_h1,sn_h2,sn_h3
      real(rt) ::   u_old(uo_l1:uo_h1,uo_l2:uo_h2,uo_l3:uo_h3,NVAR)
      real(rt) ::   u_new(un_l1:un_h1,un_l2:un_h2,un_l3:un_h3,NVAR)
      real(rt) :: src_old(so_l1:so_h1,so_l2:so_h2,so_l3:so_h3,NVAR)
      real(rt) :: src_new(sn_l1:sn_h1,sn_l2:sn_h2,sn_l3:sn_h3,NVAR)
      real(rt) :: a_old, a_new, dt

    end subroutine adjust_heat_cool

end module adjust_heat_cool_module

! unused VODE stubs if we are not doing heating/cooling
module vode_aux_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt) :: z_vode
  integer  :: i_vode, j_vode, k_vode
  integer  :: JH_vode, JHE_vode
  logical  :: firstcall
end module vode_aux_module
